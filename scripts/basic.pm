#!/usr/bin/perl -w
use strict;
use POSIX qw(strftime);
use Data::Dumper;

sub getTime {
    my ( $sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst ) = localtime();
    my $time = strftime( "%Y-%m-%d %H:%M:%S", localtime() );
    return $time;
}

sub checkPrograms {
    my @tools = @_;
    foreach my $tool (@tools) {
        my $tool_path = `which $tool`;
        print "\t$tool: $tool_path";
    }
}

sub run {
    my $cmd = shift;
    if ( !system($cmd) ) {
        print "Succeeded: $cmd\n";
        return 1;
    }
    else {
        print "Failed: $cmd\n";
        return 0;
    }
}

sub calculate_ccs_pass {
    my $ccs = shift;
    open CCS, "samtools view $ccs | " or die "Could not open file $ccs:$!\n";
    open OUT, ">ccs.pass.txt"         or die "Could not create file ccs.pass.txt:$!\n";
    while (<CCS>) {
        my @fields = split;
        my $id     = $fields[0];

        my $pass = 0;
        if (/\s+np:i:(\d+)\s+/) {
            $pass = $1;
        }
        print OUT "$id\t$pass\n" if $pass > 0;
    }
    close CCS;
    close OUT;
}

sub call_polyA_tail {
    my ( $sample, $minimap2_bam, $ccs2gene_file, $ccs_pass_file, $out, $min_pal ) = @_;

    my %ccs_pass = ();
    open PASS, $ccs_pass_file or die "Could not open file $ccs_pass_file $!\n";
    while (<PASS>) {
        chomp;
        my ( $ccs, $pass ) = split /\s+/;
        $ccs_pass{$ccs} = $pass;
    }
    close PASS;

    my %ccs2gene = ();
    open CCS2GENE, $ccs2gene_file or die "Could not open file $ccs2gene_file $!\n";
    while (<CCS2GENE>) {
        chomp;
        my ( $ccs, $gene ) = split /\s+/;
        $ccs2gene{$ccs} = $gene;
    }
    close CCS2GENE;

    open BAM, "samtools view $minimap2_bam | " or die "Could not open file $minimap2_bam $!\n";
    open OUT, ">$out" or die "Could not open file $out !\n";
    print OUT "sample\tid\tgene\tpass\tis_polyA\tA\tT\tC\tG\tnon_A\tpolyA_seq\tccs_seq\n";
    while (<BAM>) {
        chomp;
        my @fields = split /\t/;
        my ( $id, $flag, $cigar, $seq ) = @fields[ 0, 1, 5, 9 ];

        my $clip  = 0;
        my $polyA = "";
        my ( $A, $T, $C, $G, $non_A ) = ( 0, 0, 0, 0, 0 );
        if ( $flag == 0 ) {    # read forward strand
            if ( $cigar =~ /(\d+)S$/ ) {
                $clip = $1;
                next if $clip < $min_pal;
            }
            else {
                next;
            }

            $seq   = reverse($seq);
            $polyA = substr( $seq, 0, $clip );
            $A     = ( $polyA =~ tr/A/A/ );
            $T     = ( $polyA =~ tr/T/T/ );
            $C     = ( $polyA =~ tr/C/C/ );
            $G     = ( $polyA =~ tr/G/G/ );

        }
        elsif ( $flag == 16 ) {    # read reverse strand
            if ( $cigar =~ /^(\d+)S/ ) {
                $clip = $1;
                next if $clip < $min_pal;
            }
            else {
                next;
            }

            $polyA = substr( $seq, 0, $clip );
            my @polyA = split( //, $polyA );
            $A = ( $polyA =~ tr/T/T/ );
            $T = ( $polyA =~ tr/A/A/ );
            $C = ( $polyA =~ tr/G/G/ );
            $G = ( $polyA =~ tr/C/C/ );
            $seq =~ tr/ATCGNatcgn/TAGCNtagcn/;
            $polyA =~ tr/ATCGNatcgn/TAGCNtagcn/;
        }

        next if $clip < $min_pal;

        $non_A = $T + $C + $G;
        my $gene = "";
        if ( exists $ccs2gene{$id} ) {
            $gene = $ccs2gene{$id};
        }
        else {
            $gene = "Unknown";
        }
        my $pass = $ccs_pass{$id};

        my $is_polyA = 1;
        $is_polyA = 0 if $clip < $min_pal;
        $is_polyA = 0 if !( $polyA =~ /A{5}/ );
        #$is_polyA = 0 if $T / $clip >= 0.25 || $C / $clip >= 0.25 || $G / $clip >= 0.25;
        $is_polyA = 0 if $non_A / $clip >= 0.5;
		$is_polyA = 0 if $non_A >= 20;

        # seq orientation: terminal => 3UTR
        print OUT "$sample\t$id\t$gene\t$pass\t$is_polyA\t$A\t$T\t$C\t$G\t$non_A\t$polyA\t$seq\n";
    }
    close BAM;
    close OUT;
}

sub call_nonA_in_polyA_tail {
    my ( $sample, $polyA_tail_length, $min_pass, $min_pal, $limit, $type ) = @_;

    # type = 1, terminal -> 3'UTR (trunc 100 nt)
    # type = 2, 3UTR -> terminal (trunc 100 nt)
    # type = 3, terminal -> 3'UTR (scale)
    my @total       = ();
    my @A           = ();
    my @T           = ();
    my @C           = ();
    my @G           = ();
    my $bins        = 100;
    my $window_size = 1 / $bins;

    if ( $type == 1 || $type == 2 ) {
        @total = (0) x $limit;
        @A     = (0) x $limit;
        @T     = (0) x $limit;
        @C     = (0) x $limit;
        @G     = (0) x $limit;
    }
    elsif ( $type == 3 ) {
        @total = (0) x $bins;
        @A     = (0) x $bins;
        @T     = (0) x $bins;
        @C     = (0) x $bins;
        @G     = (0) x $bins;
    }

    open PA, $polyA_tail_length or die "Could not open file $polyA_tail_length $!\n";
    my $head = <PA>;
    while (<PA>) {
        chomp;
        my ( $sample_name, $id, $gene, $pass, $is_polyA, $A, $T, $C, $G, $non_A, $polyA_seq, $ccs_seq ) = split /\t/;
        next if !$is_polyA;
        next if $pass < $min_pass;
        my $polyA_len = $A + $T + $C + $G;
        next if $polyA_len < $min_pal;

        my $len = $limit;
        $len = $polyA_len if $polyA_len < $limit;

        $polyA_seq = reverse($polyA_seq) if $type == 1;
        my @polyA_seq = split( //, $polyA_seq );
        if ( $type == 1 || $type == 2 ) {

            # type = 1: 3UTR(100nt) => terminal
            # type = 2; terminal(100nt) => 3UTR
            for ( my $i = 0 ; $i < $len ; $i++ ) {
                $total[$i]++;
                if ( $polyA_seq[$i] eq "A" ) {
                    $A[$i]++;
                }
                elsif ( $polyA_seq[$i] eq "T" ) {
                    $T[$i]++;
                }
                elsif ( $polyA_seq[$i] eq "C" ) {
                    $C[$i]++;
                }
                elsif ( $polyA_seq[$i] eq "G" ) {
                    $G[$i]++;
                }
            }
        }
        elsif ( $type == 3 ) {    # orientation: 3UTR => terminal, scaled
            $polyA_seq = reverse($polyA_seq);
			@polyA_seq = split(//, $polyA_seq);
            for ( my $i = 0 ; $i < $polyA_len ; $i++ ) {
                my $relative_pos = ( $i + 1 ) / $polyA_len;
                for ( my $j = 0 ; $j < $bins ; $j++ ) {
                    if ( $relative_pos >= $j * $window_size && $relative_pos < ( $j + 1 ) * $window_size ) {
                        $total[$j]++;
                        if ( $polyA_seq[$i] eq "A" ) {
                            $A[$j]++;
                        }
                        elsif ( $polyA_seq[$i] eq "T" ) {
                            $T[$j]++;
                        }
                        elsif ( $polyA_seq[$i] eq "C" ) {
                            $C[$j]++;
                        }
                        elsif ( $polyA_seq[$i] eq "G" ) {
                            $G[$j]++;
                        }
                    }
                }
            }
        }
    }
    if ( $type == 1 ) {
        open OUT, ">$sample.3UTR_nonA.txt" or die "Could not open file $sample.3UTR_nonA.txt $!\n";
        print "Output file: $sample.3UTR_nonA.txt\n";
    }
    elsif ( $type == 2 ) {
        open OUT, ">$sample.terminal_nonA.txt" or die "Could not open file $sample.terminal_nonA.txt $!\n";
        print "Output file: $sample.terminal_nonA.txt\n";
    }
    elsif ( $type == 3 ) {
        open OUT, ">$sample.scaled_nonA.txt" or die "Could not open file $sample.scaled_nonA.txt $!\n";
        print "Output file: $sample.scaled_nonA.txt\n";
    }

    if ( $type == 1 || $type == 2 ) {
        print OUT "Pos\tA\tT\tC\tG\n";
        for ( my $i = 0 ; $i < $limit ; $i++ ) {
            my $A_freq = ( $A[$i] / $total[$i] ) * 100;
            my $T_freq = ( $T[$i] / $total[$i] ) * 100;
            my $C_freq = ( $C[$i] / $total[$i] ) * 100;
            my $G_freq = ( $G[$i] / $total[$i] ) * 100;
            my $pos    = $i + 1;
            print OUT "$pos\t$A_freq\t$T_freq\t$C_freq\t$G_freq\n";
        }
    }
    elsif ( $type == 3 ) {
        print OUT "Pos\tA\tT\tC\tG\n";
        for ( my $i = 0 ; $i < $bins ; $i++ ) {
            my $A_freq = ( $A[$i] / $total[$i] ) * 100;
            my $T_freq = ( $T[$i] / $total[$i] ) * 100;
            my $C_freq = ( $C[$i] / $total[$i] ) * 100;
            my $G_freq = ( $G[$i] / $total[$i] ) * 100;
            my $pos    = $i + 1;
            print OUT "$pos\t$A_freq\t$T_freq\t$C_freq\t$G_freq\n";
        }
    }
    close PA;
    close OUT;
}

1;
