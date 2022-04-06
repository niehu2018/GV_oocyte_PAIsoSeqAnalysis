#!/usr/bin/env Rscript
library("ggplot2")
library("plyr")

usage <- "
Description:
	This script is used to plot the distribution of polyA tail length
Usage:
	Rscript plot_polyA_tail_length.R min_pass min_ccs_per_gene level type
Examples:
	Rscript sample.polyA_tail_length.txt output_prefix 1 3 transcript density
	Rscript sample.polyA_tail_length.txt output_prefix 1 3 transcript violin
	Rscript sample.polyA_tail_length.txt output_prefix 1 3 gene density
	Rscript sample.polyA_tail_length.txt output_prefix 1 3 gene violin
Authors:
	Hu Nie, niehu@genetics.ac.cn
Versions:
	Version 0.0.1, 2019-04-27
Notes:
 The min_ccs_per_gene option has no effect on transcript level plot
"

Args <- commandArgs()
if(length(Args) != 11){
	cat(usage);
	quit()
}

# Read command line arguments
input              <- Args[6]
output_prefix      <- Args[7]
min_pass           <- as.integer(Args[8])
min_ccs_per_gene   <- as.integer(Args[9])
level              <- as.integer(Args[10])
type               <- Args[11]

data <- read.table(input , header = TRUE, stringsAsFactors = FALSE)
data <- data[data$is_polyA == 1, ]            # select transcripts with polyA tail
data <- data[data$pass >= min_pass, ]         # select transcripts with CCS pass >= min_pass
#data <- data[data$gene != "Unknown", ]        # select transcripts assigned to genes
data$length <- data$A + data$T + data$C + data$G

# Calculate the median length of polyA tail
sample <- data$sample[1]  # add sample
gene_median_pal <- ddply(data, "gene", summarise, median=median(length))
gene_median_pal$sample <- sample

ccs_per_gene <- count(data$gene)
colnames(ccs_per_gene) <- c("gene","count")
selected_gene <- ccs_per_gene[ccs_per_gene$count >= min_ccs_per_gene, ]$gene
gene_median_pal <- gene_median_pal[gene_median_pal$gene %in% selected_gene, ]
gene_median_pal <- gene_median_pal[gene_median_pal$gene != "Unknown", ]

gene_median_pal.bk <- gene_median_pal

if(length(summary(gene_median_pal$median > 200)) == 3 ){
  gene_median_pal[gene_median_pal$median > 200, ]$median = 200
}

data.bk <- data
if(length(summary(data$length > 200)) == 3){
  data[data$length > 200, ]$length = 200
}
long <- data.bk[ data.bk$length > 200 , ]

# Plot
if(type == "density" && level == "transcript"){
  out <- paste(output_prefix,"transcript.polyA_tail_length.density_plot.pdf", sep =".")
  ggplot(data = data, aes(x=length, col = sample)) + 
    geom_density() + 
    geom_vline(aes(xintercept=median(data.bk$length)), linetype="dashed", size=0.25) + 
    xlab("Poly(A) tail length (nt)") + 
    ylab("Density") + 
    scale_x_continuous(limits = c(0,200), breaks = seq(0,200, by = 20)) + 
    ggtitle( paste("Distribution of polyA tail length\n Median = ", median(data.bk$length)," nt",sep =" "))
  ggsave(out, width = 8, height = 5)
}

if(type == "violin" && level == "transcript"){
  out <- paste(output_prefix,"transcript.polyA_tail_length.violin_plot.pdf", sep =".")
  ggplot(data = data, aes(x=sample, y= length, col =sample)) + 
  	geom_violin()  +  
	stat_summary(data = data.bk, fun.y=median, geom="point", size=2, color="red") + 
    geom_hline(aes(yintercept=median(data.bk$length)), linetype="dashed", size=0.25) + 
    ylab("Poly(A) tail length (nt)") + 
	xlab("") + 
    scale_y_continuous(limits = c(0,200), breaks = seq(0,200, by = 20)) +
	ggtitle( paste("Distribution of polyA tail length\n Median = ", median(data.bk$length), " nt",sep =" "))
  ggsave(out, width=4, height=5)
}

if(type == "density" && level == "gene"){
	out <- paste(output_prefix,"gene.polyA_tail_length.density_plot.pdf", sep = ".")
	ggplot(data = gene_median_pal , aes(x = median, col = sample)) +
		  geom_density() +
          geom_vline(aes(xintercept=median(gene_median_pal.bk$median)), linetype="dashed", size=0.25) + 
		  xlab("Poly(A) tail length (nt)") +
		  ylab("Density") + 
	      scale_x_continuous(limits = c(0,200), breaks = seq(0,200, by = 20)) + 
		  ggtitle( paste("Distribution of polyA tail length\n Median = ", median(gene_median_pal.bk$median), " nt",sep =" "))
  	ggsave(out, width=8, height=5)
}

if(type == "violin" && level == "gene"){
  out <- paste(output_prefix,"gene.polyA_tail_length.violin_plot.pdf", sep = ".")
  ggplot(data = gene_median_pal , aes(x=sample, y= median, col = sample)) +
		      geom_violin() +
			  stat_summary(data = gene_median_pal.bk, fun.y=median, geom="point", size=2, color="red") +
			  geom_hline(aes(yintercept=median(gene_median_pal.bk$median)), linetype="dashed", size=0.25) + 
		      ylab("Poly(A) tail length (nt)") +
		      xlab("") +
			  scale_y_continuous(limits = c(0,200), breaks = seq(0,200, by = 20)) +
		      ggtitle( paste("Distribution of polyA tail length\n Median = ", median(gene_median_pal.bk$median), " nt",sep =" "))
  	ggsave(out, width=4, height=5)
}

# long
long_out <- paste(output_prefix, "long_polyA_tail_length_distribution.pdf", sep =".")

# breaks
breaks <- seq(200,300, by=20)
breaks <- c(breaks,350,400, 450, 500, max(long$length))

# breaks name
breaks.name=c()
for(i in 1:(length(breaks)-2)){
  breaks.name[i] <- paste(breaks[i],breaks[i+1], sep ="-")
}
breaks.name[1] <- "201-220"
breaks.name[length(breaks)-1] <- ">500"

# hist
hist <- hist(long$length, breaks = breaks, plot=FALSE)
long_polyA_length_bins <- cbind(breaks.name, hist$counts)
colnames(long_polyA_length_bins) <- c("range","counts")

# write table
table_out <- paste(sample, "long_polyA.lenght_bin.txt", sep =".");
write.table(long_polyA_length_bins, file=table_out, quote = FALSE, sep="\t", row.names = FALSE, col.names = TRUE)

# barplot
ylim <- max(hist$counts)
pdf(long_out, width = 5, height = 4)
barplot <- barplot(hist$counts, 
                   ylim = c(0, ylim), 
                   xlab ="Poly(A) tail length", 
                   ylab = "Transcripts counts",
                   xaxt = "n")
text(barplot, par("usr")[3], 
     labels = breaks.name, 
     srt = 45, adj = c(1.2,1.1), 
     xpd = TRUE, cex= 1)
	text( barplot, hist$counts + 1, 
     labels = hist$counts, 
     srt = 0, adj = c(0.5,-1), 
     xpd = TRUE, cex= 1)
invisible(dev.off())
