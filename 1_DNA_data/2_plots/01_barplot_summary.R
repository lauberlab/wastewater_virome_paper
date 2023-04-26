# packages
library(tidyverse)
library(gridExtra)
library(scales)
library(cowplot)


#
# global
#
theme_set( theme_minimal() )


#
# load data
#
theData <- read.delim( "data/virus_data_by_contig_final.tsv", sep="\t" )
head(theData)


#
# preprocess / filter data
#
theData <- theData %>% filter( ! order %in% c("Mulpavirales","uc_Pisoniviricetes","uc_Negarnaviricota","uc_Kitrinoviricota","Jingchuvirales") )


# 
# number of contigs - DNA viruses
#
thedata <- theData %>% filter( vtype=="DNAviruses" ) %>%
			group_by(order) %>% summarize( n=n(), phylum=unique(phylum) ) %>% 
			arrange(n) %>% ungroup
orderorder1 <- thedata$order
N1          <- length(orderorder1)
thedata <- thedata %>% mutate( order=factor(order,levels=orderorder1) )

p1 <- thedata %>%  ggplot( aes(x=order, y=n, fill=phylum) ) +
			geom_bar( stat="identity" ) +
			scale_y_log10( expand=c(0,0), labels=comma ) +
			coord_flip() +
			labs( x="", y="number of contigs" ) +
			theme( legend.position="none")


# 
# number of contigs - RNA viruses
#
thedata <- theData %>% filter( vtype=="RNAviruses" ) %>%
			group_by(order) %>% summarize( n=n(), phylum=unique(phylum) ) %>% 
			arrange(n) %>% ungroup
orderorder2 <- thedata$order
N2          <- length(orderorder2)
thedata <- thedata %>% mutate( order=factor(order,levels=orderorder2) )

p2 <- thedata %>%  ggplot( aes(x=order, y=n, fill=phylum) ) +
			geom_bar( stat="identity" ) +
			scale_y_log10( expand=c(0,0), labels=comma ) +
			coord_flip() +
			labs( x="", y="number of contigs" ) +
			theme( legend.position="none")


#
# contig lengths - DNA viruses
#
thedata <- theData %>% filter( vtype=="DNAviruses" ) %>% mutate( order=factor(order,levels=orderorder1) )

p3 <- thedata %>% ggplot( aes(x=order, y=contig_length, fill=phylum) ) +
			geom_violin() +
			#geom_boxplot() +
			scale_y_log10( expand=c(0,0), labels=comma, breaks=c(500,1500,5000) ) +
			coord_flip() +
			labs( x="", y="contig length" ) +
			theme( legend.position="none")


#
# contig lengths - RNA viruses
#
thedata <- theData %>% filter( vtype=="RNAviruses" ) %>% mutate( order=factor(order,levels=orderorder2) )

p4 <- thedata %>% ggplot( aes(x=order, y=contig_length, fill=phylum) ) +
			geom_violin() +
			#geom_boxplot() +
			scale_y_log10( expand=c(0,0), labels=comma, breaks=c(500,1500,5000) ) +
			coord_flip() +
			labs( x="", y="contig length" ) +
			theme( legend.position="none")


#
# sequence identity - DNA viruses
#
thedata <- theData %>% filter( vtype=="DNAviruses" ) %>% mutate( order=factor(order,levels=orderorder1) )

p5 <- thedata %>% ggplot( aes(x=order, y=identity_to_known, fill=phylum) ) +
			geom_boxplot() +
			coord_flip() +
			labs( x="", y="sequence identity [%]" ) +
			theme( legend.text=element_text(size=8) )


#
# sequence identity - RNA viruses
#
thedata <- theData %>% filter( vtype=="RNAviruses" ) %>% mutate( order=factor(order,levels=orderorder2) )

p6 <- thedata %>% ggplot( aes(x=order, y=identity_to_known, fill=phylum) ) +
			geom_boxplot() +
			coord_flip() +
			labs( x="", y="sequence identity [%]" ) +
			theme( legend.text=element_text(size=8) )


#
# abundance - DNA viruses
#
thedata <- theData %>% filter( vtype=="DNAviruses" ) %>% mutate( order=factor(order,levels=orderorder1) ) %>%
			group_by(order) %>% summarize( abundance=sum(abundance), phylum=unique(phylum) )

p7 <- thedata %>% ggplot( aes(x=order, y=(abundance+1), fill=phylum) ) +
			geom_bar( stat="identity" ) +
			scale_y_log10( expand=c(0,0), labels=comma ) +
			coord_flip() +
			labs( x="", y="abundance" ) +
			theme( legend.position="none" )


#
# abundance - RNA viruses
#
thedata <- theData %>% filter( vtype=="RNAviruses" ) %>% mutate( order=factor(order,levels=orderorder2) ) %>%
			group_by(order) %>% summarize( abundance=sum(abundance), phylum=unique(phylum) )

p8 <- thedata %>% ggplot( aes(x=order, y=(abundance+1), fill=phylum) ) +
			geom_bar( stat="identity" ) +
			scale_y_log10( expand=c(0,0), labels=comma ) +
			coord_flip() +
			labs( x="", y="abundance" ) +
			theme( legend.position="none" )


# 
# compose final multi-panel figure
#
pdf( "plot.pdf", width=12, height=8 )
#grid.arrange( p1,p7,p3,p5, p2,p8,p4,p6, ncol=4, nrow=2, widths=c(1,1,1,1.65), heights=c(N1,N2*0.9) )
plot_grid( p1,p7,p3,p5, p2,p8,p4,p6, ncol=4, nrow=2, rel_widths=c(1,1,1,1.65), rel_heights=c(N1,N2*0.9),
		align="v", axis="l", labels=c("A","B","C","D","E","F","G","H"), hjust=-1 )
dev.off()

