# packages
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(scales)

# global stuff
theme_set( theme_minimal() )

# directories
odir  <- "RNA/plots_forPaper"
DEBUG <- F 
PUNIT <- 6

# load mapping data
theData <- NULL
if ( DEBUG ){
	theData <- read.delim( "combined_mapping_results_RNA_shuf10k.tsv",  sep="\t", stringsAsFactors=F )
}else{
	theData <- read.delim( "combined_mapping_results_RNA.tsv",          sep="\t", stringsAsFactors=F )
}

# load and join total read count data
theCnt  <- read.delim( "combined_total_read_count_RNA.tsv", sep="\t", header=F, stringsAsFactors=F )
colnames(theCnt) <- c("sample_id","total_read_count")
theData <- theData %>% left_join( theCnt, by="sample_id" )
dim(theData)


# sum by date-time
theData$date <- sapply( theData$sample_id, function(x){ paste(strsplit(x,"_")[[1]][2:3],collapse="_") } )
theData$date <- sub( "mix", "", theData$date )

theTota <- theData %>% dplyr::select(date,sample_id,total_read_count) %>% unique
theTota <- theTota %>% group_by(date) %>% summarize( total_read_count = sum(total_read_count,na.rm=T) ) %>% ungroup
theData <- theData %>% group_by( date, contig_id, contig_length, tax_id, taxonomy ) %>%
			    summarize( read_count = sum(read_count,na.rm=T),
				         sample_n = n() ) %>%
			    ungroup() %>%
			    left_join( theTota, by="date" )

# extract taxonomy info
theData <- theData %>% separate( taxonomy, sep=";", remove=F,
				 into=c("superkingdom","kingdom","phylum","subphylum","class","order","suborder","family","subfamily","genus","subgenus","species") )
for ( ti in c("superkingdom","kingdom","phylum","subphylum","class","order","suborder","family","subfamily","genus","subgenus","species") ){
	theData[ which(is.na(theData[,ti])), ti]        <- "unclassified"
	theData[ which(theData[,ti]=="uc_Viruses"), ti] <- "unknown"
}

# remove contigs with missing taxonomy information
theData <- theData %>% filter( superkingdom != "unclassified", superkingdom != "unknown", superkingdom != "",
				     phylum != "unclassified",       phylum != "unknown",       phylum != "",
			   	      order != "unclassified",        order != "unknown",        order != "" )

# info per contig
theContigs <- theData %>% dplyr::select( contig_id, contig_length, tax_id, kingdom, phylum, order, family, species ) %>% unique()

# sequence identity to closest known virus
theIdent <- NULL
if ( DEBUG ){
	theIdent <- read.delim(  "combined_sequence_identity_results_RNA.tsv",  sep="\t", stringsAsFactors=F )
	theIdent <- theIdent %>% left_join( theContigs %>% dplyr::select(contig_id,phylum,order,family,species), by="contig_id" )
}else{
	theIdent <- read.delim(  "combined_sequence_identity_results_RNA.tsv",  sep="\t", stringsAsFactors=F )
	theIdent <- theIdent %>% left_join( theContigs %>% dplyr::select(contig_id,phylum,order,family,species), by="contig_id" )
}
head(theIdent)

# highest abundance
theAbu <- theData %>% group_by(contig_id) %>% 
			summarize( abundance = sum( read_count / total_read_count / contig_length * 1000000 ) %>% round(6) ) %>% 
			ungroup() %>%
			arrange( desc(abundance) ) %>%
			filter( !is.na(abundance) ) 

# add sequence identity
theAbu <- theAbu %>% left_join( theIdent %>% dplyr::select(contig_id,identity_to_known=pident) %>% unique(), by="contig_id" )

# add taxonomy
theAbu <- theAbu %>% left_join( theData  %>% dplyr::select(contig_id,contig_length,phylum,order,family,genus,species) %>% unique, by="contig_id" ) 

# discriminate RNA and DNA viruses
RNA_phyla <- c("Lenarviricota","Pisuviricota","Kitrinoviricota","Duplornaviricota","Negarnaviricota")
theAbu <- theAbu %>% mutate( vtype = ifelse( phylum %in% RNA_phyla,"RNAviruses","DNAviruses" ) )


# summarize by order
theAbuSum <- theAbu %>% group_by(order) %>% summarize( phylum=unique(phylum), vtype=unique(vtype),
							abundance = sum(abundance) %>% round(6), 
							n_contigs = n(),
							identity_to_known_q0  = quantile(identity_to_known,na.rm=T)[1] %>% round(1),
							identity_to_known_q25 = quantile(identity_to_known,na.rm=T)[2] %>% round(1),
							identity_to_known_q50 = quantile(identity_to_known,na.rm=T)[3] %>% round(1),
							identity_to_known_q75 = quantile(identity_to_known,na.rm=T)[4] %>% round(1),
							identity_to_known_q100= quantile(identity_to_known,na.rm=T)[5] %>% round(1),
							contig_length_q0  = quantile(contig_length,na.rm=T)[1] %>% round(1),
							contig_length_q25 = quantile(contig_length,na.rm=T)[2] %>% round(1),
							contig_length_q50 = quantile(contig_length,na.rm=T)[3] %>% round(1),
							contig_length_q75 = quantile(contig_length,na.rm=T)[4] %>% round(1),
							contig_length_q100= quantile(contig_length,na.rm=T)[5] %>% round(1) )


# save data tables
write.table( theAbu,    file=paste0(odir,"/virus_data_by_contig.tsv"), sep="\t", row.names=F, quote=F )
write.table( theAbuSum, file=paste0(odir,"/virus_data_by_order.tsv"),  sep="\t", row.names=F, quote=F )

q()







#
# plot virus counts
#
plot_taxa_count <- function( data, rank_plot, rank_fill, outdir, log_scale=F, t_per_h=6, vtype="RNA" ){

	# count taxa per rank
	rN <- data %>% pull(!!rank_plot) %>% unique() %>% length()

	# order by higher level
	rx <- data %>% arrange(!!sym(rank_fill),!!sym(rank_plot)) %>% pull(!!rank_plot) %>% unique()

	# by virus phylum
	p  <- data %>% mutate( !!sym(rank_plot) := factor( !!sym(rank_plot), levels=rx ) ) %>%
				ggplot( aes(x=!!sym(rank_plot),fill=!!sym(rank_fill)) ) + 
				geom_bar( position="dodge" ) +
				labs( y="# contigs",x="",title=paste("virus",rank_plot) ) +
				coord_flip() +
				theme( legend.text=element_text(size=7) ) +
				guides( fill=guide_legend(ncol=2) )

	# optionally, plot in log10-scale
	if ( log_scale ){
		p  <- p + scale_y_log10( expand=c(0,0), labels=comma )
	}
	
	# save
		pdffile <- paste0(outdir,"/virus_count_",rank_plot,"_",vtype,".pdf")
	if ( log_scale ){
		pdffile <- paste0(outdir,"/virus_count_",rank_plot,"_",vtype,"_log10.pdf")
	}
	ggsave( pdffile, p, width=8, height=(2+rN/t_per_h) )
}

plot_taxa_count( theAbuRNA, "order",  "phylum",  odir, t_per_h=PUNIT, log_scale=T, vtype="RNAviruses" )
plot_taxa_count( theAbuDNA, "order",  "phylum",  odir, t_per_h=PUNIT, log_scale=T, vtype="DNAviruses" )

q()


#
# plot taxa abundance
#
plot_taxa_abundance <- function( data, rank_plot, rank_fill, outdir, log_scale=F, t_per_h=6 ){

	# count taxa per rank
	rN <- data %>% pull(!!rank_plot) %>% unique() %>% length()

	# order by higher level
	rx <- data %>% arrange(!!sym(rank_fill),!!sym(rank_plot)) %>% pull(!!rank_plot) %>% unique()

	# adjust abundance in case of log scale plot
	if ( log_scale ){
	      data <- data %>% mutate( abundance=abundance+1 )
	}

	# by virus phylum
	p  <- data %>% mutate( !!sym(rank_plot) := factor( !!sym(rank_plot), levels=rx ) ) %>%
				ggplot( aes(x=!!sym(rank_plot),y=abundance,fill=!!sym(rank_fill)) ) + 
				geom_bar( position="dodge", stat="identity" ) +
				labs( y="abundance",x="",title=paste("virus",rank_plot) ) +
				coord_flip() +
				theme( legend.text=element_text(size=7) ) +
				guides( fill=guide_legend(ncol=2) )

	# optionally, plot in log10-scale
	if ( log_scale ){
		p  <- p + scale_y_log10( expand=c(0,0), labels=comma )
	}
	
	# save
		pdffile <- paste0(outdir,"/virus_abundance_",rank_plot,".pdf")
	if ( log_scale ){
		pdffile <- paste0(outdir,"/virus_abundance_",rank_plot,"_log10.pdf")
	}
	ggsave( pdffile, p, width=8, height=(2+rN/t_per_h) )
}

plot_taxa_abundance( theAbu, "phylum", "kingdom", odir, t_per_h=PUNIT )
plot_taxa_abundance( theAbu, "order",  "phylum",  odir, t_per_h=PUNIT )
plot_taxa_abundance( theAbu, "family", "phylum",  odir, t_per_h=PUNIT )
plot_taxa_abundance( theAbu, "order",  "phylum",  odir, t_per_h=PUNIT, log_scale=T )
plot_taxa_abundance( theAbu, "family", "phylum",  odir, t_per_h=PUNIT, log_scale=T )



#
# plot virus abundance
#
# by virus all
p1 <- theAbu %>% filter( species != "unclassified", species != "unknown" ) %>%
			ggplot( aes(x=1, y=abundance) ) +
			geom_violin( fill="black", alpha=0.4 ) + scale_y_log10(n.breaks=10,expand=c(0,0), labels=comma ) +
			labs(y="abundance",x="",title="viral contigs - all") + theme( axis.text.x=element_blank() )

# by virus top
p2 <- theAbu %>% filter( species != "unclassified", species != "unknown" ) %>%
			arrange( abundance ) %>% tail( n=50 ) %>%
			mutate( contig_id = factor(contig_id,levels=contig_id) ) %>%
			ggplot( aes(x=contig_id,y=abundance,fill=species,label=species) ) +
			geom_bar(stat="identity") + geom_label( size=2, hjust=0, nudge_y=2, alpha=0.2 ) +
			expand_limits( y=c(0,max(theAbu$abundance)*1.25) ) +
			labs(y="abundance",x="",title="viral contigs - top 50") + coord_flip() +
			theme( axis.text.y=element_text(size=5), legend.position="none" ) +
			guides( fill=guide_legend(ncol=1) )

# save
pN <- 50 
pdf( paste(odir,"virus_abundance_virus.pdf",sep="/"), width=8, height=2+(pN/PUNIT) )
grid.arrange( p1, p2, ncol=2, widths=c(0.25,0.75) )
dev.off()


# repeat without phages
pdf( paste(odir,"virus_abundance_virus_noPhages.pdf",sep="/"), width=8, height=10 )

# define virus groups
euphyla <- c("Cossaviricota","Cressdnaviricota","Nucleocytoviricota","Preplasmiviricota","Peploviricota",
             "Negarnaviricota","Duplornaviricota","Kitrinoviricota","Pisuviricota","Artverviricota")

# by virus all
p1 <- theAbu %>% filter( phylum %in% euphyla ) %>%
			ggplot( aes(x=1, y=abundance) ) +
			geom_violin( fill="black", alpha=0.4 ) + scale_y_log10(n.breaks=10,expand=c(0,0), labels=comma ) +
			labs(y="abundance",x="",title="viral contigs - all") + theme( axis.text.x=element_blank() )

# by virus top
p2 <- theAbu %>% filter( phylum %in% euphyla ) %>%
			arrange( abundance ) %>% tail( n=50 ) %>%
			mutate( contig_id = factor(contig_id,levels=contig_id) ) %>%
			ggplot( aes(x=contig_id,y=abundance,fill=species,label=species) ) +
			geom_bar(stat="identity") + geom_label( size=2, hjust=0, nudge_y=2, alpha=0.2 ) +
			expand_limits( y=c(0,max(theAbu$abundance)*1.25) ) +
			labs(y="abundance",x="",title="viral contigs - top 50") + coord_flip() +
			theme( axis.text.y=element_text(size=5), legend.position="none" ) +
			guides( fill=guide_legend(ncol=1) )

# save
pN <- 50 
pdf( paste(odir,"virus_abundance_virus_noPhages.pdf",sep="/"), width=8, height=2+(pN/PUNIT) )
grid.arrange( p1, p2, ncol=2, widths=c(0.25,0.75) )
dev.off()



#
# plot order abundance per sample
#
# by order
theAbu2 <- theData %>% mutate( abundance = read_count / contig_length / total_read_count * 1000000 ) %>%
			group_by(date,order) %>% summarize( abundance=sum(abundance) ) %>% ungroup() %>%
			arrange( desc(abundance) ) %>%
			filter( !is.na(abundance) ) %>%
			filter( grepl("virales",order) )

theAbu2   <- theAbu2 %>% left_join( theData %>% dplyr::select(order,phylum) %>% unique(), by="order" )
facetOrd2 <- theAbu2 %>% dplyr::select(order,phylum) %>% unique() %>% arrange(phylum,order) %>% pull(order)
theAbu2   <- theAbu2 %>% mutate( order=factor(order,levels=facetOrd2) )

facetN  <- theAbu2 %>% pull(order) %>% unique() %>% length()
theAbu2 %>% group_by(order) %>% summarize( abundance=sum(abundance) )

p1 <- theAbu2 %>% ggplot( aes(x=date,y=abundance,fill=phylum) ) + geom_bar(stat="identity") +
			coord_flip() +
			facet_wrap( ~order, ncol=facetN ) +
			theme( axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), legend.position="bottom" ) +
			labs( x="" ) +
			guides( fill=guide_legend(nrow=1) )

p2 <- theAbu2 %>% ggplot( aes(x=date,y=abundance,fill=phylum) ) + geom_bar(stat="identity") +
			coord_flip() +
			facet_wrap( ~order, ncol=facetN, scales="free_x" ) +
			theme( axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), legend.position="bottom" ) +
			labs( x="" ) +
			guides( fill=guide_legend(nrow=1) )

# save
pNw <- theAbu2 %>% pull(order) %>% unique() %>% length()
pNh <- theAbu2 %>% pull(date)  %>% unique() %>% length()
pdf( paste(odir,"per_sample_abundance_order.pdf",sep="/"), width=2+(pNw/PUNIT*5), height=(2+(pNh/PUNIT))*2 )
grid.arrange( p1, p2, nrow=2 )
dev.off()


#
# plot family abundance per sample
#
# by family
theAbu3 <- theData %>% mutate( abundance = read_count / contig_length / total_read_count * 1000000 ) %>%
			group_by(date,family) %>% summarize( abundance=sum(abundance) ) %>% ungroup() %>%
			arrange( desc(abundance) ) %>%
			filter( !is.na(abundance) ) %>%
			filter( grepl("viridae",family) )

theAbu3   <- theAbu3 %>% left_join( theData %>% dplyr::select(family,phylum) %>% unique(), by="family" )
facetOrd3 <- theAbu3 %>% dplyr::select(family,phylum) %>% unique() %>% arrange(phylum,family) %>% pull(family)
theAbu3   <- theAbu3 %>% mutate( family=factor(family,levels=facetOrd3) )

facetN  <- theAbu3 %>% pull(family) %>% unique() %>% length()
theAbu3 %>% group_by(family) %>% summarize( abundance=sum(abundance) )

p1 <- theAbu3 %>% ggplot( aes(x=date,y=abundance,fill=phylum) ) + geom_bar(stat="identity") +
			coord_flip() +
			facet_wrap( ~family, ncol=facetN ) +
			theme( axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), legend.position="bottom" ) +
			labs( x="" ) +
			guides( fill=guide_legend(nrow=1) )

p2 <- theAbu3 %>% ggplot( aes(x=date,y=abundance,fill=phylum) ) + geom_bar(stat="identity") +
			coord_flip() +
			facet_wrap( ~family, ncol=facetN, scales="free_x" ) +
			theme( axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), legend.position="bottom" ) +
			labs( x="" ) +
			guides( fill=guide_legend(nrow=1) )

# save
pNw <- theAbu3 %>% pull(family) %>% unique() %>% length()
pNh <- theAbu3 %>% pull(date)   %>% unique() %>% length()
pdf( paste(odir,"per_sample_abundance_family.pdf",sep="/"), width=2+(pNw/PUNIT*5), height=(2+(pNh/PUNIT))*2 )
grid.arrange( p1, p2, nrow=2 )
dev.off()


#
# plot protein sequence identity to closest known virus by order
#
# by order
theIdent2 <- theIdent  %>% mutate( order=factor(order,levels=facetOrd2) ) %>% filter( !is.na(order) )
theIdent2 <- theIdent2 %>% mutate( pident = if_else(is.na(pident),0,pident) )

p1 <- theIdent2 %>% ggplot( aes(x=order, y=pident, fill=phylum) ) +
			geom_boxplot() +
			coord_flip() +
			labs( x="", y="protein sequence identity to reference" ) +
			theme( legend.text=element_text(size=7) ) +
			guides( fill=guide_legend(ncol=2) )

# save
pN <- theIdent2 %>% pull(order) %>% unique() %>% length()
pdf( paste(odir,"sequence_identity_order.pdf",sep="/"), width=6, height=2+(pN/PUNIT) )
p1
dev.off()


#
# plot protein sequence identity to closest known virus by family
#
# by order
theIdent3 <- theIdent  %>% mutate( family=factor(family,levels=facetOrd3) ) %>% filter( !is.na(family) )
theIdent3 <- theIdent3 %>% mutate( pident = if_else(is.na(pident),0,pident) )

p1 <- theIdent3 %>% ggplot( aes(x=family, y=pident, fill=phylum) ) +
			geom_boxplot() +
			coord_flip() +
			labs( x="", y="protein sequence identity to reference" ) +
			theme( legend.text=element_text(size=7) ) +
			guides( fill=guide_legend(ncol=2) )

# save
pN <- theIdent3 %>% pull(family) %>% unique() %>% length()
pdf( paste(odir,"sequence_identity_family.pdf",sep="/"), width=6, height=2+(pN/PUNIT) )
p1
dev.off()


#
# viral vs nonviral read count
#
# barplotting
theVirCnt <- theData %>% filter( superkingdom != "unclassified",
				 superkingdom != "unknown",
				 superkingdom != "" ) %>%		
			group_by(date) %>%
			summarize(read_count=sum(read_count)) %>%
			ungroup() %>%
			left_join( theData %>% dplyr::select(date,total_read_count) %>% unique, by="date" ) %>%
			mutate( diff=total_read_count-read_count ) %>%
			mutate( diff=diff/1e6, read_count=read_count/1e6 ) %>%
			dplyr::select( date, viral=read_count, unknown=diff ) %>%
			gather( key="origin", value="read_count", -date )

p1 <- theVirCnt %>% ggplot( aes(x=date, y=read_count, fill=origin) ) +
			geom_bar( stat="identity" ) +
			scale_fill_manual( values=c("lightsteelblue","steelblue") ) +
			coord_flip() +
			theme( legend.position="top" ) +
			labs(x="", y="read count [millions]")

# save
pN <- theVirCnt %>% pull(date) %>% unique() %>% length()
pdf( paste(odir,"viral_vs_nonviral_reads.pdf",sep="/"), width=6, height=2+(pN/PUNIT) )
p1
dev.off()


#
# plot contig length by order
#
# by order
p1 <- theData %>% filter( grepl("virales",order) ) %>% mutate( order=factor(order,levels=facetOrd2) ) %>%
			ggplot( aes(x=order, y=contig_length, fill=phylum) ) +
			geom_boxplot() +
			coord_flip() +
			scale_y_log10( n.breaks=10,expand=c(0,0), labels=comma ) +
			labs( x="", y="contig length" ) +
			theme( axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
			       legend.text=element_text(size=7) ) +
			guides( fill=guide_legend(ncol=2) )

# save
pN <- length(facetOrd2)
pdf( paste(odir,"sequence_length_order.pdf",sep="/"), width=6, height=2+(pN/PUNIT) )
p1
dev.off()


