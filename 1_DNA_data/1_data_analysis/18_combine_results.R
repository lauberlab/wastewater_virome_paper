# packages
library(tidyverse)

# parameters
what="DNA"

#
# A) Mapping results
#

# read count data
theData <- NULL
for ( cfile in list.files( path=paste0(what,"/mapping/"), pattern="*.counts" ) ){
	dat <- read.delim( paste0(what,"/mapping/",cfile), sep="\t", header=F, stringsAsFactors=F )
	colnames(dat) <- c("contig_id","read_count")
	tmp <- strsplit( cfile, "_" )[[1]]
	id  <- paste0( rev(rev(tmp)[-c(1:2)]), collapse="_" )
	dat <- dat %>% mutate( sample_id=id )
	theData <- rbind( theData, dat %>% dplyr::select(sample_id,contig_id,read_count) )
}
dim(theData)

# read taxonomy data
tax <- read.delim( paste0(what,"/mmseqs/mmseqs_taxonomy.tsv"), sep="\t", header=F, stringsAsFactors=F )
colnames(tax) <- c( "orf_id", "tax_id", "rank", "taxon", "taxonomy" )
tax$contig_id <- sapply( tax$orf_id, function(x){ paste0( rev(rev( strsplit(x, "_")[[1]] )[-c(1)]), collapse="_" ) } ) %>% as.character()
dim(tax)

# duplicated IDs, keep only one
tax$ucount <- sapply( tax$taxonomy, function(x){ sum( grepl("unknown",strsplit(x,";")[[1]]) ) } )
tax        <- tax %>%	group_by( contig_id ) %>%
			filter( ucount == min(ucount)[1] ) %>%
			slice(1) %>%
			ungroup()
tax        <- tax %>%	dplyr::select( contig_id, tax_id, taxonomy ) %>% unique()
dim(tax)

# contig lengths
theData$contig_length <- sapply( theData$contig_id, function(x){ as.integer( strsplit(x,"_")[[1]][4] ) } )

# combine
theData <- theData %>% left_join( tax, by="contig_id" )
dim(theData)

# save result to file
write.table( theData, file=paste0("combined_mapping_results_",what,".tsv"), sep="\t", row.names=F, quote=F )


#
# A) Diamond results
#
theIdent <- theData %>% dplyr::select( contig_id, contig_length ) %>% unique()

# read diamond results
ddir = paste0( what, "/uniqueness" )
for ( dfile in list.files(ddir,pattern="*_diamond.tsv") ){
	did = strsplit( dfile, "_" )[[1]][1]
	dda = read.delim( paste(ddir,dfile,sep="/"), sep="\t", header=F )
	colnames(dda) <- c("contig_id","contig_length","sseqid",did,"evalue","cigar")
	toselect <- c("contig_id",did)
	theIdent <- theIdent %>% left_join( dda %>% dplyr::select(any_of(toselect)), by="contig_id" )
}
theIdent <- theIdent %>% mutate( pident = apply( theIdent[,-c(1,2)], 1, max, na.rm=T ) )
theIdent <- theIdent %>% dplyr::select( contig_id, contig_length, pident )
theIdent$pident[ which(theIdent$pident == -Inf) ] <- NA
dim(theIdent)
sum(is.na(theIdent$pident))

# save result to file
write.table( theIdent, file=paste0("combined_sequence_identity_results_",what,".tsv"), sep="\t", row.names=F, quote=F )

