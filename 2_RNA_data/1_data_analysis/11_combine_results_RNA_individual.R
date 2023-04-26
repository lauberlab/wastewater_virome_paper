# packages
library(tidyverse)

# parameters
what="RNA"

#
# A) Mapping results
#

# read count data
theData <- NULL
for ( cfile in list.files( path=paste0(what,"/mapping_individual/"), pattern="*.counts" ) ){
	dat <- read.delim( paste0(what,"/mapping_individual/",cfile), sep="\t", header=F, stringsAsFactors=F )
	colnames(dat) <- c("contig_id","read_count")
	cfile2 <- sub( "_L001", "", cfile )
	tmp <- strsplit( cfile2, "_" )[[1]]	
	id  <- paste0( rev(rev(tmp)[-1]), collapse="_" )
	dat <- dat %>% mutate( sample_id=id )
	theData <- rbind( theData, dat %>% dplyr::select(sample_id,contig_id,read_count) )
}
dim(theData)


# read taxonomy data
tax <- NULL
for ( tdir in list.files( path=paste0(what,"/assemblies/"), pattern="*" ) ){
	tax0 <- read.delim( paste0(what,"/assemblies/",tdir,"/mmseqs/mmseqs_taxonomy.tsv"), sep="\t", header=F, stringsAsFactors=F )
	tmp1 <- sub( "_L001", "", tdir )
	tmp2 <- strsplit( tmp1, "_" )[[1]]
	id   <- paste0( rev(rev(tmp2)[-1]), collapse="_" )
	tax0 <- cbind( tax0, rep(id,nrow(tax0)) )
	tax  <- rbind(tax,tax0)
}
colnames(tax) <- c( "orf_id", "tax_id", "rank", "taxon", "taxonomy", "sample_id" )
tax$contig_id <- sapply( tax$orf_id, function(x){ paste0( rev(rev( strsplit(x, "_")[[1]] )[-c(1)]), collapse="_" ) } ) %>% as.character()
dim(tax)

# duplicated IDs, keep only one
tax$ucount <- sapply( tax$taxonomy, function(x){ sum( grepl("unknown",strsplit(x,";")[[1]]) ) } )
tax        <- tax %>%	group_by( contig_id ) %>%
			filter( ucount == min(ucount)[1] ) %>%
			slice(1) %>%
			ungroup()
tax        <- tax %>%	dplyr::select( sample_id, contig_id, tax_id, taxonomy ) %>% unique()
dim(tax)

# contig lengths
theData$contig_length <- sapply( theData$contig_id, function(x){ as.integer( strsplit(x,"_")[[1]][4] ) } )

# combine
tax     <-     tax %>% mutate( uid=paste(sample_id,contig_id,sep=":") )
theData <- theData %>% mutate( contig_id2 = sub( "_s\\d+$", "", contig_id ) ) %>% 
		       mutate( uid=paste(sample_id,contig_id2,sep=":") ) %>% dplyr::select(-contig_id2)
theData <- theData %>% left_join( tax %>% dplyr::select(-contig_id,-sample_id), by="uid" )
theData <- theData %>% dplyr::select(-uid)
dim(theData)

print( "data in tax:" )
sum( theData$uid %in% tax$uid )
print( "tax in data:" )
sum( tax$uid %in% theData$uid )

# save result to file
write.table( theData, file=paste0("combined_mapping_results_",what,".tsv"), sep="\t", row.names=F, quote=F )


#
# A) Diamond results
#
theIdent <- theData %>% dplyr::select( contig_id, contig_length ) %>% unique()

# read diamond results
ddir = paste0( what, "/uniqueness_individual" )
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

