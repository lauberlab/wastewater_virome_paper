# packages
library( tidyverse )
library( cowplot )
library( phytools )
library( ggtree )

# input
treefile <- "data/BunyaArticulaMononega_RdRp_SRA_hits_added_known_selected_trim_LG.G4.UCLN.Yule.5Mio_annotated.nwk"
taxfile  <- "data/BunyaArticulaMononega_RdRp_trim_reference_accessions_taxonomy.tsv"
minPAT   <- 0.05

# colors
col_groups <- c( "Mononegavirales"="black", "Articulavirales"="chocolate4",
		 "Arenaviridae"="goldenrod", "Hantaviridae"="tomato", "Fimoviridae"="lightskyblue", "Nairoviridae"="plum", "Leishbuviridae"="purple",
		 "Peribunyaviridae"="darkblue", "Phasmaviridae"="steelblue", "Phenuiviridae"="chartreuse3", "Tospoviridae"="lightsalmon",
		 "new_MDC"="cyan2", "new_SRA"="magenta" )

# condense tree
tree0 <- read.tree( treefile )
system( paste("Rscript helper/tree_condense_by_distance.R", treefile, "tree_condensed.nwk", minPAT), intern=F )

# load tree
tree <- read.tree( "tree_condensed.nwk" )
tree <- drop.tip( tree, "NODE_291_length_3610_cov_2.374684_2" ) # manually remove duplicated sequence
write( tree$tip.label, file="tree_condensed.ids" )
tree0$tip.label[ ! tree0$tip.label %in% tree$tip.label ]
print( paste( length(tree0$tip.label)-length(tree$tip.label), "tips removed" )  )

# midpoint-pseudoroot tree
#tree <- midpoint.root(tree)

# read taxonomy information
tax <- read.delim( taxfile, sep="\t", header=F )
colnames(tax) <- c("accession","taxon","rank")
tax <- rbind( tax, data.frame(accession="YP_009553313.1",taxon="Phasmaviridae",rank="family") ) 

# color tree branches by host
hosts <- tree %>% as_tibble() %>% mutate( scientific_name = "internal" ) %>%
				  mutate( scientific_name = if_else( grepl("NODE_",label), "new_MDC", scientific_name ) ) %>%
				  mutate( scientific_name = if_else( grepl( "^SRR",label), "new_SRA", scientific_name ),
					  scientific_name = if_else( grepl( "^ERR",label), "new_SRA", scientific_name ),
					  scientific_name = if_else( grepl( "^DRR",label), "new_SRA", scientific_name ) )
for ( id in tax %>% filter(rank=="family") %>% pull(accession) ){
	fam   <- tax %>% filter( accession==id, rank=="family" ) %>% pull(taxon)
	fam   <- if_else( fam=="", "unclassified", fam )
	hosts <- hosts %>% mutate( scientific_name = if_else( grepl(id,label), fam, scientific_name ) )
}
for ( id in tax %>% filter(taxon=="Mononegavirales") %>% pull(accession) ){
	hosts <- hosts %>% mutate( scientific_name = if_else( grepl(id,label), "Mononegavirales", scientific_name ) )
}
for ( id in tax %>% filter(taxon=="Articulavirales") %>% pull(accession) ){
	hosts <- hosts %>% mutate( scientific_name = if_else( grepl(id,label), "Articulavirales", scientific_name ) )
}
head(hosts)
table(hosts$scientific_name)

hosts %>% filter(scientific_name=="internal",grepl("virus",label)) %>% pull(label)


grouping <- list()
for ( hh in unique(hosts$scientific_name) ){
	grouping[[ hh ]] <- hosts %>% filter(scientific_name==hh) %>% pull(label)
}

tree <- groupOTU( tree, grouping )

# produce tree plot
p1 <- ggtree(tree, aes(color=group), layout="circular" ) +
	scale_color_manual( values=col_groups )
p2 <- ggtree(tree, aes(color=group), layout="circular", branch.length='none' ) +
	scale_color_manual( values=col_groups )

# save tree plot
ggsave( p1, file="tree_condensed_circular_phylo.pdf",  limitsize=F, width=7, height=6 )
ggsave( p2, file="tree_condensed_circular_dendro.pdf", limitsize=F, width=7, height=6 )
