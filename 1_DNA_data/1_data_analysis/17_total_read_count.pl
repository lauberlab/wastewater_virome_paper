#!/usr/bin/perl

chdir( "DNA/fastq" );
open(OUT,">../../combined_total_read_count_DNA.tsv");
foreach my $i ( glob("*_R1_trim.fastq.gz") ){
	my $id =  $i;
	   $id =~ s/_R1_trim\.fastq\.gz//;
	   $id =~ s/_L\d+$//;
	   $id =~ s/_S\d+$//;
	my $cn = `echo \$(zcat $i | wc -l)/4 | bc`;
	chomp($cn);
	printf OUT "%s\t%d\n", $id, $cn;
}
close(OUT);


chdir( "../../RNA/fastq" );
open(OUT,">../../combined_total_read_count_RNA.tsv");
foreach my $i ( glob("*_R1_trim.fastq.gz") ){
	my $id =  $i;
	   $id =~ s/_R1_trim\.fastq\.gz//;
	   $id =~ s/_L\d+$//;
	   $id =~ s/_S\d+$//;
	my $cn = `echo \$(zcat $i | wc -l)/4 | bc`;
	chomp($cn);
	printf OUT "%s\t%d\n", $id, $cn;
}
close(OUT);

