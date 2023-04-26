#!/usr/bin/perl

# input
my $data    = shift; # DNA or RNA
my $threads = shift; # number of threads
my $ram     = shift; # RAM per thread
my $errcorr = shift; # do or do not run error correction in spades
my $kmers   = shift; # k-mer sizes for spades

# read datasets
my @datasets = ();
open( DS, "<datasets_${data}.txt" );
while ( my $line = <DS> ){
	chomp($line);  next if $line eq "";
	my @v = split( " ", $line );
	push( @datasets, $v[1] );
}
close(DS);

# init spades command
my $cmd = "spades --checkpoints last ";

# k mers
my @ks = split( /,/, $kmers );
my $spades_k = "-k";
foreach ( @ks ){ $spades_k .= " $_"; }
$cmd .= " $spades_k";

# compile spades command
my $ddir = "${data}/fastq";
for ( my $i=0; $i<=$#datasets; $i++ ){
	$cmd .= sprintf " --pe-1 %d %s/%s", $i+1, $ddir, $datasets[$i]."_R1_trim.fastq.gz";
	$cmd .= sprintf " --pe-2 %d %s/%s", $i+1, $ddir, $datasets[$i]."_R2_trim.fastq.gz";
}
if ( $errcorr == 1 ){
	$cmd .= sprintf " -t %d -m %d -o %s/spades", $threads, $ram*$threads, $data;
}else{
	$cmd .= sprintf " -t %d -m %d -o %s/spades_NoErrCor --only-assembler", $threads, $ram*$threads, $data;
}

# run spades
`$cmd`;

