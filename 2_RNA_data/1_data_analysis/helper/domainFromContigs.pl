#!/usr/bin/perl

# ------------ #
# load modules
# ------------ #
use warnings;
use strict;
use Data::Dumper;


# --------- #
# parameter
# --------- #
my $getorf  = "/home/lauber/lib/EMBOSS-6.5.7/emboss/getorf";
my $minorf  = 150;
my $threads =   6;
my $E_th    = 100;


# --------------- #
# usage and input
# --------------- #
if ( $#ARGV < 3 ){
	die("\nusage: virusutils-SortByDomain-nucl.pl <contig_file> <profile_file> <outdir> <min_contig_len> [-e=<E_cutoff>]\n\n");
}
my $cFile  = shift;
my $pFile  = shift;
my $outdir = shift;
my $L_th   = shift;

while ( $#ARGV > -1 ){
	my $option = shift;
	if ( $option =~ /-e=(.*)/   ){  $E_th = $1; }
	#if ( $option =~ /-rsv=(.*)/ ){  $reserv = $1; }
}


# ------- #
# init
# ------- #
#my $timestamp = `date +%y%m%d%H%M%S`;  chomp($timestamp);
#my $outdir0   = $outdir;
#   $outdir   .= "_".$timestamp;

mkdir( "$outdir" )  if ( ! -d $outdir );
#`rm -rf $outdir/1_specific $outdir/2_grouped`;
mkdir( "$outdir/1_specific" )  if ( ! -d "$outdir/1_specific" );
`rm -rf $outdir/2_grouped`;
mkdir( "$outdir/2_grouped" );
mkdir( "$outdir/2_grouped/unclassified" );


# --------------------- #
# get all profile names
# --------------------- #
my @pNames = ();
my $pcmd   = "grep NAME $pFile | awk '{ print \$2}'";
open( P, "$pcmd |" ) or die( "Can't execute '$pcmd': $!\n" );
while ( my $line = <P> ){
	chomp( $line );
	push( @pNames, $line );
}
close(P);

foreach my $pn ( @pNames ){
	mkdir( "$outdir/2_grouped/$pn" );
}
#print Dumper( @pNames );


# ---------------------------- #
# read contigs from fasta file
# ---------------------------- #
my %contigs = ();
my $cid;
open( C, "<$cFile" ) or die( "Can't open file '$cFile': $!\n" );
while ( my $line = <C> ){
	chomp( $line );
	next if $line eq "";
	if ( $line =~ />(.*)/ ){
		$cid = clean_seq_names( $1 );
		$contigs{ $cid }  = "";
	}else{
		$contigs{ $cid } .= $line;
	}
}
close(C);

foreach $cid ( keys %contigs ){
	delete $contigs{ $cid }  if ( length( $contigs{ $cid } ) < $L_th );
}



# ------------------ #
# do for each contig
# ------------------ #

# HMMer searches
my $calculated = 0;
my $noorfs     = 0;
foreach $cid ( sort keys %contigs ){
	# dir
	my $cdir = "$outdir/1_specific/$cid";

	# test if already done
	next if ( -e "$cdir/$cid.orfs" and -e "$cdir/hmmsearch.tsv" );

	# create dir
	mkdir( $cdir );
	
	# write indiviual contig to file
	open(  CS, ">$cdir/$cid.fas" ) or die( "Can't write to file '$cdir/$cid.fas': $!\n" );
	printf CS ">%s\n%s\n", $cid, $contigs{ $cid };
	close( CS );

	# extract ORFs
	`$getorf -sequence $cdir/$cid.fas -outseq $cdir/$cid.orfs -minsize $minorf 2> $cdir/$cid.orfs.err`;
	$noorfs++  if (   -e "$cdir/$cid.orfs" and -s "$cdir/$cid.orfs" != 0 );
	next       if ( ! -e "$cdir/$cid.orfs"  or -s "$cdir/$cid.orfs" == 0 );
	
	# run hmmsearch
	`hmmsearch -E $E_th --cpu $threads --domtblout $cdir/hmmsearch.tsv $pFile $cdir/$cid.orfs > $cdir/hmmsearch.txt`;
	$calculated++ if ( -e "$cdir/hmmsearch.tsv" );
}

# report what was done
printf "\nhmmsearch results newly calculated for %d sequences\n", $calculated;
printf "No ORFs >%d nt found for %d sequences\n", $minorf, $noorfs;



# collect results
my %summary  = ();
my $missing  = 0;
my $multiple = 0;
foreach $cid ( sort keys %contigs ){
	# dir
	my $cdir = "$outdir/1_specific/$cid";
	
	# init
	$summary{ $cid } = ();

	# test if HMMer search completed
	if ( ! -e "$cdir/hmmsearch.tsv" ){
		$missing++;
		next;
	}	

	# gather best hit for each ORF
	my %bestE = ();
	open( HMA, "<$cdir/hmmsearch.tsv" ) or die( "Can't read file '$cdir/hmmsearch.tsv': $!\n" );
	while ( my $line = <HMA> ){
		chomp( $line );
		next if $line eq "";
		next if $line =~ /^#.*/;
		my @v = split( /\s+/, $line );
		my $orfid = $v[0];
		my $hmmid = $v[3];
		if ( ! exists $summary{ $cid }{ $orfid } ){
			$summary{ $cid }{ $orfid } = ();
			foreach my $pn ( @pNames ){ $summary{ $cid }{ $orfid }{ $pn } = "NA"; }
			$summary{ $cid }{ $orfid }{ $hmmid } = $v[6];
		}
		if    ( "NA" eq $summary{ $cid }{ $orfid }{ $hmmid } ){
			$summary{ $cid }{ $orfid }{ $hmmid } = $v[6];
		}elsif( $v[6] < $summary{ $cid }{ $orfid }{ $hmmid } ){
			$summary{ $cid }{ $orfid }{ $hmmid } = $v[6];
		}
	}
	close(HMA);
	
	# read aa seq from file
	open( AA, "<$cdir/$cid.orfs" ) or die( "Can't read file '$cdir/$cid.orfs': $!\n" );
	my %orfs = ();
	my $orfid;
	while ( my $line = <AA> ){
		chomp( $line );
		next if $line eq "";
		if ( $line =~ /^>([^ ]+) .*/ ){
			$orfid = $1;
			$orfs{ $orfid }  = "";
		}else{
			$orfs{ $orfid } .= $line;
		}
	}
	close(AA);

	# sort ORF according to best profile hit
	if ( scalar keys %{$summary{ $cid }} == 0 ){
		open(  RNT, ">>$outdir/2_grouped/unclassified/nucleotide.fas" );
		printf RNT ">%s\n%s\n", $cid, $contigs{ $cid };
		close( RNT );
	}else{
		my %nodupdup = ();
		my $orfhitN  = 0;
		foreach my $orfid ( sort keys %{$summary{ $cid }} ){
			$orfhitN++;
			$multiple++ if ( $orfhitN == 2 );

			my ($bestP,$bestE) = ("",999999);
			foreach my $hmmid ( keys %{$summary{ $cid }{ $orfid }} ){
				next if ( $summary{ $cid }{ $orfid }{ $hmmid } eq "NA" );
				if ( $summary{ $cid }{ $orfid }{ $hmmid } < $bestE ){
					$bestE = $summary{ $cid }{ $orfid }{ $hmmid };
					$bestP = $hmmid;
				}
			}
			
			if ( ! exists $nodupdup{ $bestP } ){
				open(  RNT, ">>$outdir/2_grouped/$bestP/nucleotide.fas" );
				printf RNT ">%s\n%s\n", $cid, $contigs{ $cid };
				close( RNT );
				$nodupdup{ $bestP } = 1;
			}else{
				$nodupdup{ $bestP }++;
			}
		
			my $outid = $cid;
			   $outid = sprintf "%s_%d", $outid, $orfhitN if( $orfhitN > 1 );
			open(  RAA, ">>$outdir/2_grouped/$bestP/protein.fas" );
			printf RAA ">%s\n%s\n", $outid, $orfs{ $orfid };
			close( RAA );
		}
	}
}

# report finishing status
printf "%d sequences received hits against two or more ORFs.\n", $multiple;
printf "hmmsearch results are missing for %d sequences.", $missing;
printf "  TRY RUNNING ME AGAIN ?" if ($missing > 0);
printf "\n\n";



# ---------------------- #
# compile result summary
# ---------------------- #
#print Dumper( %summary );
my $sFile = "$outdir/3_summary.tsv";
open( SUM, ">$sFile" ) or die( "Can't write to file '$sFile': $!\n" );
printf SUM "contig";
foreach ( sort @pNames ){ printf SUM "\t%s", $_; }
printf SUM "\n";
foreach $cid ( sort keys %summary ){
	my $orfhitN = 0;
	foreach my $orfid ( sort keys %{$summary{ $cid }} ){
		$orfhitN++;

		my $outid = $cid;
		   $outid = sprintf "%s_%d", $outid, $orfhitN if( $orfhitN > 1 );

		printf SUM "%s", $outid;
		foreach ( sort @pNames ){
			my $eval = $summary{ $cid }{ $orfid }{ $_ };
			   $eval = "" if ( $eval eq "NA" );
			printf SUM "\t%s", $eval;
		}
		printf SUM "\n";
	}
}
close(SUM);


# -------- #
# finalize
# -------- #
#`rm -rf $outdir0`;
#`mv $outdir $outdir0`;


# ------------ #
# helper stuff
# ------------ #
sub clean_seq_names{
	my $id = shift;

	#$id =~ s/, complete genome//g;
	#$id =~ s/, /_/g;
	#$id =~ s/,/_/g;
	#$id =~ s/; /_/g;
	#$id =~ s/;/_/g;
	#$id =~ s/: /_/g;
	#$id =~ s/:/_/g;
	#$id =~ s/\(/_/g;
	#$id =~ s/\)/_/g;
	#$id =~ s/\\/_/g;
	#$id =~ s/\//_/g;
	#$id =~ s/'/_/g;
	#$id =~ s/"/_/g;
	#$id =~ s/ /_/g;
	$id =~ s/ .*//;
	#$id =~ s/\.\d+$//;
	
	return( $id );
}


### END ###

