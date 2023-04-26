#!/usr/bin/perl

# load modules
use warnings;
use strict;
#use Data::Dumper;

# input
if ( $#ARGV < 1 or $#ARGV > 2 ){
	die( "\nvirusutils-extractFromFasta.pl <fasta_file> <id_file> [-invert -maxL=<int>]\n\n" );
}
my $infile = shift;
my $idfile = shift;
my $invert = 0;
my $maxLen = 1e12;

while ( $#ARGV >=0 ){
	my $option = shift;
	$invert = 1 if $option eq "-invert";
	if( $option =~ /-maxL=(\d+)/ ){ $maxLen = $1; }
}

# read ids
my %ids = ();
open( IDS, "<$idfile" );
while ( my $line = <IDS> ){
	chomp($line);  next if $line eq ""; 
	$line       =~ s/ .*//g;
	$ids{$line} = 1;
}
close(IDS);
printf STDERR "read %d IDs\n", scalar keys %ids;


# read fasta
my %fas = ();
my %ids_orig = ();
my $id  = ""; my $id0 = "";
open( FAS, "<$infile" );
while ( my $line = <FAS> ){
	chomp($line);  next if $line eq "";
	if ( $line =~ />(.*)/ ){

		if ( exists $fas{$id} ){
			delete $fas{$id}  if( length($fas{$id}) > $maxLen );
		}

		$id        = $1;
		$id0       = $id;
		#$id        =~ s/\s//g;
		$id        =~ s/ .*//g;  #$id =~ s/\..*//g;
		#$id        =~ s/ $//g;  #$id =~ s/\..*//g;
		$fas{$id}  = "";
		$ids_orig{$id} = $id0;
	}else{
		$fas{$id} .= $line;
	}
}
close(IDS);
printf STDERR "read %d sequences\n", scalar keys %fas;


# extract sequences and write to screen
foreach $id ( sort { length($fas{$b}) <=> length($fas{$a}) } keys %fas ){
	next if ( $invert == 1 and   exists $ids{$id} );
	next if ( $invert == 0 and ! exists $ids{$id} );
	#printf ">%s\n%s\n", $id, $fas{$id};
	printf ">%s\n%s\n", $ids_orig{$id}, $fas{$id};
}

