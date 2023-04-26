#!/usr/bin/perl

use Data::Dumper;

# read node IDs
my $i = 0;
my $nodes = {};
open( ID, "grep '>' clustered_0_rep_seq.fasta |");
while ( my $line = <ID> ){
	chomp($line);
	$line =~ s/>//g;
	$line =~ s/ *$//g;
	$i++;
	$nodes->{$line} = {};
	$nodes->{$line}->{'id'} = $i;
	$nodes->{$line}->{'li'} = [];
}
close(ID);
#printf "%d\n", scalar keys %{$nodes};

# read clustering
open( CL, "<clustered_1_cluster.tsv" );
while ( my $line = <CL> ){
	chomp($line);
	$line =~ /(.*)\t(.*)/;
	my $id1 = $1;
	my $id2 = $2;
	next if ( $id1 eq $id2 );

#	if ( ! exists $nodes->{$id1} ){
#		printf "'$id1' not found\n";
#	}

	push( @{ $nodes->{$id1}->{'li'} }, $id2 );
}
close(CL);
#printf "%d\n", scalar keys %{$nodes};

# write node info to file
open(  OUTA, ">nodes.csv" );
foreach my $i ( keys %{$nodes} ){
	my $cat = 1;
	   $cat = 2 if ( $i =~ /novel__.*/ );
	printf OUTA "%d,%d,%s\n",  $nodes->{$i}->{'id'}, $cat, $i;
}
close( OUTA );

# write edge info to file
open(  OUTB, ">edges.csv" );
foreach my $i ( keys %{$nodes} ){
	next if ( scalar @{ $nodes->{$i}->{'li'} } == 0 );
	foreach my $j  ( @{ $nodes->{$i}->{'li'} } ){
		printf OUTB "%d,%d,1\n",  $nodes->{$i}->{'id'}, $nodes->{$j}->{'id'};
	}
}
close( OUTB );

