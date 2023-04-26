#!/usr/bin/perl

#
# parameters
#
my $infile = "DNA/spades/scaffolds.fasta";
my $wdir   = "DNA/hmmsearch";
my $parts  = 1000;


#
# run
#

# prepare directories and files
mkdir( "$wdir" )       if ( ! -d "$wdir" );
mkdir( "$wdir/split" ) if ( ! -d "$wdir/split") ;

# translate ORFs
`getorf -sequence $infile -outseq $wdir/getorf.fasta -minsize 300` if ( ! -e "$wdir/getorf.fasta" );
my %seqs = ();
readFasta( "$wdir/getorf.fasta", \%seqs );

# calculate sequence number per part
my $snum  = scalar keys %seqs;
   $parts = $snum if ( $snum < $parts );
my $SperP = int( $snum / $parts + 0.5 ); # round to closest integer
print "will split Fasta file into $parts parts with $SperP sequences each\n";

# split
`mv $wdir/getorf.fasta $wdir/split/getorf.fasta`;
splitFasta( "$wdir/split/getorf.fasta", \%seqs, $SperP, $parts, "$wdir/split_files.txt" );
`mv $wdir/split/getorf.fasta $wdir/getorf.fasta`;


# read Fasta file into memory
sub readFasta{
        my $cFile = shift;
        my $cntgs = shift;
        my $contigID = "";
        open( CNTG, "<$cFile" ) or die ( "Can't open file '$cFile': $!\n" );
        while ( my $line = <CNTG> ){
                chomp($line);  next if $line eq "";
                if ( $line =~ />(.*)/ ){
                        $contigID = $1;
                        $contigID =~ s/ .*//g;
                        $$cntgs{ $contigID } = {};
                        $$cntgs{ $contigID }->{'seq'}  = "";
                }else{
                        $$cntgs{ $contigID }->{'seq'} .= $line;
                }
        }
        close(CNTG);
}


# split Fasta file into parts
sub splitFasta{
        my $fasin = shift;
        my $cntgs = shift;
        my $CperP = shift;
	my $partN = shift;
	my $partF = shift;
        my ($ci, $pi) = (0, 0);
        my $pfas = "";
	open(PART,">$partF") or die("Can't open file '$partF': $!\n");
        foreach my $cid ( keys %{$cntgs} ){
                $ci++;
                if ( $ci > $CperP and $pi < ($partN-1) ){
                        $ci = 1;
                        $pi++;
                        close( PFAS ) or die ("Can't close file '$pfas': $!\n");
                }
                if ( $ci == 1 ){
                        $pfas = $fasin.".".$pi;
                        open( PFAS, ">$pfas" ) or die ("Can't open file '$pfas': $!\n");
			printf PART "%d %s\n", $pi, $pfas;
                }
                printf PFAS ">%s\n%s\n", $cid, $$cntgs{$cid}->{'seq'};
        }
        close( PFAS ) or die ("Can't close file '$pfas': $!\n");
        close( PART ) or die ("Can't close file '$partF': $!\n");
} # end splitFasta


