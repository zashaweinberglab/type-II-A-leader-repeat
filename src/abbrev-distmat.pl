use strict;

# perl abbrev-distmat.pl II-A-clustered70.distmat II-A-clustered70.phylip-distmat II-A-clustered70.phylip-tree II-A-clustered70-myNames.phylip-tree


my $inFile=shift;
my $outFile=shift;
my $phylipTreeFile=shift;
my $nameTreeFile=shift;

#my $inFile="II-A-clustered70.distmat";
#my $outFile="II-A-clustered70.phylip-distmat";
#my $phylipTreeFile="II-A-clustered70.phylip-tree";
#my $nameTreeFile="II-A-clustered70-myNames.phylip-tree";

my %abbrevToName=();
my $idlen=11;
my $i=0;
if (!open(OUT,">$outFile")) {
    die "cannot open $outFile";
}
if (!open(IN,$inFile)) {
    die "cannot open $inFile";
}
while (<IN>) {
    s/[\r\n]//g;
    my $line=$_;
    if (/^[0-9]+$/) {
	# number of seqs, leave it alone
	print OUT "$line\n";
    }
    else {
	my @fields=split /\s+/,$line;
	my $name=shift @fields;
	my $abbrev="NAME$i";
	my $spaces=" " x ($idlen-length($abbrev));
	my $padi=$abbrev.$spaces;
	print OUT join(" ",($padi,@fields))."\n";
	$abbrevToName{$abbrev}=$name;
	$i++;
    }
}
close(IN);
close(OUT);


if (-e $phylipTreeFile) {
    if (!open(IN,$phylipTreeFile)) {
	die "cannot open $phylipTreeFile";
    }
    if (!open(OUT,">$nameTreeFile")) {
	die "cannot open $nameTreeFile";
    }
    while (<IN>) {
	s/[\r\n]//g;

	# super brute force -- just to find&replace on each line
	for my $abbrev (keys %abbrevToName) {
	    my $name=$abbrevToName{$abbrev};
	    s/$abbrev:/\n$name:/;
	}
	print OUT "$_\n";
    }
    close(OUT);
    close(IN);
}
