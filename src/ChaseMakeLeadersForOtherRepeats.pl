use strict;

use FetchSeqNcbi;
use Stockholm;

my $actualFile=shift;
my $allRepeatFile=shift;

my $isHeaderRow=1;
if (!open(F,$actualFile)) {
    die "cannot open $actualFile";
}
if (!open(OUT,">$allRepeatFile")) {
    die "cannot open $allRepeatFile";
}
while (<F>) {
    s/[\r\n]//g;
    my $line=$_;
    if ($isHeaderRow) {
	$isHeaderRow=0;
	next;
    }
    
    my ($name,$leaderSeq,$repeatSeq)=split /\t/,$line;

    $leaderSeq =~ tr/T/U/;
    $repeatSeq =~ tr/T/U/;

    my $leaderLen=length($leaderSeq);

    my ($seqId,$s,$e);
    if ($name =~ /^([A-Z0-9_]+)-CRT_([0-9]+)_([0-9]+)-([FR])$/) {
	$seqId=$1;
	$s=$2;
	$e=$3;
	my $fr=$4;
	if ($fr eq "R") {
	    ($s,$e)=($e,$s);
	}
    }
    else {
	die "name \"$name\" didn't conform to my regex";
    }

    my $dir=$s<$e ? +1 : -1;
    $s -= $dir*300;
    $e += $dir*300;
    my ($seq,$regionStart,$regionEnd)=FetchSeqNcbi($seqId,$s,$e);

    $seq=uc($seq);
    $seq =~ tr/T/U/;

    my @repeatIndexList=();
    my $startIndex=0;
    while (1) {
	my $foundIndex=index($seq,$repeatSeq,$startIndex);
	if ($foundIndex==-1) {
	    last;
	}
	push @repeatIndexList,$foundIndex;
	
	$startIndex=$foundIndex+1; # get next
    }

    my $n=scalar @repeatIndexList;
    print "$name : found $n repeats, including the leader's repeat\n";
    shift @repeatIndexList; # get rid of first repeat
    for my $repeatIndex (@repeatIndexList) {
	my $thisLeaderSeq=substr($seq,$repeatIndex-$leaderLen,$leaderLen);
	print OUT join("\t",("$name-repeat$repeatIndex",$thisLeaderSeq,$repeatSeq))."\n";
    }
}
close(OUT);
close(F);
