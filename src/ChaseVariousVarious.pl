use strict;
use Getopt::Long;

use ProcessPipeline_Utils;

my $leaderSeqLenListStr="0,100,80";
my $datasets="II-A--training II-C--training"; # can't be in a dir, due to limitations in my interest in fighting with bash
my $datasetsDesc="";
my $fixBug1;
if (!GetOptions(
	 "datasetsFiles=s" => \$datasets,
	 "datasetsDesc=s" => \$datasetsDesc,
	 "leaderSeqLenList=s" => \$leaderSeqLenListStr,
	 "fixBug1" => \$fixBug1,
    )) {
    die "problem with commandline";
}

my @leaderSeqListStr=split /,/,$leaderSeqLenListStr;
for my $leaderSeqLen (@leaderSeqListStr) {
    for my $munge (
	["normal",""],
	["noWithinRepeat","--no-within-repeat"],
	["extendRepeat15-noWithinRepeat","--extend-repeat-by 15 --no-within-repeat"],
	) {
	my ($fileAdd,$flags)=@$munge;
	if ($leaderSeqLen>0) {
	    $fileAdd .= "-leader$leaderSeqLen";
	    $flags .= " --max-leader-len $leaderSeqLen";
	}
	if ($fixBug1) {
	    $flags .= " --fix-bug1";
	}
	if (length($flags)==0) {
	    next; # did this already
	}

	my $touchFile="various-touch--$fileAdd--$datasetsDesc";
	print "touchFile=$touchFile\n";
	::DoIfNotTouchOrForce
	    ($touchFile,
	     sub {
		 my $cmd="~/bliss/code/motifs_2007/rpl/ChaseVarious.sh $fileAdd \"$flags\" \"$datasets\" \"$datasetsDesc\"";
		 print "\n\n\n\n\n\n------------\n\n\n$cmd\n\n\n";
		 my $r=system($cmd);
		 if ($r!=0) {
		     die "problem with $cmd: $? $!";
		 }
	     },
	    0);
    }
}

