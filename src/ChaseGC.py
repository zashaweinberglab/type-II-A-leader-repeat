import sys
#print ('\n'.join(sys.path))
import numpy as np
import RNA
import string
import argparse
import re

from ChaseLibrary import LoadLeaderRepeatFile,ParseSsToArrayMap,GetPairProbFromProbArray

def GC (seq):
    l=len(seq)
    gcless=re.sub(r'[GCgc]','',seq)
    gc=l-len(gcless)
    return float(gc)/float(l)

parser = argparse.ArgumentParser('quick & dirty G+C content calculations')
parser.add_argument("inFileName", help="input file name (tab-delimited, fields are: name,leaderSeq,repeatSeq)",type=str)
args = parser.parse_args()

leaderRepeatsFileName=args.inFileName

moveNucsFromLeaderToRepeat=0
maxLeaderLen=-1
leaderRepeatList=LoadLeaderRepeatFile(leaderRepeatsFileName,moveNucsFromLeaderToRepeat,maxLeaderLen)

for leaderRepeat in leaderRepeatList:
    (name,leaderSeq,repeatSeq)=leaderRepeat
    print("%s : %g , %g" % (name,GC(leaderSeq+repeatSeq),GC(repeatSeq)))
