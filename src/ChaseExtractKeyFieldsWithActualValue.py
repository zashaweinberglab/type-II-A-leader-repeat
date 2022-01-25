import sys
import os
#print ('\n'.join(sys.path))
import numpy as np
import RNA
import string
import argparse
import re
import csv

from ChaseLibrary import LoadLeaderRepeatFile,ParseSsToArrayMap,GetPairProbFromProbArray,ParseSsToPairList,ChaseFindLeaderAndRepeatMfe,ChaseCalcBppAsDigitStr,ChaseDrawLeaderRepeat

def GC (seq):
    l=len(seq)
    gcless=re.sub(r'[GCgc]','',seq)
    gc=l-len(gcless)
    return float(gc)/float(l)

parser = argparse.ArgumentParser('extract the most-important fields from ChaseFolding output')
parser.add_argument("leaderRepeatsFileName", help="input file name (tab-delimited, fields are: name,leaderSeq,repeatSeq)",type=str)
parser.add_argument("foldingFileName",help="output of ChaseFolding",type=str)
parser.add_argument("outFileName",help="output",type=str)
parser.add_argument("--showMfe",help="output MFE structure",action="store_true")
parser.add_argument("--showBpp",help="show base-pair probabilities",action="store_true")
parser.add_argument("--sortActualBpp",help="sort results about the actual value of base pairing probs",action="store_true")
parser.add_argument("--max-leader-len",help="if a leader is longer than this amount, chop nucs from the 5' end to truncate it to this length",type=int,default=-1,dest="maxLeaderLen")
parser.add_argument("--draw-dir",help="draw everything with R2R using oneseq, putting PDFs into this directory.  r2r must be in the path.",type=str,dest="drawDir")
args = parser.parse_args()

if args.drawDir:
    os.makedirs(args.drawDir,exist_ok=True)
    args.showMfe=True # we need these data, so turn this on
    args.showBpp=True

statsShortlistList=['avgPairProbs','helix-0-8','helix-1-11']

moveNucsFromLeaderToRepeat=0
maxLeaderLen=args.maxLeaderLen
leaderRepeatList=LoadLeaderRepeatFile(args.leaderRepeatsFileName,moveNucsFromLeaderToRepeat,maxLeaderLen)

nameToRepeatGc={}
nameToLeaderRepeatGc={}
nameToValues={}
nameToActualValues={}
nameToLeaderRepeatSeq={}

for leaderRepeat in leaderRepeatList:
    (name,leaderSeq,repeatSeq)=leaderRepeat
    nameToLeaderRepeatGc[name]=GC(leaderSeq+repeatSeq)
    nameToRepeatGc[name]=GC(repeatSeq)
    nameToLeaderRepeatSeq[name]=[leaderSeq,repeatSeq]

md = RNA.md()

with open(args.foldingFileName,'r') as inputFile,open(args.outFileName,'w') as outfile:
    reader=csv.reader(inputFile,delimiter='\t')
    statNameList=next(reader, None) # read stat names
    statNameList=statNameList[1:] # but actually, we ignore the first column

    statNameToCol={}
    for col in range(len(statNameList)):
        statNameToCol[statNameList[col]]=col
    statsShortlistColList=[statNameToCol[statName] for statName in statsShortlistList]

    for row in reader:
        rowName=row[0]
        if re.search(r'params=',rowName):
            continue # ignore these lines

        statList=[float(x) for x in row[1:]]
        isActualValue=re.search(r'actual-value$',rowName) is not None
        rowName=re.sub(r'actual-value$','',rowName)
        if isActualValue:
            nameToActualValues[rowName]=statList
        else:
            nameToValues[rowName]=statList

    if args.sortActualBpp:
        nameList=sorted(nameToValues,key=lambda x : nameToActualValues[x][statNameToCol['avgPairProbs']])
    else:
        nameList=sorted(nameToValues,key=lambda x : nameToValues[x][statNameToCol['helix-0-7']] + nameToValues[x][statNameToCol['avgPairProbs']])

    writer=csv.writer(outfile,delimiter='\t',quoting=csv.QUOTE_MINIMAL)
    writer.writerow(['name'] + ['p-value('+x+')' for  x in statsShortlistList] + ['actualValue('+x+')' for x in statsShortlistList] + ['leader-repeat-GC','repeat-GC'])
    for name in nameList:
        values=nameToValues[name]
        actualValues=nameToActualValues[name]
        leaderRepeatGc=nameToLeaderRepeatGc[name]
        repeatGc=nameToRepeatGc[name]
    
        writer.writerow([name] + [values[col] for col in statsShortlistColList]  + [actualValues[col] for col in statsShortlistColList] + [leaderRepeatGc,repeatGc])

        if args.showMfe or args.showBpp:
            (leaderSeq,repeatSeq)=nameToLeaderRepeatSeq[name]
            #print("%s / %s" % (leaderSeq,repeatSeq))
            seq=leaderSeq+repeatSeq

        if args.showMfe:
            leaderSs,repeatSs=ChaseFindLeaderAndRepeatMfe(leaderSeq,repeatSeq,md)
            writer.writerow(['',leaderSeq,repeatSeq])
            writer.writerow(['',leaderSs,repeatSs])

        if args.showBpp:
            digitStr=ChaseCalcBppAsDigitStr(leaderSeq,repeatSeq,md)
            writer.writerow(['',digitStr[0:len(leaderSeq)],digitStr[len(leaderSeq):]])

        if args.drawDir:
            drawPdfFileName=args.drawDir+"/"+name+".pdf" # hopefully 'name' is okay for a file
            drawStoFileName=args.drawDir+"/"+name+".sto"
            drawMetaFileName=args.drawDir+"/"+name+".r2r_meta"
            ChaseDrawLeaderRepeat(name,leaderSeq,repeatSeq,leaderSs,repeatSs,digitStr,drawPdfFileName,drawStoFileName,drawMetaFileName)
