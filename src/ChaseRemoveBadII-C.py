import sys
#print ('\n'.join(sys.path))
import numpy as np
import RNA
import csv
import string
import time
import operator
import argparse
import random
import re

from ChaseLibrary import LoadLeaderRepeatFile,ParseSsToArrayMap,GetPairProbFromProbArray

srcTrainingFile="II-C_Leader_firsrepeat_4May.tab.csv" # just for checking if we got everything.  it's not really the training, since we eliminate redundant seqs
srcCombinedFile="train-test-II-C/II-C--combined-test.tab"
destCombinedFile="train-test-II-C/II-C--combined-test-stricter.tab"
#srcCombinedFile="II-C_Leader_firsrepeat_4May.tab.csv"
#destCombinedFile="II-C_Leader_firsrepeat_OmerRejects.tab.csv"
omerRejectFile="Genomes_were_II-C_And_now_become_II-A.tab_Final_26Nov.tab"

moveNucsFromLeaderToRepeat=0
maxLeaderLen=-1
srcList=LoadLeaderRepeatFile(srcCombinedFile,moveNucsFromLeaderToRepeat,maxLeaderLen)
trainingList=LoadLeaderRepeatFile(srcTrainingFile,moveNucsFromLeaderToRepeat,maxLeaderLen)

rejectList=[]
with open(omerRejectFile) as srcFile:
    reader=csv.reader(srcFile, delimiter='\t', quotechar='\"')
    for row in reader:
        name=row[0]
        if name!='ID':
            rejectList.append(name)

rejectCount=0
acceptCount=0
rejectSet={x for x in rejectList}
usedRejectSet=set()

destList=[]
for x in srcList:
    (name,leaderSeq,repeatSeq)=x
    if name in rejectSet:
        rejectCount += 1
        usedRejectSet.add(name)
    else:
        acceptCount += 1
        destList.append(x)

print("accept=%d , reject=%d" % (acceptCount,rejectCount))

# check to see which to-be-rejected names are in the training set
for x in trainingList:
    (name,leaderSeq,repeatSeq)=x
    if name in rejectSet:
        usedRejectSet.add(name)
for x in rejectSet:
    if x not in usedRejectSet:
        print("%s is in reject set, but not found" % (x))
        
with open(destCombinedFile,'w',newline='') as csvfile:
    outFile = csv.writer(csvfile, delimiter='\t', quoting=csv.QUOTE_MINIMAL)
    for x in destList:
        outFile.writerow(x)
