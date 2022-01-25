import sys
import csv
import argparse
import os
import re
import random
import logging
from ChaseLibrary import LoadLeaderRepeatFile

# data used to run CD-HIT to cluster the sequences
cdhitWordSizeList=[(0.4,2),(0.5,3),(0.6,4),(0.7,5)]

# function GetCdhitWordSize: to use CD-HIT to cluster sequences, find the appropriate word size for a given percent identity
def GetCdhitWordSize (desiredPercentId):
    desiredWordSize=None
    for cdhitWordSize in cdhitWordSizeList:
        thisPercentId,thisWordSize=cdhitWordSize
        if desiredPercentId>=thisPercentId:
            desiredWordSize=thisWordSize

    if desiredWordSize is None:
        raise AssertionError('that percentId is not possible with cd-hit, or I didn\'t know it was')
    return desiredWordSize

# function RunCdhit:  run CD-HIT to cluster sequences
def RunCdhit (inFastaFileName,outFileName,percentId) :
    wordSize=GetCdhitWordSize(percentId)
    cmd="cd-hit -i %s -o %s -c %f -n %d -d 0 > /dev/null" % (inFastaFileName,outFileName,percentId,wordSize)
    print("Running: %s" % (cmd))
    os.system(cmd)

# function LoadLeaderRepeatList: load leader/repeat pairs
def LoadLeaderRepeatList(inFileName):
    return LoadLeaderRepeatFile(inFileName,0,-1)

# function LoadCdhitClusters: load the clusters in the output of CD-HIT, after it clustered the sequences
def LoadCdhitClusters(clustersFileName):
    clusterList=[]
    currClusterList=[]
    with open(clustersFileName) as clustersFile:
        for line in clustersFile:
            newClusterMatch = re.match(r'^>Cluster', line)
            if newClusterMatch:
                if currClusterList:
                    clusterList.append(currClusterList)
                    currClusterList=[]
            memberMatch = re.match(r'^[0-9]+\s+[0-9]+aa,\s+>(.+)[.][.][.] ',line)
            if memberMatch:
                name=memberMatch.group(1)
                currClusterList.append(name)
            if not newClusterMatch and not memberMatch:
                raise AssertionError('could not figure out what line %s is' % (line))
        if currClusterList:
            clusterList.append(currClusterList)
    return clusterList

# function OutputNewTabFile: create a new .tab file with a subset of leader/repeat seqs.  The subsets will correspond to training, validation, test data
def OutputNewTabFile (dataset,baseFileName,clusterList,nameToSeqs):
    outFileName=baseFileName+"--"+dataset+".tab"
    count=0
    with open(outFileName,'w') as outFile:
        out = csv.writer(outFile, delimiter='\t', quoting=csv.QUOTE_MINIMAL)
        for cluster in clusterList:
            for member in cluster:
                if member not in nameToSeqs:
                    raise AssertionError('cannot find member')
                (leaderSeq,repeatSeq)=nameToSeqs[member]
                out.writerow([member,leaderSeq,repeatSeq])
                count += 1
    print("for %s got %d members" % (dataset,count))

# parse command line & determine parameter values

parser = argparse.ArgumentParser('split dataset based on cd-hit clustering')
parser.add_argument("inFileName", help="input file name of systems to cluster (tab-delimited, fields are: name,leaderSeq,repeatSeq)",type=str)
parser.add_argument("outFileBase", help="base of output files with different partitions",type=str)
parser.add_argument("firstPercentId",help="coarsest percent id to use (e.g 50%).  Each member of one data set (training, validation or test) must be at most this much similar to any member of a different data set.  If you set this number too high (e.g. close to 100%), the training and test sets could potentially be so similar to one another that we're effectively going to test on the training data.  If you set this number too low (e.g. 0%), it might not be possible to partition the data at all.",type=float)
parser.add_argument("secondPercentId",help="finer percent id to use (e.g. 70%).  Each memer of some data set (training, validation or test) must be at most this much similar to any other member of the same data set.  If you set this number too high (e.g. close to 100%), the members of the test data set (and also the other test sets) will potentially be highly correlated with each other.  Techniques like Fisher's Method that assume that the samples are independent could be compromised because the samples are very much not independent.  If you make this number too low (e.g. close to 0%), your data sets will be too small to do any meaningful analysis.",type=float)
parser.add_argument("--force-in-training-file",help="force these seqs to be in training.  These seqs have previously been analyzed, possibly manually by Chase, so we should treat them as contaminated and put them in the training data set.",type=str,default="",dest="inTrainingFile")
parser.add_argument("--train-valid-test-nums",help="number of first-level clusters for training, validation and test sets.  Exact numbers must be given, comma sep, e.g. \"8,4,5\", and they must sum to the actual number of clusters remaining after removing the forced-to-be-in-training members",type=str,dest="trainValidTestNums")
parser.add_argument("--seed",help="random number seed",type=int,default=203040)
args = parser.parse_args()

inFileName=args.inFileName
outFileBase=args.outFileBase
tempFastaFileName=outFileBase + "--temp.fasta"
tempFirstCdhitOutFileName=outFileBase + "--temp-first-cdhit"
tempSecondCdhitOutFileName=outFileBase + "--temp-second-cdhit"

firstPercentId=args.firstPercentId
secondPercentId=args.secondPercentId

if firstPercentId>1 or secondPercentId>1:
    raise AssertionError("the percent id values must be expressed as a fraction between 0 and 1")

trainValidTestNums=args.trainValidTestNums
if not args.trainValidTestNums:
    trainValidTestNums="0,0,0" # doesn't matter
(numClustersTraining,numClustersValidation,numClustersTest)=[int(x) for x in trainValidTestNums.split(",")] # this is the target number of members in each of the data sets.  you need to run the program once to see what numbers work.

random.seed(args.seed)
logging.basicConfig(filename='clustering.log',level=logging.DEBUG)

# optionally load a list of systems that should be forced into the training data set, because they've been analyzed manually already.  we match systems based on their names.
forceTrainingSet={}
if len(args.inTrainingFile)>0:
    lrList=LoadLeaderRepeatList(args.inTrainingFile)
    for row in lrList:
        (name,leaderSeq,repeatSeq)=row
        forceTrainingSet[name]=1

# load the systems we'll partition
leaderRepeatList=LoadLeaderRepeatList(inFileName)
print("read %d systems from file %s" % (len(leaderRepeatList),inFileName))


# first, make .fasta file of the leader seqs (not repeats) for cd-hit to analyze.  Thus, cluster based on leader seqs
with open(tempFastaFileName,'w') as tempFastaFile:
    for leaderRepeat in leaderRepeatList:
        (name,leaderSeq,repeatSeq)=leaderRepeat
        print ("qber %s,%s,%s",(name,leaderSeq,repeatSeq))
        tempFastaFile.write(">{}\n".format(name))
        tempFastaFile.write(leaderSeq+"\n")

# run the clustering by the two percent identity threaholds
RunCdhit(tempFastaFileName,tempFirstCdhitOutFileName,firstPercentId)
RunCdhit(tempFastaFileName,tempSecondCdhitOutFileName,secondPercentId)

tempFirstClustersFileName=tempFirstCdhitOutFileName+".clstr"
tempSecondClustersFileName=tempSecondCdhitOutFileName+".clstr"

firstClusterList=LoadCdhitClusters(tempFirstClustersFileName)
print("got %d clusters at percent id %g" % (len(firstClusterList),firstPercentId))
secondClusterList=LoadCdhitClusters(tempSecondClustersFileName)
print("got %d clusters at percent id %g" % (len(secondClusterList),secondPercentId))

if False:
    for c in firstClusterList:
        print("CLUSTER:")
        for m in c:
            print("\t%s" % (m))

# pick representative of each cluster in second clustering -- we'll output these
forcedMembersThatCantBeTakenList=[]
selectedMemberSet={}
for cluster in secondClusterList:
    forcedMember=None # must pick a forced member, if any
    selectedMember=None
    for member in cluster:
        if not selectedMember:
            selectedMember=member
        if member in forceTrainingSet:
            if forcedMember:
                forcedMembersThatCantBeTakenList.append(member)
            else :
                forcedMember=member
    if forcedMember:
        selectedMember=forcedMember # take the forced member, if we have one
    selectedMemberSet[selectedMember]=1

print("the following forced-in-training members could not be taken because they were too similar to others: ")
print(forcedMembersThatCantBeTakenList)

# remove everything in the first list, except for the selected members
selectedClusterList=[]
for cluster in firstClusterList:
    newCluster=[]
    for member in cluster:
        if member in selectedMemberSet:
            newCluster.append(member)
    selectedClusterList.append(newCluster)

# start to work out the final cluster lists

trainingClusterList=[]
validationClusterList=[]
testClusterList=[]

# first take care of things forced into training

numForcedIntoTraining=0
selectedClusterListNoForced=[]
for cluster in selectedClusterList:
    forceTraining=False
    for member in cluster:
        if member in forceTrainingSet:
            forceTraining=True
    if forceTraining:
        #logging.debug('forced into training cluster:')
        #for member in cluster:
        #    logging.debug(member)
        trainingClusterList.append(cluster)
        numForcedIntoTraining += 1
    else:
        selectedClusterListNoForced.append(cluster)

print("forced %d clusters into training" % (numForcedIntoTraining))
print("%d clusters remain" % (len(selectedClusterListNoForced)))
print("selecting %d,%d,%d for training,validation,test" % (numClustersTraining,numClustersValidation,numClustersTest))

random.shuffle(selectedClusterListNoForced)

if numClustersTraining + numClustersValidation + numClustersTest != len(selectedClusterListNoForced):
    raise AssertionError('the numbers for training, validation, test do not add up to the number of clusters remaining')

trainingClusterList.extend(selectedClusterListNoForced[0:numClustersTraining])
validationClusterList.extend(selectedClusterListNoForced[numClustersTraining:numClustersTraining+numClustersValidation])
testClusterList.extend(selectedClusterListNoForced[numClustersTraining+numClustersValidation:numClustersTraining+numClustersValidation+numClustersTest])

# now we can output them into files

# populate dict to lookup the seq info
leaderRepeat_nameToSeqs={}
for leaderRepeat in leaderRepeatList:
    (name,leaderSeq,repeatSeq)=leaderRepeat
    leaderRepeat_nameToSeqs[name]=(leaderSeq,repeatSeq)

OutputNewTabFile("training",outFileBase,trainingClusterList,leaderRepeat_nameToSeqs)
OutputNewTabFile("validation",outFileBase,validationClusterList,leaderRepeat_nameToSeqs)
OutputNewTabFile("test",outFileBase,testClusterList,leaderRepeat_nameToSeqs)


