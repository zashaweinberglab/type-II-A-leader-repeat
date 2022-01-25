import sys
#print ('\n'.join(sys.path))
import numpy as np
import RNA # ViennaRNA library
import csv
import string
import time
import operator
import argparse
import random
import re
from multiprocessing import Process,Queue

# import libraries I wrote
from altschulEriksonDinuclShuffle import dinuclShuffleUnnormalizedString
from ChaseLibrary import LoadLeaderRepeatFile,ParseSsToArrayMap,GetPairProbFromProbArray

jobQueue=Queue() # for parallel processing with multiple processes (library 'multiprocessing')
resultQueue=Queue()

d=False # whether to print debugging messages


if sys.version_info < (3,0,0):
    print(__file__ + ' requires Python 3, while Python ' + str(sys.version[0] + ' was detected. Terminating. '))
    sys.exit(1)

# function CalcListOfPairsAndNumBadForSs: analyze the number of base pairs that are consecutive or allowing bulges/internal loops.
# purpose: given a secondary structure, find the distances between all base pairs in the repeat.  get these distances in terms of the NumBad (e.g. 1-nt bulge --> NumBad==1, contiguous base pairs --> NumBad==0). return these distances as a list.  later we'll process them to find the number of base pairs we can get, with a maximum total NumBad
# input:
#   ss = secondary structure.  
#   firstRepeatSeqPos,lastRepeatSeqPos = location of repeat.  NOTE: the code assumes that the repeat is to the right of everything interesting.  This will require some modification for type V-A systems
#   doAllowWithinRepeat: if false, then completely ignore base pairs that are contained within the repeat
# output:
#   list of NumBad-style distances between base pairs.
# problem: doesn't distinguish multistem junctions: however, the maxNumBad should take care of it
# problem: should I accept pairs within repeatSeq or not?  It's a parameter.  If we allow within-repeat stems, we should still require only looking at right pairs, to avoid counting left and then right pairs.
# problem: we should perhaps impose a maximum bulge/internal loop length for things to be considered together.
# problem: must propagate best values.  numBad-->numPairs should be the most number of pairs we can get with **UP TO** numBad skips.
def CalcListOfPairsAndNumBadForSs(ss,firstRepeatSeqPos,lastRepeatSeqPos,doAllowWithinRepeat):
    pairAndNumBadList=[]

    pairArrayMap=ParseSsToArrayMap(ss)

    prevR=-1
    
    for currR in range(firstRepeatSeqPos,lastRepeatSeqPos):

        currL=pairArrayMap[currR]
        if currL>=0 and currL<currR: # valid base pair for us -- it's a base pair, and we're on the right
            if doAllowWithinRepeat or currL<firstRepeatSeqPos: # make sure it's not within the repeat, or we're taking that
                thisNumBad=0 # if no previous pair (i.e., this is the first one), then we can just take it
                if prevR!=-1:
                    thisNumBad=max(abs(prevL-currL),abs(prevR-currR))-1 # -1 because a distance of 1 corresponds to consecutive base pairs.  the 'max' function is to take care of bulges and mismatches/internal loops, which I kind of put on an equivalence
                pairAndNumBadList.append(thisNumBad) # we could include currL,currR, but there's no clear need

                prevL=currL
                prevR=currR

    #print(pairAndNumBadList)
    return pairAndNumBadList

# function AddCountsForNumBadAndNumPairs: further analyze what kind of helices (number of base pairs, allowed number of bulges/internal loop nucleotides)
def AddCountsForNumBadAndNumPairs (countsByNumBadAndNumPairs,maxNumBad,maxNumPairs,ss,firstRepeatSeqPos,lastRepeatSeqPos,doAllowWithinRepeat,fixBug1):

    pairAndNumBadList=CalcListOfPairsAndNumBadForSs(ss,firstRepeatSeqPos,lastRepeatSeqPos,doAllowWithinRepeat)

    #print("LIST")
    #print(pairAndNumBadList)

    numPosInList=len(pairAndNumBadList) # we don't have to worry about actual nucleotide positions within repeatSeq, just in the pairs that are there, and the numBad between them

    # dyn programming -- I think this is probably fastest given that maxNumBad should be fairly small (like not more than 3), and I think it's easier to program. it seems to fit into a nice series of relatively simple problems.  Obviously it'd be easier to fix a length & numBad, and just compute that
    mostPairsAtNumBadAndPos = np.zeros((maxNumBad+1,numPosInList),dtype=int) # This is the dyn prog table.  first dimension is the NumBad we can allow, and the second dimension is the index (0-based) into the list pairAndNumBadList. BTW, we only need current and previous column, but who cares

    listPos=0
    for thisPairInfo in pairAndNumBadList: # go through each element in the list, from index=0 to the end
        numBadAtPos=thisPairInfo
    
        for numBad in range(0,maxNumBad+1): # start with NumBad==0, then increase
            thisMostPairs=1 # can always just take this pair at numBad=0, or can do more bad even if it wasn't necessary
            if numBad>0:
                thisMostPairs=max(thisMostPairs,mostPairsAtNumBadAndPos[numBad-1,listPos]) # can always add a badness, even if unnecessary, except if numBad==0 (obviously can't be negative)

            if numBad>=numBadAtPos and listPos>0: # can we extend from the previous pair?  doesn't work if we're too bad, or if we're at the first list pos (in which case there's nowhere to go)
                thisMostPairs=max(thisMostPairs,mostPairsAtNumBadAndPos[numBad-numBadAtPos,listPos-1] + 1) # +1: we get another pair
                
            mostPairsAtNumBadAndPos[numBad,listPos]=thisMostPairs

        listPos += 1

    #print("dyn tab")
    #print(mostPairsAtNumBadAndPos)

    # now we go through for each numBad, find the mostPairs, and we can set everything from [0..mostPairs] pairs
    for numBad in range(0,maxNumBad+1):
        mostPairs=1
        for listPos in range(0,numPosInList):
            mostPairs=max(mostPairs,mostPairsAtNumBadAndPos[numBad,listPos])
            #print("numBad,listPos=%d,%d , mostPairs=%d (table=%d)" % (numBad,listPos,mostPairs,mostPairsAtNumBadAndPos[numBad,listPos]))

        if fixBug1:
            for numPairs in range(0,min(mostPairs,maxNumPairs)+1): # maxNumPairs is just to make sure we don't exceed array bounds, although in practice it should be easy.  +1: 'range' uses a half-open interval
                countsByNumBadAndNumPairs[numBad,numPairs] += 1
        else:
            for numPairs in range(0,min(mostPairs,maxNumPairs+1)): # buggy original code: missing +1
                countsByNumBadAndNumPairs[numBad,numPairs] += 1
    #print(countsByNumBadAndNumPairs)
    #print(countsByNumBadAndNumPairs[0,7])

# for given leader,repeat sequences, calculate all the statistics we're interested in.  the sequences could be natural or randomized
def CalcStats (leaderSeq,repeatSeq,numBoltzmannSamples,maxNumBad,minNumPairs,maxNumPairs,doAllowWithinRepeat,dumpBoltzmannActual,fixBug1) :

    seq=leaderSeq + repeatSeq

    statMap={}

    # folding from ViennaRNA
    md = RNA.md()
    #md.temperature=37
    fc = RNA.fold_compound(seq,md)
    (propensity, ensemble_energy) = fc.pf()
    basepair_probs = fc.bpp()

    pairProbArray=GetPairProbFromProbArray(basepair_probs,seq)

    firstRepeatSeqPos=len(leaderSeq)
    lastRepeatSeqPos=len(leaderSeq)+len(repeatSeq)
    
    # calc 
    logProbFullyUnpaired=1
    sumOfPairProbs=0
    for i in range(firstRepeatSeqPos,lastRepeatSeqPos):
        p=pairProbArray[i]
        if p<0 or p>1:
            raise AssertionError('p should be a probability')
        unpairProb=1-p
        logProbFullyUnpaired *= np.log(unpairProb) # I hope there aren't any zeroes.  I have a bad feeling about this.
        sumOfPairProbs += p
        
    avgPairProbs=sumOfPairProbs / (lastRepeatSeqPos-firstRepeatSeqPos)
    logGeomMeanFullyUnpaired=logProbFullyUnpaired / (lastRepeatSeqPos-firstRepeatSeqPos)

    statMap['avgPairProbs']=avgPairProbs
    statMap['logGeomMeanFullyUnpaired']=logGeomMeanFullyUnpaired

    # check out the helix stats with Boltzmann samples
    md = RNA.md()
    # activate unique multibranch loop decomposition
    md.uniq_ML = 1
    fc = RNA.fold_compound(seq, md)
    (mfe_ss, mfe) = fc.mfe()
    # rescale Boltzmann factors according to MFE
    fc.exp_params_rescale(mfe)
    # compute partition function to fill DP matrices
    fc.pf()
    

    countsByNumBadAndNumPairs=np.zeros((maxNumBad+1,maxNumPairs+1),dtype=int)

    if dumpBoltzmannActual or d:
        print("%s%s" % (str.lower(leaderSeq),str.upper(repeatSeq)))

    # estimate probabilities of structural observations (e.g. helix with 8 consecutive base pairs) by sampling from the Boltzmann distribution.
    for ss in fc.pbacktrack(numBoltzmannSamples):

        # for testing: countsByNumBadAndNumPairs=np.zeros((maxNumBad+1,maxNumPairs+1),dtype=int)
        AddCountsForNumBadAndNumPairs (countsByNumBadAndNumPairs,maxNumBad,maxNumPairs,ss,firstRepeatSeqPos,lastRepeatSeqPos,doAllowWithinRepeat,fixBug1)

        if dumpBoltzmannActual:
            leaderSs=ss[0:len(leaderSeq)].translate(str.maketrans('()','<>'))
            repeatSs=ss[firstRepeatSeqPos:lastRepeatSeqPos].translate(str.maketrans('()','[]'))
            contextySs=leaderSs + repeatSs
            print(contextySs)
        
        if d: # debugging dumps of the structures and numBad
            leaderSs=ss[0:len(leaderSeq)].translate(str.maketrans('()','<>'))
            repeatSs=ss[firstRepeatSeqPos:lastRepeatSeqPos].translate(str.maketrans('()','[]'))
            contextySs=leaderSs + repeatSs
            print(contextySs)
            print(countsByNumBadAndNumPairs)

    for numBad in range(0,maxNumBad+1):
        for numPairs in range(minNumPairs,maxNumPairs+1):
            f = countsByNumBadAndNumPairs[numBad,numPairs] / numBoltzmannSamples
            statMap['helix-{}-{}'.format(numBad,numPairs)]=f
            
    return statMap

# convenience to transform a job in Python's multi-process stuff into a function call
def CalcStatsViaWork (work):
    (leaderSeq,repeatSeq,numBoltzmannSamples,maxNumBad,minNumPairs,maxNumPairs,doAllowWithinRepeat,dumpBoltzmannActual,fixBug1)=work
    return CalcStats(leaderSeq,repeatSeq,numBoltzmannSamples,maxNumBad,minNumPairs,maxNumPairs,doAllowWithinRepeat,dumpBoltzmannActual,fixBug1)

# parallel processing
def WorkerFunction ():
    random.seed() # try to make sure that each worker is using different numbers
    while True:
        pill=jobQueue.get()
        (isPoisonPill,work)=pill
        if isPoisonPill:
            break
        result=CalcStatsViaWork(work)
        resultQueue.put(result)
    # fall through to our death

# parallel processing
def MakeCalcStatsWork(leaderSeq,repeatSeq,numBoltzmannSamples,maxNumBad,minNumPairs,maxNumPairs,doAllowWithinRepeat,dumpBoltzmannActual,fixBug1):
    
    negSampleLeaderSeq=leaderSeq
    negSampleRepeatSeq=repeatSeq
    if doShuffleRepeatAlso:
        # shuffle leader & repeat together, then extract same-size sequences
        combinedSeq=dinuclShuffleUnnormalizedString(leaderSeq+repeatSeq)
        negSampleLeaderSeq=combinedSeq[0:len(leaderSeq)]
        negSampleRepeatSeq=combinedSeq[len(leaderSeq):len(combinedSeq)]
    else:
        # just shuffle the leader seq
        negSampleLeaderSeq=dinuclShuffleUnnormalizedString(leaderSeq)

    return (negSampleLeaderSeq,negSampleRepeatSeq,numBoltzmannSamples,maxNumBad,minNumPairs,maxNumPairs,doAllowWithinRepeat,dumpBoltzmannActual,fixBug1)

def WriteRowTextAndHtml (outFile,htmlFile,fieldList):
    outFile.writerow(fieldList)
    htmlFile.write("<tr>"+"".join(["<td>"+str(x)+"</td>" for x in fieldList])+"</tr>")

def ProcessShuffledResults(statMapSample,outFile,statActualValueList,countVector):
    statSampleValueList=[value for key,value in sorted(statMapSample.items())]
            
    if doPrintAllValues:
        WriteRowTextAndHtml(outFile,htmlFile,[name+'sample-value']+statSampleValueList)

    addVector=[-int(sampleValue>=actualValue) for sampleValue,actualValue in zip(statSampleValueList,statActualValueList)]
    countVector -= addVector # hacking around the annoying fact that the '+' symbol means append, by negating and using subtraction.  this also means that we can set the caller's variable with Python's call-by-pointer, which wouldn't work if we assigned to countVector
#    countVector = list(map(operator.add,countVector,addVector)) # somewhat annoying that + is overloaded as append, even though there's already append.  I could also use numpy arrays.
#print(countVector)



# main function

random.seed()

# the following code is mostly about setting values of configurable parameters based on the command line

maxNumBad=2
minNumPairs=5
maxNumPairs=20
doPrintAllValues=False

params=" ".join(sys.argv)

parser = argparse.ArgumentParser('folding stats for type II-A CRISPR leader/repeat seqs')
parser.add_argument("inFileName", help="input file name (tab-delimited, fields are: name,leaderSeq,repeatSeq)",type=str)
parser.add_argument("-o",dest="outFile",help="output .tab file name",type=str,default="out.tab")
parser.add_argument("--num-sample",dest="numSample",help="num background samples",type=int,default=100)
parser.add_argument("--num-boltzmann",dest="numBoltzmann",help="num Boltzmann samples to estimate helix probabilies",type=int,default=100)
parser.add_argument("--print-actual-value",dest="printActualValue",help="print the actual values of the stats, not just their p-values",action="store_true")
parser.add_argument("--no-within-repeat",dest="noWithinRepeat",help="ignore pairs contained within repeat seq",action="store_true")
parser.add_argument("--shuffle-repeat",dest="shuffleRepeat",help="shuffle repeat along with leader for null model",action="store_true")
parser.add_argument("--test-dinuc", help="test dinuc shuffling stuff",action="store_true",dest="testDinuc")
parser.add_argument("--dump-actual-boltzmann",help="dump boltzmann structures for the actual sequence (as opposed to the shuffled ones).  goes to stdout.",action="store_true",dest="dumpActualBoltzmann")
parser.add_argument("--cpu",help="run multithreaded with this many CPUs",type=int,default=0)
parser.add_argument("--extend-repeat-by",help="move this many nucs on the 3' end of the leader into the 5' end of repeat, to look for hairpins that are just beyond the repeat",type=int,default=0,dest="moveNucsFromLeaderToRepeat")
parser.add_argument("--max-leader-len",help="if a leader is longer than this amount, chop nucs from the 5' end to truncate it to this length",type=int,default=-1,dest="maxLeaderLen")
parser.add_argument("--fix-bug1",help="fix bug where I didn't add 1 to the range, so it doesn't necessarily work if there's exactly the right number of base pairs",action="store_true",dest="fixBug1")
args = parser.parse_args()

if not args.fixBug1:
    raise SystemExit("Are you sure you don't want use --fix-bug1?")

if args.testDinuc:
    # test the dinuc-aware shuffling code for randomized sequences
    seq="AGUCUACUGAUCGAUGAUGCAUGCAGACU"
    list=[seq]
    for i in range(0,10):
        list.append(dinuclShuffleUnnormalizedString(seq))
    for seq in list:
        dinucToCount={}
        for i in range(0,len(seq)-1):
            dinuc=seq[i:i+2]
            dinucToCount[dinuc] = dinucToCount.get(dinuc,0) + 1
        print(seq)
        print(sorted(dinucToCount.items()))
    print("--test-dinuc : done")
    sys.exit(0)

leaderRepeatsFileName=args.inFileName
numShuffledLeaderSeqs=args.numSample
numBoltzmannSamples=args.numBoltzmann
doPrintActualValues=args.printActualValue
doAllowWithinRepeat=not args.noWithinRepeat
doShuffleRepeatAlso=args.shuffleRepeat
outFileName=args.outFile
htmlFileName=outFileName+".html"
dumpBoltzmannActual=args.dumpActualBoltzmann
moveNucsFromLeaderToRepeat=args.moveNucsFromLeaderToRepeat
maxLeaderLen=args.maxLeaderLen

useThreads=args.cpu>0
numCpus=args.cpu

# numShuffled=100, numBoltzmann=100, 4 seqs, 1 CPU took 24 secs

# for parallel processing: set up worker processes (called "threads" because I was originally planning to use Python's multithreading, before I realized that the standard Python implementation doesn't allow more than one thread to run at one time)
threadList=[]
if useThreads:
    if numCpus==1:
        raise AssertionError('--cpu 1 means zero workers, so no work will get done')
    for i in range(0,numCpus):
        #t = threading.Thread(target=WorkerFunction)
        t=Process(target=WorkerFunction)
        t.start()
        threadList.append(t)

# load input file
leaderRepeatList=LoadLeaderRepeatFile(leaderRepeatsFileName,moveNucsFromLeaderToRepeat,maxLeaderLen)

# open output files and start processing
startTime=time.time()
with open(outFileName, 'w', newline='') as csvfile,open(htmlFileName,'w') as htmlFile:
    outFile = csv.writer(csvfile, delimiter='\t', quoting=csv.QUOTE_MINIMAL)
    htmlFile.write("<html><body><table>\n")

    didHeaderRow=False # did we write out the column headers in the output file yet?

    # process each leader-repeat pair
    for leaderRepeat in leaderRepeatList:
        (name,leaderSeq,repeatSeq)=leaderRepeat

        print("processing %s (lens: %d,%d)" % (name,len(leaderSeq),len(repeatSeq)))

        # calculate various statistics on the real biological sequences
        # even with multiple threads, just call this directly -- it's much easier that way
        #(actual_fullyUnpairProb,actual_unpairProbSum,actual_
        statMapActual=CalcStats(leaderSeq,repeatSeq,numBoltzmannSamples,maxNumBad,minNumPairs,maxNumPairs,doAllowWithinRepeat,dumpBoltzmannActual,args.fixBug1)

        if dumpBoltzmannActual:
            continue;

        if not didHeaderRow:
            statKeyList=['name']+sorted(statMapActual.keys())
            WriteRowTextAndHtml(outFile,htmlFile,statKeyList)
            didHeaderRow=True

            
        statActualValueList=[value for key,value in sorted(statMapActual.items())]
        if doPrintAllValues or doPrintActualValues:
            WriteRowTextAndHtml(outFile,htmlFile,[name+'actual-value']+statActualValueList)

        countVector=np.zeros((len(statActualValueList)))

        # calculate statistics on randomized sequences in order to get p-values.  This uses parallel processing, since it's the most time-intensive step by far
        dumpBoltzmannActual=False # for the negative samples, this must be false
        if useThreads:
            for sample in range(0,numShuffledLeaderSeqs):

                work=MakeCalcStatsWork(leaderSeq,repeatSeq,numBoltzmannSamples,maxNumBad,minNumPairs,maxNumPairs,doAllowWithinRepeat,dumpBoltzmannActual,args.fixBug1)
                isPoisonPill=False
                jobQueue.put((isPoisonPill,work))

            for sample in range(0,numShuffledLeaderSeqs):
                statMapSample=resultQueue.get()
                ProcessShuffledResults(statMapSample,outFile,statActualValueList,countVector)
                
        else:
            for sample in range(0,numShuffledLeaderSeqs):

                work=MakeCalcStatsWork(leaderSeq,repeatSeq,numBoltzmannSamples,maxNumBad,minNumPairs,maxNumPairs,doAllowWithinRepeat,dumpBoltzmannActual,args.fixBug1)
                statMapSample=CalcStatsViaWork(work)
                ProcessShuffledResults(statMapSample,outFile,statActualValueList,countVector)

        empiricalPvalueVector = [c / numShuffledLeaderSeqs for c in countVector]
        WriteRowTextAndHtml(outFile,htmlFile,[name]+empiricalPvalueVector)

    if dumpBoltzmannActual:
        print("--dump-actual-boltzmann: early exit, skip negative samples")

    outFile.writerow(["params=%s\n" % (params)])
    htmlFile.write("</table>\n")
    htmlFile.write("<p>params=%s</p>\n" % (params))
    htmlFile.write("</body></html>\n")

endTime=time.time()
elapsedTime=endTime-startTime
print("took %s secs" % (elapsedTime))


print("joining threads")
if useThreads:
    # first give all workers poison pills (oops, we can't give one pill, then wait for the worker, because we don't know which worker will get the next item on the queue -- might not be the one we wanted)
    for t in threadList:
        dummyWork=None
        isPoisonPill=True
        poisonPill=(isPoisonPill,dummyWork)
        jobQueue.put(poisonPill)
    # now wait for the workers to die, err, retire
    for t in threadList:
        t.join()
