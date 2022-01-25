import sys
import os
import csv
import re
import argparse
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.backends.backend_pdf
import scipy.stats
import random

from ChaseLibrary import LoadLeaderRepeatFile

matplotlib.use('pdf') # valid strings are ['GTK3Agg', 'GTK3Cairo', 'MacOSX', 'nbAgg', 'Qt4Agg', 'Qt4Cairo', 'Qt5Agg', 'Qt5Cairo', 'TkAgg', 'TkCairo', 'WebAgg', 'WX', 'WXAgg', 'WXCairo', 'agg', 'cairo', 'pdf', 'pgf', 'ps', 'svg', 'template']


parser = argparse.ArgumentParser('summarize results of multiple input files')
parser.add_argument('file', type=str, nargs='+')
parser.add_argument("--out-dir",help="directory to store results",type=str,default="summary",dest="outDir")
parser.add_argument("--actual-value",help="plot actual values, not p-values",action="store_true",dest="plotActualValue")
args = parser.parse_args()

outDir=args.outDir
inputFileNameList=args.file
plotActualValue=args.plotActualValue

os.makedirs(outDir,mode=0o777,exist_ok=True)

statsShortlistList=['avgPairProbs','helix-0-7','helix-1-11']

statNameAndFileNameToValueList=None

logicalFileNameList=[] # includes virtual files that we make from the other repeats
inputFileNameToChaseFoldingParams={}

for inputFileName in inputFileNameList:
    with open(inputFileName,'r') as inputFile:

        params=None

        reader=csv.reader(inputFile,delimiter='\t')
        statNameList=next(reader, None) # read stat names
        statNameList=statNameList[1:] # but actually, we ignore the first column

        statNameToCol={}
        for col in range(len(statNameList)):
            statNameToCol[statNameList[col]]=col

        nameAndStatsShortlist=[]
        statsShortlistColList=[statNameToCol[statName] for statName in statsShortlistList]

        statNameToValueList={}
        for statName in statNameList:
            statNameToValueList[statName]=[]

        prevOtherRepeatBaseName=None
        otherRepeatBaseCount=0
        isOtherRepeats=False
        if re.search(r'other-repeat',inputFileName) is not None:
            isOtherRepeats=True
        skipOtherRepeats_statNameToValueList={}
        for statName in statNameList:
            skipOtherRepeats_statNameToValueList[statName]=[]
        skipOtherRepeats_fileName=inputFileName+"-skip-first-repeats"

        numSamples=0
        for row in reader:
            rowName=row[0]
            statList=[float(x) for x in row[1:]]
            isActualValue=re.search(r'actual-value$',rowName) is not None

            if re.search(r'^params=',rowName) is not None:
                params=rowName
                params=params.rstrip()
                simplerParams=params
                simplerParams=re.sub(r'^params=.*ChaseFolding.py','',simplerParams)
                inputFileNameToChaseFoldingParams[inputFileName]=simplerParams

            if ((isActualValue and plotActualValue) or (not isActualValue and not plotActualValue)) and re.search(r'^params=',rowName) is None: # we don't want the values, we want the empirical p-values, and we want to skip the params, of course
                for i in range(0,len(statList)):
                    statNameToValueList[statNameList[i]].append(statList[i])

                statShortlist=[]
                for col in statsShortlistColList:
                    statShortlist.append(statList[col])
                statShortlist.append(rowName)
                nameAndStatsShortlist.append(statShortlist)

                if isOtherRepeats:
                    baseNameMatch=re.match(r"(.*)-repeat[0-9]+",rowName)
                    if not baseNameMatch:
                        raise AssertionError("based on the file name, I think we're doing other-repeats, but I couldn't match rowName %s" % (rowName))
                    baseName=baseNameMatch.group(1)
                    if baseName==prevOtherRepeatBaseName:
                        if otherRepeatBaseCount>=2:
                            for i in range(0,len(statList)):
                                skipOtherRepeats_statNameToValueList[statNameList[i]].append(statList[i])
                        otherRepeatBaseCount += 1
                    else: # new repeat
                        prevOtherRepeatBaseName=baseName
                        otherRepeatBaseCount=1 # for the next round of the loop

        # put this in our master array
        if not statNameAndFileNameToValueList:
            statNameAndFileNameToValueList={}
            for statName in statNameList:
                statNameAndFileNameToValueList[statName]={}
        for statName in statNameToValueList:
            statNameAndFileNameToValueList[statName][inputFileName]=statNameToValueList[statName]

        logicalFileNameList.append(inputFileName)

        if isOtherRepeats:
            # make the other repeats data, skipping the first couple of repeats
            for statName in statNameToValueList:
                statNameAndFileNameToValueList[statName][skipOtherRepeats_fileName]=skipOtherRepeats_statNameToValueList[statName]

            logicalFileNameList.append(skipOtherRepeats_fileName)

        # look up seqs
        inputFileMatch=re.search(r'ChaseFolding.py ([^ ]+[.]tab) ',params)
        if not inputFileMatch:
            inputFileMatch=re.search(r'ChaseFolding.py.* ([^ ]+[.]tab)$',params)
            if not inputFileMatch:
                raise AssertionError("cannot find the input file from params %s" % (params))
        leaderRepeatFileName=inputFileMatch.group(1)
        leaderRepeatList=LoadLeaderRepeatFile(leaderRepeatFileName,0,0)
        nameToLeaderRepeat={}
        for leaderRepeat in leaderRepeatList:
            (name,leaderSeq,repeatSeq)=leaderRepeat
            nameToLeaderRepeat[name]=(leaderSeq,repeatSeq)
        outFileName=outDir + '/' + os.path.basename(inputFileName)+'.tab'
        with open(outFileName, 'w', newline='') as csvfile:
            outFile = csv.writer(csvfile, delimiter='\t', quoting=csv.QUOTE_MINIMAL)
            outFile.writerow(['name','leader sequence','repeat sequence'] + statsShortlistList)
            for nameAndStats in nameAndStatsShortlist:
                name=nameAndStats.pop()
                if name not in nameToLeaderRepeat:
                    raise AssertionError('cannot find %s in leaderRepeat seqs' % (name))
                (leaderSeq,repeatSeq)=nameToLeaderRepeat[name]
                row=[name,leaderSeq,repeatSeq] + nameAndStats
                outFile.writerow(row)

numBins=20
binList=[x/numBins for x in range(0,numBins+1)] # set bins explicitly, just in case we don't observe any values of 0 or 1, so that the X-axis of each graph is comparible
if plotActualValue:
    binList=None
matplotlib.rcParams['figure.max_open_warning'] = 200

pdf = matplotlib.backends.backend_pdf.PdfPages(outDir + "/summary.pdf")
figNum=1

summaryTabFileName=outDir + "/summary.tab"
with open(summaryTabFileName,'w',newline='') as summaryTabFile:
    summaryTabCsv=csv.writer(summaryTabFile,delimiter='\t',quoting=csv.QUOTE_MINIMAL)

    for statName in statNameAndFileNameToValueList:
        #for statName in ['helix-1-11','helix-1-12']:

        useUserFileOrder=True
        if useUserFileOrder:
            fileNameList=logicalFileNameList
        else:
            fileNameList=sorted(statNameAndFileNameToValueList[statName])
            #valueListList=[value for key,value in sorted(statNameAndFileNameToValueList[statName].items())]

        valueListList=[statNameAndFileNameToValueList[statName][x] for x in fileNameList]
        weightListList=[np.ones_like(valueList)/len(valueList) for valueList in valueListList] # make the histogram sum to 1 in the obvious way, not the silly integral definition that plt.hist uses (okay, I guess it's only silly for this problem)

        # calculate the significance of this many p-values beyind < 0.05
        fileNameToPvalue={}
        fileNameToAntiPvalue={}
        if plotActualValue:
            labelList=[x+" (n="+str(len(statNameAndFileNameToValueList[statName][x]))+")" for x in fileNameList]
        else:
            for i in range(0,len(valueListList)):

                # calculate aggregated p-value for tendency to fold, and to not fold
                for forFolding in True,False:

                    method='fisher'
                    #method='stouffer'
                    indivPvalueList=valueListList[i]

                    if not forFolding:
                        # we're doing anti-folding
                        indivPvalueList=[1.0-x for x in indivPvalueList]

                    if method=='fisher':
                        indivPvalueList=[max(0.001,pvalue) for pvalue in indivPvalueList] # don't allow p-values of zero, since Fisher's method wants to take the log.  0.001 makes sense since we normally use 1000 samples.  With Laplace's Rule, we'd get 1/1002 as the estimated probability
                        stat,pvalue=scipy.stats.combine_pvalues(indivPvalueList,method=method)

                        if forFolding:
                            fileNameToPvalue[fileNameList[i]]=pvalue
                        else:
                            fileNameToAntiPvalue[fileNameList[i]]=pvalue

                pvalueFor=fileNameToPvalue[fileNameList[i]]
                pvalueAnti=fileNameToAntiPvalue[fileNameList[i]]
                inputFileName=fileNameList[i]
                params='?'
                if inputFileName in inputFileNameToChaseFoldingParams:
                    params=inputFileNameToChaseFoldingParams[inputFileName]
                summaryTabCsv.writerow([outDir,fileNameList[i],params,statName,pvalueFor,pvalueAnti])
                avgPvalue=sum(indivPvalueList)/len(indivPvalueList)
                thresh=0.01
                numLessThan05=len([x for x in indivPvalueList if x<thresh])
                #print("statName=%s,fileName=%s,average=%g,fisher=%g,num<%g=%d,frac<%g=%g" % (statName,fileNameList[i],avgPvalue,pvalue,thresh,numLessThan05,thresh,numLessThan05/len(indivPvalueList)))
                if statName=='helix-0-7' and fileNameList[i]=='out-II-A--training-noWithinRepeat.tab':
                    for i in range(0,20):
                        l=indivPvalueList
                        random.shuffle(l)
                        l=l[0:30]
                        stat,pvalue=scipy.stats.combine_pvalues(l,method=method)
                        #print("pvalue with 30=%g" % (pvalue))

            # calculate the number of samples as just the list len
            labelList=[x+" (n="+str(len(statNameAndFileNameToValueList[statName][x]))+") p<(" + str(fileNameToPvalue[x]) + ")" for x in fileNameList]

        if False: # the gray line works too, so this is unnecessary
            # and make fake data that's uniformly distributed.  easier than drawing a line
            valueListList.append([float(x)/(numBins)+1/(2*numBins) for x in range(0,numBins)])
            weightListList.append([1/(numBins) for x in range(0,numBins)])
            fileNameList.append('random')

        fig=plt.figure(figsize=(4,8))
        plt.title('histogram for '+statName)
        #weights = np.ones_like(valueListList) / len(valueListList)
        n,bins,patches = plt.hist(valueListList,bins=binList,weights=weightListList,label=labelList,histtype='bar')
        #n,bins,patches = plt.hist([[1,1,1,1,1,1,2,3],[1,2,3,3],[1,1,2,2,3,3,3,3]],normed=True,label=['A','B','C'])
        if not plotActualValue:
            plt.plot([0,1],[0.05,0.05],color='gray',linestyle='--')
        plt.legend()
        #plt.xticks(rotation='vertical')
        if plotActualValue:
            plt.xlabel('value of statistic')
            plt.ylabel('fraction of repeats with value in this bin')
        else:
            plt.xlabel('p-value (in bins of width 0.05)')
            plt.ylabel('fraction of repeats with p-value in this bin')
        #plt.show()
        plt.savefig(outDir+'/'+statName+'.pdf',format='pdf')
        pdf.savefig(figNum)

        figNum += 1

pdf.close()
