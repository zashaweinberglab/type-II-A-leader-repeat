import sys
import subprocess
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

# function ParseSsToPairList: 'ss' is a string with brackets (angle brackets, square brackets, round brackets or French braces) indicating a pseudoknot-free secondary structure.  This function returns a list of base pairs.  Each element of the returned list is a 2-tuple i,j that correspond to the string index of the left & right nucleotides in a base pair
def ParseSsToPairList (ss):
    pairList=[]
    stack=[]
    for j in range(len(ss)):
        if ss[j]=='(' or ss[j]=='[' or ss[j]=='<' or ss[j]=='{':
            stack.append(j)
        if ss[j]==')' or ss[j]==']' or ss[j]=='>' or ss[j]=='}':
            if not stack:
                raise AssertionError('parentheses are unbalanced {} (no left pair at pos={})'.format(ss,j))
            i=stack.pop();
            pairList.append([i,j])
    if stack:
        raise AssertionError('parentheses are unbalanced {} (at the 3\' end of the molecule, I realize there\'s too many left pairs)'.format(ss))

    return pairList

# function ParseSsToArrayMap: same as ParseSsToPairList, but return the list as an array, where each item is the string index of its pair partner, or -1 if the nucleotide is unpaired
def ParseSsToArrayMap (ss):
    pairList=ParseSsToPairList(ss)
    arrayMap=np.full((len(ss)), -1)
    for pair in pairList:
        (i,j)=pair
        arrayMap[i]=j
        arrayMap[j]=i
    return arrayMap

# function HasDegenNucOrJunk: we require sequences to consist of the normal 4 nucleotides (but T and U are both accepted).  This function says whether a sequence has another letter in it (which mostly represent degenerate nucleotides)
def HasDegenNucOrJunk (seq):
    return re.search(r'[^ACGTUacgtu]',seq) # has something that's not a legal & normal nuc

# function LoadLeaderRepeatFile: load a file in the tab-delimited format that one gets by copying the Excel spreadsheets that Omer sends into a text file.  leaderRepeatsFileName: the file name. moveNucsFromLeaderToRepeat: after loading leader & repeat seqs, move some of the leader nucs to the repeat, so it's as if they were part of the repeat. maxLeaderLen: truncate leaders longer than this amount.  In the loaded file, a line with column headers is optional, but if it exists, the first column must be labeled  "name" or "ID" (case sensitive).  The fields are: name, leader sequence, repeat sequence.  Additional fields are allowed, but ignored.
def LoadLeaderRepeatFile(leaderRepeatsFileName,moveNucsFromLeaderToRepeat,maxLeaderLen):
    leaderRepeatList=[]
    with open(leaderRepeatsFileName) as leaderRepeatsFile:
        reader=csv.reader(leaderRepeatsFile, delimiter='\t', quotechar='\"')
        for row in reader:
            (name,leaderSeq,repeatSeq)=row[0:3]  # allow additional fields

            leaderSeq=leaderSeq.upper()
            repeatSeq=repeatSeq.upper()

            # detect optional column headers
            if name=="name" or name=="ID":
                #print("bloop")
                continue

            leaderSeq=re.sub(r'^N+','', leaderSeq)
            if HasDegenNucOrJunk(leaderSeq) or HasDegenNucOrJunk(repeatSeq):
                print("warning: sequence for %s has N nucleotides in the middle -- skipping" % (name))
                continue                    

            if maxLeaderLen>=0: # if it's set
                if len(leaderSeq)>maxLeaderLen:
                    print("%s before: %s  %s" % (name,leaderSeq,repeatSeq))
                    if maxLeaderLen==0: # the [-maxLeaderLen:] syntax seems to not work when maxLeaderLen is zero, which is annoying
                        leaderSeq=""
                    else:
                        leaderSeq=leaderSeq[-maxLeaderLen:]
                    print("%s  after: %s  %s" % (name,leaderSeq,repeatSeq))

            if moveNucsFromLeaderToRepeat>0:
                #print("%s before: %s  %s" % (name,leaderSeq,repeatSeq))
                remainingLeaderSeqLen=len(leaderSeq)-moveNucsFromLeaderToRepeat
                repeatSeq=leaderSeq[remainingLeaderSeqLen:] + repeatSeq
                leaderSeq=leaderSeq[0:remainingLeaderSeqLen]
                #print("%s  after: %s  %s" % (name,leaderSeq,repeatSeq))

            leaderRepeatList.append([name,leaderSeq,repeatSeq])

    return leaderRepeatList

# function GetPairProbFromProbArray: given an array of base-pair probabilities returned by the ViennaRNA library, calculate the probabilities that each nucleotide is involved in some base pair (with any other nucleotide)
def GetPairProbFromProbArray(basepair_probs,seq):
    pairProbArray=np.zeros(len(seq))
    
    for i in range(1, len(basepair_probs)):
        for j in range(i+1,len(basepair_probs[i])):
            pairProbArray[i-1] += basepair_probs[i][j] # remember that Vienna's arrays are 1-based, but I'm doing arrays zero-based
            pairProbArray[j-1] += basepair_probs[i][j]

    return pairProbArray

# function GetPairProbFromProbArray_OptionalAllowWithinRepeat: same as previous function, but optionally ignore base pairs within the repeat, so that we are looking at the probabilities that repeat nucleotides pair with the leader sequence.  This function is used to generate the blue-dot figure, where a scale from 0 (white) to 1 (blue) is the probability of pairing with the leader for each nucleotide in the repeat
def GetPairProbFromProbArray_OptionalAllowWithinRepeat(basepair_probs,seq,firstRepeatPos,doAllowWithinRepeat):
    pairProbArray=np.zeros(len(seq))
    
    for i in range(1, len(basepair_probs)):
        for j in range(i+1,len(basepair_probs[i])):
            if doAllowWithinRepeat or (i<firstRepeatPos and j>=firstRepeatPos):
                pairProbArray[i-1] += basepair_probs[i][j] # remember that Vienna's arrays are 1-based, but I'm doing arrays zero-based
                pairProbArray[j-1] += basepair_probs[i][j]

    return pairProbArray

# function ChaseFindLeaderAndRepeatMfe: paste the leader and repeat sequence together, and find the MFE structure using the ViennaRNA library.  Base pairs that are partially or fully within the repeat sequence are changed to use square brackets (instead of the round brackets returned by ViennaRNA)
def ChaseFindLeaderAndRepeatMfe(leaderSeq,repeatSeq,md):
    seq=leaderSeq+repeatSeq
    fc = RNA.fold_compound(seq, md)
    (mfe_ss, mfe) = fc.mfe()
    ss=mfe_ss
    ssAsList=list(ss)
    pairList=ParseSsToPairList(ss)
    for i,j in pairList:
        if j>=len(leaderSeq):
            ssAsList[i]='['
            ssAsList[j]=']'
            ss="".join(ssAsList)
            leaderSs=ss[0:len(leaderSeq)]
            repeatSs=ss[len(leaderSeq):]
    return leaderSs,repeatSs

# function ChaseCalcBppAsDigitStr: for each nucleotide in the combined leader-repeat seq, calculate the probability that the nucleotide pairs with some other nucleotide (using GetPairProbFromProbArray), and convert these probabilities into a string of numbers where the numbers 0,1,...,9 mean >= 0%, 10%, ..., >=90%.  An asterisk means 100%.  This function is used for debugging, since it's easy to look at in text.
def ChaseCalcBppAsDigitStr(leaderSeq,repeatSeq,md):
    seq=leaderSeq+repeatSeq
    fc = RNA.fold_compound(seq,md)
    (propensity, ensemble_energy) = fc.pf()
    basepair_probs = fc.bpp()
    pairProbArray=GetPairProbFromProbArray(basepair_probs,seq)
    digitList=[]
    for i in range(0,len(seq)):
        p=pairProbArray[i]
        ip=int(p*10.0)
        if ip>=0 and ip<10:
            digit=chr(ord('0')+ip)
        elif ip==10:
            digit='*'
        else:
            raise AssertionError("got %d from %g" % (ip,p))
        digitList.append(digit)
    digitStr="".join(digitList)
    return digitStr

# function ChaseDrawLeaderRepeat: make a PDF image with the combined leader/repeat sequence with colors indicating the probability of base pairing.  digitStr is the result of the ChaseCalcBppAsDigitStr function.  The drawing is performed using R2R.  This function was just for visualization.
def ChaseDrawLeaderRepeat(name,leaderSeq,repeatSeq,leaderSs,repeatSs,digitStr,drawPdfFileName,drawStoFileName,drawMetaFileName):
    with open(drawStoFileName,'w') as drawStoFile,open(drawMetaFileName,'w') as drawMetaFile:
        repeatIndicator_leader="." * len(leaderSeq)
        repeatIndicator_repeat="x" * len(repeatSeq)
        repeatIndicator=repeatIndicator_leader + repeatIndicator_repeat
        print("# STOCKHOLM 1.0\n%s %s%s\n#=GC R2R_XLABEL_repeat %s\n#=GC SS_cons %s%s\n#=GC R2R_XLABEL_bpp %s\n" % (name,leaderSeq,repeatSeq,repeatIndicator,leaderSs,repeatSs,digitStr),file=drawStoFile)
        print("#=GF R2R shade_along_backbone repeat:x rgb:192,192,192\n",file=drawStoFile)
        print("#=GF R2R SetDrawingParam nucShrinkWithCircleNuc 1 indicateOneseqWobblesAndNonCanonicals false defaultOneseqLabeling false\n",file=drawStoFile)

        minrgb=[255,255,255]
        maxrgb=[0,0,255]

        minrgb=[255,216,216] # try paler colors
        maxrgb=[200,216,255]
        
        for digit in range(0,10):
            if str(digit) in digitStr: # otherwise R2R will report label not found
                x=digit/10.0
                rgb=[]
                for i in range(0,3):
                    low=minrgb[i]
                    high=maxrgb[i]
                    rgb.append(int((high-low)*x+low))
                print("#=GF R2R circle_nuc bpp:%d rgb:%d,%d,%d width -1" % (digit,rgb[0],rgb[1],rgb[2]),file=drawStoFile)
        if '*' in digitStr:
            print("#=GF R2R circle_nuc bpp:* rgb:%d,%d,%d width -1" % (maxrgb[0],maxrgb[1],maxrgb[2]),file=drawStoFile)

        print("//\n",file=drawStoFile)

        print("%s\toneseq\t%s\n" % (drawStoFileName,name),file=drawMetaFile)

    result=subprocess.call("r2r --disable-usage-warning %s %s" % (drawMetaFileName,drawPdfFileName),shell=True)
    if result!=0:
        raise AssertionError("r2r command returned result %d",result)
