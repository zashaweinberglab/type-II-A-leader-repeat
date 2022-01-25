import sys
#print ('\n'.join(sys.path))
import math
import numpy as np
import RNA
import csv
import string
import time
import operator
import argparse
import random
import re
from operator import itemgetter
from matplotlib import pyplot as plt
from matplotlib import collections  as mc
import cairo
from cairo import PDFSurface
#holy crap, this Biopython stuff doesn't want to enable a recursive traversal: from Bio import Phylo
from skbio import TreeNode

from ChaseLibrary import ParseSsToArrayMap,LoadLeaderRepeatFile,GetPairProbFromProbArray,GetPairProbFromProbArray_OptionalAllowWithinRepeat

# [zasha@wl chase]$ python3 ~/bliss/code/motifs_2007/rpl/ChaseFigureFromFolding.py figure-I-E/I-E-split--training.tab -o t.tab --use-prob --no-within-repeat --phylip-tree I-E-clustered-myNames.phylip-tree

if sys.version_info < (3,0,0):
    print(__file__ + ' requires Python 3, while Python ' + str(sys.version[0] + ' was detected. Terminating. '))
    sys.exit(1)

def DumpTree(node,leafNameList):
    
    if node.is_tip():
        leafNameList.append(node.name)
    else:
        for child in node.children:
            DumpTree(child,leafNameList)

def TreeLongestLen(node):
    thisLength=node.length
    if thisLength is None:
        thisLength=0
    if node.is_tip():
        return node.length
    else:
        return thisLength + max(TreeLongestLen(child) for child in node.children)

def NucAndSsToHtmlChar (nuc,dotBracket):
    if dotBracket=='.':
        return nuc
    if dotBracket=='(' or dotBracket==')':
        #return "<u>"+nuc+"</u>"
        return "<b>"+nuc+"</b>"
    raise AssertionError("unexpected dotBracket character %s" % (dotBracket))

def DrawTree(node,context,nameToY,leafX,lengthFactor,x0,y0): # return our x,y coordinate, so that the parent can lay itself out relative to us
    if node.is_tip():
        # figure out coordinate based on which leaf we are
        y= - nameToY[node.name]
        x=-leafX
        return x,y
    else:
        # find out coords of children
        childCoordList=[]
        for child in node.children:
            childX,childY=DrawTree(child,context,nameToY,leafX,lengthFactor,x0,y0)
            length=child.length*lengthFactor
            childX -= length
            childCoordList.append([childX,childY,length])
        if len(childCoordList)!=2:
            raise AssertionError()
        x1,y1,l1=childCoordList[0]
        x2,y2,l2=childCoordList[1]
        #print("%f,%f,%f,%f" % (x1,y1,x2,y2))
        if abs(x1-x2)>lengthFactor*1e-3:
            raise AssertionError()

        # draw lines
        context.move_to(x1+l1-x0,y1-y0)
        context.line_to(x1-x0,y1-y0)
        context.line_to(x2-x0,y2-y0)
        context.line_to(x2+l2-x0,y2-y0)
        context.stroke()

        # calc our coordinates
        x=x1
        y=(y1+y2)/2.0
        return x,y

def DrawFigure (context,x0,y0,plotX,plotY,plotC,circDiameter,totalTreeWidth,leafX,maxX,maxY,stepX,nameAndYList,moveNucsFromLeaderToRepeat,trainingSet,experimentedSet,tree):

    # gray lines along the sequences
    if False: # can't see the lines because they're blocked by the circles. Could do something with alpha, but I think the lines would just be distracting.
        context.set_source_rgb(0.6,0.6,0.6)
        context.set_line_width(0.5)
        for name,y in nameAndYList:
            context.move_to(0-x0,-y-y0)
            context.line_to(maxX-x0,-y-y0)
            context.stroke()

    # separator between leader-nucs-moved-to-repeat and actual repeat nucs
    if moveNucsFromLeaderToRepeat>0:
        context.set_source_rgb(0,0,0)
        context.set_line_width(0.5)
        #print("%f,%d,%f" % (stepX,moveNucsFromLeaderToRepeat,maxY))
        context.move_to(stepX*moveNucsFromLeaderToRepeat-x0,0-y0)
        context.line_to(stepX*moveNucsFromLeaderToRepeat-x0,-maxY-y0)
        context.stroke()

    # marking training set thingies
    if False:
        context.set_source_rgb(0,0,0)
        context.set_line_width(0.5)
        for name,y in nameAndYList:
            if name not in trainingSet:
                context.move_to(maxX+circDiameter/2-x0,-y-circDiameter/2-y0)
                context.line_to(maxX+circDiameter/2-x0,-y+circDiameter/2-y0)
                context.stroke()

    context.set_line_width(0)
    for point in zip(plotX,plotY,plotC):
        (x,y,c)=point
        #print("%g,%g,%g" % (x,-y,circDiameter))
        context.arc(x-x0,-y-y0,circDiameter/2.0,0,2.0*math.pi)
        context.set_source_rgb(1-c, 1-c  ,1)
        context.fill()

    context.set_font_size(circDiameter)
    context.select_font_face("Arial",cairo.FONT_SLANT_NORMAL,cairo.FONT_WEIGHT_NORMAL)
    context.set_source_rgb(0,0,0)

    if tree:
        treeLongestLen=TreeLongestLen(tree)
        targetTreeLen=totalTreeWidth
        #print("treeLongestLen=%f , target=%f" % (treeLongestLen,targetTreeLen))
        lengthFactor=targetTreeLen/treeLongestLen

        nameToY={name:y for name,y in nameAndYList}
        context.set_source_rgb(0,0,0)
        context.set_line_width(0.5)
        rootX,rootY=DrawTree(tree,context,nameToY,leafX,lengthFactor,x0,y0)
        # draw the root's branch
        context.move_to(rootX-x0,rootY-y0)
        context.line_to(rootX-circDiameter/2.0-x0,rootY-y0)
        context.stroke()

    if False: # label each line with its name (this was for syncing the tree with the dots, which is no longer necessary now that I draw the tree too)
        maxY=None
        for name,y in nameAndYList:
            #print("%s,%d" % (name,y))
            label=name
            if name in experimentedSet:
                label += " ****"
            context.move_to(maxX+circDiameter-x0,-(y-circDiameter/2.0)-y0)
            #print(label)
            context.show_text(label)
            if maxY is None or y>maxY:
                maxY=y

    if False: # draw the gradient from blue to white
        rectLeft=0
        rectRight=0+circDiameter*10
        rectTop=maxY+3*circDiameter
        rectBot=maxY+3*circDiameter+2*circDiameter
        lg = cairo.LinearGradient(rectLeft,rectTop,rectRight,rectTop)
        lg.add_color_stop_rgb(0,1,1,1)
        lg.add_color_stop_rgb(1,0,0,1)

        context.rectangle(rectLeft-x0,-rectTop-y0,rectRight-rectLeft,rectBot-rectTop)
        context.set_source(lg)
        context.fill()


params=" ".join(sys.argv)
parser = argparse.ArgumentParser('folding stats for type II-A CRISPR leader/repeat seqs')
parser.add_argument("inFileName", help="input file name (tab-delimited, fields are: name,leaderSeq,repeatSeq)",type=str)
parser.add_argument("-o",dest="outFile",help="output .tab file name",type=str,default="out.tab")
parser.add_argument("--use-mea",help="Use MEA structures -- DOES NOT WORK",dest="useMea",action="store_true")
parser.add_argument("--use-prob",help="plot base pair probs",dest="useProb",action="store_true")
parser.add_argument("--no-within-repeat",dest="noWithinRepeat",help="ignore pairs contained within repeat seq",action="store_true")
parser.add_argument("--extend-repeat-by",help="move this many nucs on the 3' end of the leader into the 5' end of repeat, to look for hairpins that are just beyond the repeat",type=int,default=0,dest="moveNucsFromLeaderToRepeat")
parser.add_argument("--tree-order",help="tree file output by clustalo that gives the desired order of the inputs",type=str,dest="treeOrderFileName")
parser.add_argument("--phylip-tree",help="use this tree and plot the tree itself",type=str,dest="phylipTreeFileName")
parser.add_argument("--max-leader-len",help="if a leader is longer than this amount, chop nucs from the 5' end to truncate it to this length",type=int,default=-1,dest="maxLeaderLen")
parser.add_argument("--mark-training",help="mark examples that were part of training data, i.e. in the given file",type=str,dest="markTraining")
parser.add_argument("--mark-experimented",help="mark the examples in the given file as being experimented upon",type=str,dest="markExperimented")
args = parser.parse_args()

leaderRepeatsFileName=args.inFileName
outFileName=args.outFile
htmlFileName=outFileName+".html"
pdfFileName=outFileName+".pdf"
pdfDrawFileName=outFileName+".draw.pdf"
useMea=args.useMea
useProb=args.useProb
disableMatplotlib=True # doesn't give enough control
doAllowWithinRepeat=not args.noWithinRepeat
moveNucsFromLeaderToRepeat=args.moveNucsFromLeaderToRepeat
treeOrderFileName=args.treeOrderFileName
phylipTreeFileName=args.phylipTreeFileName
maxLeaderLen=args.maxLeaderLen
markTrainingFileName=args.markTraining
markExperimentedFileName=args.markExperimented

leaderRepeatList=LoadLeaderRepeatFile(leaderRepeatsFileName,moveNucsFromLeaderToRepeat,maxLeaderLen)

trainingSet=set()
if markTrainingFileName:
    markTrainingLeaderRepeatList=LoadLeaderRepeatFile(markTrainingFileName,0,0)
    trainingSet=set([name for name,leaderSeq,repeatSeq in markTrainingLeaderRepeatList])

experimentedSet=set()
if markExperimentedFileName:
    markExperimentedLeaderRepeatList=LoadLeaderRepeatFile(markExperimentedFileName,0,0)
    experimentedSet=set([name for name,leaderSeq,repeatSeq in markExperimentedLeaderRepeatList])
    print(experimentedSet)

tree=None
orderNameList=[]
if phylipTreeFileName:
    print("loading Newick tree")
    tree = TreeNode.read(phylipTreeFileName,format='newick',convert_underscores=False)

    leafNameList=[]
    DumpTree(tree,leafNameList)
    #print(leafNameList)

    orderNameList=leafNameList

if treeOrderFileName:
    if orderNameList:
        raise ValueError("multiple ways of determining leaf order")

    with open(treeOrderFileName) as treeOrderFile:
        for line in treeOrderFile:
            # highly specific to output of clustalo.  I assume that lines beginning with a letter or number are sequences, and punctuation symbols refers to the structure of the tree or branch lengths
            if re.search(r'^[A-Za-z0-9]',line):
                line=line.rstrip()
                line=re.sub(r':.*$','',line) # remove branch length, if any
                orderNameList.append(line)

if orderNameList:
    print("setting name order")
    leaderRepeatMap={name:(name,leaderSeq,repeatSeq) for name,leaderSeq,repeatSeq in leaderRepeatList}
    #print(sorted(leaderRepeatMap))
    #print(orderNameList)
    old_leaderRepeatList=leaderRepeatList
    leaderRepeatList=[]
    for name in orderNameList:
        if name not in leaderRepeatMap:
            # this can happen because of N nucleotides, so we don't have any internal check
            print("warning: name %s not found in list" % (name))
            continue
            #raise AssertionError("something's wrong: name %s not found" % (name))
        leaderRepeatList.append(leaderRepeatMap[name])
    

stepX=8
stepY=16
plotX=[]
plotY=[]
plotC=[]
listY=[]
minX=0
maxX=0
maxY=0
circDiameter=3
totalTreeWidth=11.0/25.4*72.0 # 11 mm, convert to inches then to points
leafX=0.15/25.4*72.0
startX=circDiameter/2.0
startY=circDiameter/2.0
currY=startY

stepY=circDiameter # more compact
stepX=circDiameter

print("writing")
with open(outFileName, 'w', newline='') as csvfile,open(htmlFileName,'w') as htmlFile:
    
    outFile = csv.writer(csvfile, delimiter='\t', quoting=csv.QUOTE_MINIMAL)
    htmlFile.write("<html><body><table>\n")

    outFile.writerow(['shortest MFE pair dist to repeat', 'longest MFE pair dist to repeat','pairs in repeat'])
    htmlFile.write("<tr><th>Longest MFE pair dist to repeat</th><th>Shortest MFE pair dist to repeat</th><th>Sequence</th></tr>")

    nameAndYList=[]
    
    for leaderRepeat in leaderRepeatList: # we used to try to sort it ourselves: sorted(leaderRepeatList, key=itemgetter(2)):
        (name,leaderSeq,repeatSeq)=leaderRepeat

        currY -= stepY # this goes down the page, so the order's consistent
        listY.append(currY)
        maxY=min(maxY,currY)

        nameAndYList.append((name,currY))

        seq=leaderSeq+repeatSeq
        seq=seq.translate(str.maketrans('acgtT','ACGUU'))
                
        md = RNA.md()
        fc = RNA.fold_compound(seq, md)

        if useProb:
            firstRepeatSeqPos=len(leaderSeq)
            lastRepeatSeqPos=len(leaderSeq)+len(repeatSeq)

            (propensity, ensemble_energy) = fc.pf()
            basepair_probs = fc.bpp()
            firstRepeatSeqPos_forProbArray=firstRepeatSeqPos
            # was: (I think this was a mistake, although maybe I want to consider pairs within the leader-as-repeat region and the real repeat.  But by setting them equally, I'm doing the same thing as when I run ChaseFolding.py
            #firstRepeatSeqPos+moveNucsFromLeaderToRepeat # allow pairs from the part of the leader that got moved into the repeat, and the actual repeat
            pairProbArray=GetPairProbFromProbArray_OptionalAllowWithinRepeat(basepair_probs,seq,firstRepeatSeqPos_forProbArray,doAllowWithinRepeat)
            for right in range(firstRepeatSeqPos,lastRepeatSeqPos):

                plotX.append(startX+stepX*(right-firstRepeatSeqPos))
                plotY.append(currY)
                plotC.append(pairProbArray[right])

            minX=0
            maxX=max(maxX,stepX*(lastRepeatSeqPos-firstRepeatSeqPos))

        else:

            ss=""
            if useMea:
                raise AssertionError("MEA doesn't work.  The library function has a different from from the example code, and it doesn't seem worth bothering worth")
                (MEA_struct, MEA) = fc.MEA()
                ss=MEA_struct
            else:
                (mfe_ss, mfe) = fc.mfe()
                ss=mfe_ss

                repeat_ss=ss[len(leaderSeq):]
                repeatSeqHtml = "".join(NucAndSsToHtmlChar(nuc,dotBracket) for nuc,dotBracket in zip(repeatSeq,repeat_ss))

            firstLeaderPos=0
            lastLeaderPos=len(leaderSeq)
            firstRepeatPos=lastLeaderPos
            lastRepeatPos=firstRepeatPos+len(repeatSeq)
                
            # we could probably do this with list comprehensions and min/max, but I think it's a bit more complicated since we have to exclude pairing within the repeatSeq
            # also, we could use the no-pseudoknot property and search strings for brackets, but I'm concerned this will get unnecessarily complicated
            minCoordOfPairInLeader=lastLeaderPos
            maxCoordOfPairInLeader=0
            pairArrayMap=ParseSsToArrayMap(ss)
            for right in range(firstRepeatPos,lastRepeatPos):
                left=pairArrayMap[right]
                if left!=-1:
                    plotX.append(startX+stepX*(right-firstRepeatPos))
                    plotY.append(currY)
                    maxX=max(maxX,stepX*(lastRepeatPos-firstRepeatPos))
                    if left>=firstLeaderPos and left<lastLeaderPos:
                        minCoordOfPairInLeader=min(minCoordOfPairInLeader,left)
                        maxCoordOfPairInLeader=max(maxCoordOfPairInLeader,left)
                        
                        shortestDistFromRepeat = firstRepeatPos - maxCoordOfPairInLeader
                        longestDistFromRepeat = firstRepeatPos - minCoordOfPairInLeader

            debugHtml="" # "<td>%s : %d,%d</td>" % (ss,minCoordOfPairInLeader,maxCoordOfPairInLeader)
            htmlFile.write("<tr><td>%d</td><td>%d</td><td>%s</td>%s</tr>\n" % (shortestDistFromRepeat,longestDistFromRepeat,repeatSeqHtml,debugHtml))


    outFile.writerow(["params=%s\n" % (params)])
    htmlFile.write("</table>\n")
    htmlFile.write("<p>params=%s</p>\n" % (params))
    htmlFile.write("</body></html>\n")

    if not disableMatplotlib:
        print("doing matplotlib")
        fig=plt.figure(figsize=(4,8))
        plt.title('Nucleotides that pair in MFE')
        if useProb:
            print("\twith useProb")
            plotCFlip=[1-c for c in plotC]
            plt.scatter(plotX,plotY,c=plotCFlip,cmap='gray',s=20)
        else:
            plt.plot(plotX,plotY,'bo')
            plt.xticks([])
            plt.yticks([])
            lineList=[[(minX,y),(maxX,y)] for y in listY]
            lineCollection=mc.LineCollection(lineList)
            fig.axes[0].add_collection(lineCollection)
            #plt.show()
        plt.savefig(pdfFileName)

    rectWidth=5
    rectHeight=5
    #print("%g,%g" % (maxX,-currY))


    # I think it's nicer if the PDF gets automatically sized, s.t. I don't have to adjust things.  But if I just use RecordingSurface the normal way, I can't figure out how to set the background color.  Either I do this as the first thing, when the bounding box is unknown, and I get weird values from the RecordingSurface for the width and height.  Or I set the background afterwards, in which case it seems to get ignored.  So, I'll just do it in two phases
    virtualSurface = cairo.RecordingSurface(cairo.CONTENT_COLOR, None)
    context = cairo.Context(virtualSurface)

    print("drawing figure")

    DrawFigure (context,0,0,plotX,plotY,plotC,circDiameter,totalTreeWidth,leafX,maxX,maxY,stepX,nameAndYList,moveNucsFromLeaderToRepeat,trainingSet,experimentedSet,tree)

    x0,y0,width,height=virtualSurface.ink_extents()
    #print("%g,%g,%g,%g" % (x0,y0,width,height))

    surface=PDFSurface(pdfDrawFileName,width,height)
    pdfContext = cairo.Context(surface)

    pdfContext.rectangle(0, 0, width, height)
    pdfContext.set_source_rgb(1, 1, 1)
    pdfContext.fill()

    DrawFigure (pdfContext,x0,y0,plotX,plotY,plotC,circDiameter,totalTreeWidth,leafX,maxX,maxY,stepX,nameAndYList,moveNucsFromLeaderToRepeat,trainingSet,experimentedSet,tree)
        
    surface.show_page()
