
        prevNumBad=0
        prevMostNumPairs=0
        for numBad,mostNumPairs in sorted(numBadToMaxLen.items()):

            isLastOne=False

            if numBad>=maxNumBad:
                numBad=maxNumBad
                isLastOne=True

            if mostNumPairs>=maxNumPairs:
                mostNumPairs=maxNumPairs:
                isLastOne=True

            for thisNumBad in range(prevNumBad,numBad):
                for thisNumPairs in range(prevMostNumPairs,mostNumPairs):
                    countsByNumBadAndNumPairs[thisNumBad,thisNumPairs] += 1
            
            prevNumBad=numBad
            prevNumPairs=mostNumPairs

            if isLastOne:
                break



    for numBad in range(0,maxNumBad):
        for numPairs in range (0,maxNumPairs):
            statMap["helix-{}-{}".format(numBad,numPairs)] = countsByNumBadAndNumPairs[numBad,numPairs]/numBoltzmannSamples





                                currNumPairs += 1 # we get one more pair...
                    prevL=currL
                    prevR=currR
                    
                currR += 1

            if numBad not in numBadToMaxLen or numBadToMaxLen[numBad]<currNumPairs:
                #print("currNumPairs=%d" % (currNumPairs))
                numBadToMaxLen[numBad]=currNumPairs

    # propagate values.  We interpret it as **UP TO** numBad skips.  Anything we can do with up to B numBad we can also do with B+1
    bestNumPairsWithPrevNumBad=0;
    for numBad in sorted(numBadToMaxLen):
        if bestNumPairsWithPrevNumBad > numBadToMaxLen[numBad]:
            numBadToMaxLen[numBad] = bestNumPairsWithPrevNumBad
        bestNumPairsWithPrevNumBad = numBadToMaxLen[numBad]
    return numBadToMaxLen

