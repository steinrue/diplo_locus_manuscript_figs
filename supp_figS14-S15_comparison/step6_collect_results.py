import numpy
import pandas
import pathlib
import matplotlib
import matplotlib.pyplot as plt
import seaborn




# where to find the results of the different methods
RESULTS_DIR = {
    'diplo_locus' : pathlib.Path ("results_diplo_locus"),
    'lls' : pathlib.Path ("results_lls"),
    'wfabc' : pathlib.Path ("results_wfabc"),
    'bmws' : pathlib.Path ("results_bmws"),
}
SLURMIFY_DIR = 'slurmify_analyze'
SLURMIFY_SKELETON = 'slurm_analyze_%s.sbatch'
# max time jobs were allowed to run in seconds
# this is 24 hours
MAX_TIME = 86400
# this is 8 hours
# MAX_TIME = 28800
# 2000 for each of 6 selection coefficients
NUM_REPLICATES = 12000




def load_diplo_locus_mles (resultsDir, scenarios):

    assert (resultsDir.is_dir())

    mlesDiploLocus = {}
    for (scIdx, thisScenario) in enumerate(scenarios):
        print (thisScenario)

        thisResultsFile = pathlib.Path (resultsDir, f'{thisScenario}_results_off-grid_maxLLs.txt')
        # if file doesn't exist, we likely didn't compute it
        if (not thisResultsFile.is_file()):
            # but we should have computed all
            assert (False)
            print (f"Results file not found, assuming that analysis not performed: {str(thisResultsFile)}")
            continue

        thisFrame = pandas.read_csv (thisResultsFile, sep='\t', comment='#')
        # thisFrame.columns = ['id','ongrid_s2hat', 'ongrid_maxLogLikelihood','s2Hat', 'maxLogLikelihood', 'MLR']
        thisFrame.rename (columns={'ID' : 'id', 's2hat' : 's2Hat'}, inplace=True)

        # save it
        mlesDiploLocus[thisScenario] = thisFrame
        # mlesDiploLocus[thisScenario] = thisFrame['s2Hat'].to_numpy()

    return mlesDiploLocus


def load_lls_mles (resultsDir, scenarios):

    assert (resultsDir.is_dir())

    mlesLLS = {}
    for (scIdx, thisScenario) in enumerate(scenarios):
        print (thisScenario)

        thisResultsFile = pathlib.Path (resultsDir, f'{thisScenario}_lls_results.table')
        # if file doesn't exist, we likely didn't compute it
        if (not thisResultsFile.is_file()):
            print (f"Results file not found, assuming that analysis not performed: {str(thisResultsFile)}")
            continue

        thisFrame = pandas.read_csv (thisResultsFile, sep='\t', comment='#')
        thisFrame.rename (columns={'sEstimate' : 's2Hat', 'ID' : 'id'}, inplace=True)

        # save it
        mlesLLS[thisScenario] = thisFrame
        # mlesLLS[thisScenario] = thisFrame['sEstimate'].to_numpy()

    return mlesLLS


def load_wfabc_mles (resultsDir, scenarios):

    assert (resultsDir.is_dir())

    mlesWFABC = {}
    for (scIdx, thisScenario) in enumerate(scenarios):
        print (thisScenario)

        # we also need the ids, because they indicate the true value
        thisIdFile = pathlib.Path (resultsDir, f'{thisScenario}.ids')
        assert (thisIdFile.is_file()), str(thisIdFile)
        # should just be a list
        ifs = open (thisIdFile, 'r')
        # get rid of newlines
        theseIds = [x.strip() for x in ifs.readlines()]
        ifs.close()

        # here is the file (only need the s file)
        thisResultsFile = pathlib.Path (resultsDir, f'{thisScenario}_loci_posterior_s.txt')
        # if file doesn't exist, we likely didn't compute it
        if (not thisResultsFile.is_file()):
            print (f"Results file not found, assuming that analysis not performed: {str(thisResultsFile)}")
            continue

        # each row is simulated replicate
        # each column is sample from posterior of s
        # separated by whitespace
        posteriorSamples = numpy.loadtxt (thisResultsFile)
        # since some didn't finish, we can only expect less
        assert (len(theseIds) >= posteriorSamples.shape[0]), (len(theseIds), posteriorSamples.shape[0])

        # get posterior mean as estiamte of s
        sEstimates = numpy.mean (posteriorSamples, axis=1)

        # put them in a frame
        thisFrame = pandas.DataFrame ({
            "id" : theseIds[:len(sEstimates)],
            "s2Hat" : sEstimates
        })

        # save it
        mlesWFABC[thisScenario] = thisFrame
        # mlesWFABC[thisScenario] = thisFrame['sHat'].to_numpy()

    return mlesWFABC


def load_bmws_mles (resultsDir, scenarios):

    assert (resultsDir.is_dir())

    mlesBMWS = {}
    for (scIdx, thisScenario) in enumerate(scenarios):
        print (thisScenario)

        # if files for certain scenario don't exist, just skip it, we likely didn't run it
        thisResultsFile = pathlib.Path (resultsDir, f'{thisScenario}.results')
        if (not thisResultsFile.is_file()):
            print (f"Results file not found, assuming that analysis not performed: {str(thisResultsFile)}")
            continue

        # we have a results file, so proceed to collect results
        colNames = ['chr', 'pos', 'id', 'ref', 'alt', 'freq', 's2Hat', 'sl1', 'sl2']
        resultsFrame = pandas.read_csv (thisResultsFile, sep='\t', header=None, names = colNames)
        # sign of selection coefficient is reversed for some reason
        resultsFrame['s2Hat'] *= -1

        # save it
        mlesBMWS[thisScenario] = resultsFrame
        # mlesBMWS[thisScenario] = resultsFrame['sHat'].to_numpy()

    return mlesBMWS


def theBoxPlots (scenarios, allMleS, nanFraction):

    # seaborn is *^&% awesome
    seaborn.set_theme (style="whitegrid", font_scale=0.8)

    # only plot some s
    sToPlot = {
        'g4000' : [f'{x:.4f}' for x in [0, 0.002, 0.005]],
        'g160' : [f'{x:.4f}' for x in [0, 0.02, 0.05]],
    }

    # some list and names for the methods
    methodList = ['diplo_locus', 'lls', 'wfabc', 'bmws']
    methodNames = {
        'diplo_locus' : 'diplo-locus',
        'lls' : 'LLS',
        'wfabc' : 'WFABC',
        'bmws' : 'bmws',
    }
    methodColors = {
        'diplo_locus' : 'tab:blue',
        'lls' : 'tab:orange',
        'wfabc' : 'tab:green',
        'bmws' : 'tab:red',
    }

    # one boxplot for each scenario
    # actually, one boxplot per {g160, g4000}
    for thisTotalGen in ['g160', 'g4000']:
        print (f'===== {thisTotalGen}')
        # keep order specified in scenarios
        theseScenarios = [x for x in scenarios if (x.strip().split('_')[0] == thisTotalGen)]

        # adjust limits appropriately
        if (thisTotalGen == 'g4000'):
            theYLims = (-0.01, 0.015)
        elif (thisTotalGen == 'g160'):
            theYLims = (-0.1, 0.15)
        else:
            assert (False), "Tertium non datur."


        # put all for this total time into one plot
        mosaicNames = ['A', 'B', 'C', 'D']
        Arrangement = [mosaicNames]
        mosaicFig, mosaicAx = plt.subplot_mosaic(Arrangement, figsize=(6.5, 3.8),
                                    gridspec_kw={'width_ratios': [4, 2, 3, 1]})
        # plt.rcParams.update({'font.size': 5})

        # seaborn.set_theme(style="whitegrid")
        # this is completely broken
        # seaborn.set_theme(style="whitegrid", font_scale=0.5)

        # add them all to the grid
        for (scIdx, thisScenario) in enumerate(theseScenarios):

            # which axis to put them on?
            thisAx = mosaicAx[mosaicNames[scIdx]]

            print (f'+++ {thisScenario}')
            scenarioSplits = thisScenario.strip().split('_')

            plotFrame = None
            # get the methods that have results for this scenario
            thisPalette = []
            numAboveYLimit = {}
            numBelowYLimit = {}
            plotOrder = []
            for testMethod in methodList:
                if (thisScenario in allMleS[testMethod]):

                    numAboveYLimit[testMethod] = {}
                    numBelowYLimit[testMethod] = {}
                    thisPalette.append (methodColors[testMethod])

                    # add the correct frame for this method to the frame for plotting
                    thisMLEFrame = allMleS[testMethod][thisScenario]

                    # how about the estimates
                    thisEstimates = thisMLEFrame['s2Hat'].to_numpy()
                    isNan = numpy.isnan (thisEstimates)

                    # how about the truth
                    getTrueS = (lambda idString : idString.strip().split('_')[1])
                    allTrueS = numpy.vectorize(getTrueS)(thisMLEFrame['id'].to_numpy())

                    # only up to 500 estiamtes
                    plotTruths = []
                    plotEstimates = []
                    maxNumPlot = 500
                    for thisTrueS in sToPlot[thisTotalGen]:
                        thisTrueSMask = (allTrueS == thisTrueS)
                        preEstimates = thisEstimates[~isNan & thisTrueSMask]
                        if (len(preEstimates) > maxNumPlot):
                            preEstimates = preEstimates[:maxNumPlot]
                        numAboveYLimit[testMethod][thisTrueS] = numpy.sum(preEstimates > theYLims[1])
                        numBelowYLimit[testMethod][thisTrueS] = numpy.sum(preEstimates < theYLims[0])
                        plotOrder.append ((testMethod, thisTrueS))
                        plotEstimates.extend (preEstimates)
                        plotTruths.extend ([thisTrueS]*len(preEstimates))

                    # make a frame
                    thisPlotFrame = pandas.DataFrame({
                        'method' : methodNames[testMethod],
                        'trueS' : plotTruths,
                        's2Hat' : plotEstimates,
                    })

                    # add to plotFrame
                    # hooray for pandas
                    if (plotFrame is None):
                        plotFrame = thisPlotFrame
                    else:
                        plotFrame = pandas.concat ([plotFrame, thisPlotFrame], ignore_index=True)

            # hooray for seaborn
            # put this earlier
            # plt.rcParams.update({'font.size': 7})
            # thisAx.rcParams.update({'font.size': 7})
            if ((scenarioSplits[1], scenarioSplits[2]) == ("h0.5", "const")):
                firstPlot = True
            else:
                firstPlot = False

            # numMethods = len(numpy.unique(plotFrame['method'].to_numpy()))
            # if (firstPlot):
            #     addBase = .5
            # else:
            #     addBase = 0
        
            # (fig, ax) = plt.subplots(figsize=(addBase+0.5+1.7*(numMethods/4), 3))
            # (fig, ax) = plt.subplots(figsize=(3,5))

            # make the boxplot
            # seaborn.boxplot (data=plotFrame, x='trueS', y='s2Hat', hue='method', ax=ax, fliersize=0.5, linewidth=0.5, palette=thisPalette) 
            seaborn.boxplot (data=plotFrame, x='trueS', y='s2Hat', hue='method', ax=thisAx, fliersize=0.5, linewidth=0.5, palette=thisPalette) 

            # add lones for truth
            for truth in sToPlot[thisTotalGen]:
                # ax.axline ((0,float(truth)), slope=0, color='red', ls='--', lw=0.5)
                thisAx.axline ((0,float(truth)), slope=0, color='red', ls='--', lw=0.5)
            # plt.ylim (theYLims)
            thisAx.set_ylim (theYLims)
            
            # set ticks sometimes
            yTicks = numpy.linspace (theYLims[0], theYLims[1], 6)
            theTruths = list(sToPlot[thisTotalGen])
            if (firstPlot):
                yTickLabels = [f"{x:.3f}" for x in yTicks]
                xTickLabels = theTruths
            else:
                yTickLabels = [""]*len(yTicks)
                xTickLabels = [" "]*len(theTruths)
            # plt.yticks (yTicks, yTickLabels)
            # plt.xticks (numpy.arange(len(theTruths)), xTickLabels)
            thisAx.set_yticks (yTicks, yTickLabels)
            thisAx.set_xticks (numpy.arange(len(theTruths)), xTickLabels)

            # see about putting number of outliers into plot
            # this should end up having same order as plotOrder
            xMeans = []
            # for c in ax.get_children():
            for c in thisAx.get_children():
                if type(c) == matplotlib.patches.PathPatch:
                    xMeans.append (0.5*c.get_extents().x0 + 0.5*c.get_extents().x1)
            assert (len(xMeans) == len(plotOrder)), (len(xMeans), len(plotOrder))
            # annotate things
            # this is horrible
            # inv = ax.transData.inverted()    
            inv = thisAx.transData.inverted()    
            # put in some dummies, so spacing is alway right
            genericParams = dict(arrowprops=dict(arrowstyle="->", color="#66666622"),
                                    bbox=dict(boxstyle="round", pad=0.125, edgecolor='#aaaaaa66', fc='#ffffff88'),
                                    ha='center', va='center', color='black', annotation_clip=True)
            belowFactor = 0.28
            aboveFactor = 0.19
            # below
            # dummy = ax.annotate(text=str(42), xy=(0, theYLims[0]),
            #     xytext=(0, theYLims[0] - belowFactor*(theYLims[1] - theYLims[0])))
            dummy = thisAx.annotate(text=str(42), xy=(0, theYLims[0]),
                xytext=(0, theYLims[0] - belowFactor*(theYLims[1] - theYLims[0])))
            dummy.set_alpha(0)
            # above
            # dummy = ax.annotate(text=str(42), xy=(0, theYLims[1]),
            #     xytext=(0, theYLims[1] + aboveFactor*(theYLims[1] - theYLims[0])))
            dummy = thisAx.annotate(text=str(42), xy=(0, theYLims[1]),
                xytext=(0, theYLims[1] + aboveFactor*(theYLims[1] - theYLims[0])))
            dummy.set_alpha(0)
            # and more, if needed      
            for (aIdx, thisX) in enumerate(xMeans):
                # that's how we get it =(
                numBelow = numBelowYLimit[plotOrder[aIdx][0]][plotOrder[aIdx][1]]
                numAbove = numAboveYLimit[plotOrder[aIdx][0]][plotOrder[aIdx][1]]
                (realX, dummy) = inv.transform((thisX, 0))
                if (numBelow > 0):
                    # ax.annotate(text=str(numBelow), xy=(realX, theYLims[0]),
                    #     xytext=(realX, theYLims[0] - belowFactor*(theYLims[1] - theYLims[0])), **genericParams)
                    thisAx.annotate(text=str(numBelow), xy=(realX, theYLims[0]),
                        xytext=(realX, theYLims[0] - belowFactor*(theYLims[1] - theYLims[0])), **genericParams)
                if (numAbove > 0):
                    # ax.annotate(text=str(numAbove), xy=(realX, theYLims[1]),
                    #     xytext=(realX, theYLims[1] + aboveFactor*(theYLims[1] - theYLims[0])), **genericParams)
                    thisAx.annotate(text=str(numAbove), xy=(realX, theYLims[1]),
                        xytext=(realX, theYLims[1] + aboveFactor*(theYLims[1] - theYLims[0])), **genericParams)

            # random stuff
            # plt.ylabel ("")
            # plt.xlabel (r"True $s_{AA}$", fontsize=9)
            thisAx.set_ylabel ("")
            if (firstPlot):
                # thisAx.set_xlabel (r"True $s_{AA}$", fontsize=basicFontSize)
                thisAx.set_xlabel (r"True $s_{AA}$")
            else:
                thisAx.set_xlabel ("")
            theTitle = ""
            if (scenarioSplits[1] == 'h1.0'):
                theTitle += r'$h=1$'
            elif (scenarioSplits[1] == 'h0.5'):
                theTitle += r'$h=\frac{1}{2}$'
            else:
                assert (False)
            # plt.title (theTitle + f" ({scenarioSplits[2][0]})")
            # plt.tight_layout()
            thisAx.set_title (theTitle + f" ({scenarioSplits[2][0]})")
            # thisAx.set_title (theTitle + f" ({scenarioSplits[2][0]})", fontsize=basicFontSize)
            # comes ater when savng the file
            # thisAx.tight_layout()
            if (firstPlot):
                thisAx.legend().set_title("")

            else:
                # ax.get_legend().remove()
                thisAx.get_legend().remove()

            # # save it to file
            # plt.savefig (f"box_{thisScenario}.pdf")

            # break

        mosaicFig.tight_layout()
        # save it to file
        mosaicFig.savefig (f"box_{thisTotalGen}.pdf")


def load_times (slurmifyDir, slurmifySkeleton, allMleS):

    # see what files we have
    runTimes = {}
    benchmarkTime = None
    for batchName in sorted(list(pathlib.Path(slurmifyDir).glob(slurmifySkeleton % '*'))):

        print (f"+++++ {batchName}")


        # get info about method and scenario from sbatch script
        # read the script
        ifs = open (batchName, 'r')
        readLines = ifs.readlines()
        ifs.close ()

        firstPreCmdIdx = numpy.where(["PRE_CMD" in x for x in readLines])[0][0]
        cmdString = readLines[firstPreCmdIdx + 3]
        cmdPieces = cmdString.strip().split()
        # print (cmdPieces)
        if (cmdPieces[0] == 'Rscript'):
            # LLS
            thisMethod = 'lls'
            # get the scenario from the input file
            thisScenario = "_".join(pathlib.Path(cmdPieces[1]).stem.split('_')[:3])
        elif ((cmdPieces[0] == 'bmws') or ((cmdPieces[0] == 'python') and (cmdPieces[1] == 'bmws/bmws_change.py'))):
            # bmws
            thisMethod = 'bmws'
            # get the scenario from the input file
            offset = int(cmdPieces[0] == 'python')
            thisScenario = pathlib.Path(cmdPieces[2 + offset]).stem
        elif (cmdPieces[0] == 'DiploLocus-likelihood'):
            # diplo-locus
            thisMethod = 'diplo_locus'
            # get the scenario from the input file
            thisScenario = "_".join(pathlib.Path(cmdPieces[2]).stem.split('_')[:3])
        elif (cmdPieces[0] == './WFABC_v1.1/binaries/Linux/wfabc_1'):
            # WFABC
            thisMethod = 'wfabc'
            # get the scenario from the input file
            thisScenario = "_".join(pathlib.Path(cmdPieces[7]).stem.split('_')[:3])
        else:
            print (cmdPieces)
            assert (False), f'Unkown method: {cmdPieces[0]}'
        print (thisMethod, thisScenario)


        # get info on runtime from the out-file
        outputFileTest = list(pathlib.Path(batchName.parent).glob(f"{batchName.stem}_j*.out"))
        assert (len(outputFileTest) == 1), "Multiple *.out files for one slurm script is not allowed, because it is ambiguous which ones is the correct one."
        outFile = outputFileTest[0]
        ifs = open (outFile, 'r')
        outLines = ifs.readlines()
        ifs.close()

        # where are the times?
        # print (outLines)
        preBenchLines = numpy.where(["PRE_BENCH" in x for x in outLines])[0]
        # print (preBenchLines)
        assert (len(preBenchLines) == 2), "Unexpected output in file."
        preCmdLines = numpy.where(["PRE_CMD" in x for x in outLines])[0]
        # print (preCmdLines)
        assert (len(preCmdLines) == 2), "Benchmark command did not run."
        postCmdLines = numpy.where(["POST_CMD" in x for x in outLines])[0]
        # print (postCmdLines)
        assert (len(postCmdLines) <= 2), '"POST_CMD" occurs too often in outfile'
        runFinished = (len(postCmdLines) == 2)
        if (runFinished):
            assert (allMleS[thisMethod][thisScenario].shape[0] == NUM_REPLICATES), "Results missing."

        # get the times
        preBenchTime = float (outLines[preBenchLines[0]+1].strip())
        preCmdTime = float (outLines[preCmdLines[0]+1].strip())
        thisBenchTime = preCmdTime - preBenchTime
        # store ONE benchmarkTime, if we haven't already
        if (benchmarkTime is None):
            benchmarkTime = thisBenchTime
        print (thisBenchTime)

        # for the actual runtime, see what happend
        print (allMleS[thisMethod][thisScenario].shape)
        if (runFinished):
            postCmdTime = float (outLines[postCmdLines[0]+1].strip())
            thisCmdTime = postCmdTime - preCmdTime
        else:
            # actually did not finish, so estimate how long it would take from how much has been run
            elapsedTime = MAX_TIME - thisBenchTime
            numFinished = allMleS[thisMethod][thisScenario].shape[0]
            thisCmdTime = (NUM_REPLICATES/numFinished) * elapsedTime

        # now, normalize the time via the benchmark
        adjustedRuntime = benchmarkTime * (thisCmdTime/thisBenchTime)
        print (adjustedRuntime)

        # and store it
        if (thisMethod not in runTimes):
            runTimes[thisMethod] = {}
        runTimes[thisMethod][thisScenario] = adjustedRuntime


    # and return it
    return runTimes


def makeTable (tableFilename, runTimes):
    # apparently, pandas make tables
    methodList = ['diplo_locus', 'lls', 'wfabc', 'bmws']
    methodNames = {
        'diplo_locus' : '\\texttt{diplo-locus}',
        'lls' : '\\texttt{LLS}',
        'wfabc' : '\\texttt{WFABC}',
        'bmws' : '\\texttt{bmws}',
    }
    timeFrame = pandas.DataFrame(columns=['\# Gener.', '$h$', 'Var. $s$'] + [methodNames[x] for x in methodList])
    for totalGensKey in ["g160", "g4000"]:
        for fixedH in ["0.5", "1.0"]:
            for toChange in ["const", "var"]:
                thisScenarioKey = f"{totalGensKey}_h{fixedH}_{toChange}"
                rowList = [totalGensKey[1:], fixedH, "No" if (toChange == 'const') else "Yes"]
                for methodName in methodList:
                    if (thisScenarioKey in runTimes[methodName]):
                        # write it in hours
                        rowList.append (runTimes[methodName][thisScenarioKey] / 3600)
                    else:
                        rowList.append ('-')
                timeFrame.loc[len(timeFrame.index)] = rowList
    timeFile = tableFilename
    ofs = open (timeFile, 'w')
    ofs.write (timeFrame.to_latex (index=False, float_format="%.3f", column_format='cccrrrr'))
    ofs.close()


def main():

    # which scenarios do we have -- hopefully the same for all methods
    # this also determines the boxplot order
    scenarios = []
    for totalGensKey in ["g4000", "g160"]:
        for fixedH in ["0.5", "1.0"]:
            for toChange in ["const", "var"]:
                scenarios.append (f"{totalGensKey}_h{fixedH}_{toChange}")

    # get the mles
    allMleS = {}
    print ("[LOAD_DIPLO_LOCUS]")
    diploString = 'diplo_locus'
    allMleS[diploString] = load_diplo_locus_mles (RESULTS_DIR[diploString], scenarios)
    print ("[LOAD_DIPLO_LOCUS_DONE]")
    print ("[LOAD_LLS]")
    llsString = 'lls'
    allMleS[llsString] = load_lls_mles (RESULTS_DIR[llsString], scenarios)
    print ("[LOAD_LLS_DONE]")
    print ("[LOAD_WFABC]")
    wfabcString = 'wfabc'
    allMleS[wfabcString] = load_wfabc_mles (RESULTS_DIR[wfabcString], scenarios)
    print ("[LOAD_WFABC_DONE]")
    print ("[LOAD_BMWS]")
    bmwsString = 'bmws'
    allMleS[bmwsString] = load_bmws_mles (RESULTS_DIR[bmwsString], scenarios)
    print ("[LOAD_BMWS_DONE]")


    # how about nans
    nanFraction = {}
    for (thisMethod, thisMles) in allMleS.items():
        # print (f"+++++ {thisMethod}")
        nanFraction[thisMethod] = {}
        for (thisScenario, thisEstimates) in thisMles.items():
            # print (f"{thisScenario}: {thisEstimates.shape}")
            numEstimates = thisEstimates.shape[0]
            print (thisMethod, thisScenario, numEstimates)
            numNanEstimates = numpy.sum(numpy.isnan(thisEstimates['s2Hat']))
            nanFraction[thisMethod][thisScenario] = numNanEstimates/numEstimates


    # get the times
    print ("[LOAD_TIMES]")
    runTimes = load_times (SLURMIFY_DIR, SLURMIFY_SKELETON, allMleS)
    print ("[LOAD_TIMES_DONE]")

    for (thisMethod, thisTimes) in runTimes.items():
        print (f"+++++ {thisMethod}")
        for (thisScenario, thisTime) in thisTimes.items():
            floatString = f"{thisTime:12.3f}"
            scenarioString = f"{thisScenario}"
            bufferString = " "*(30-len(scenarioString)-len(floatString))
            print (f"{scenarioString}:{bufferString}{floatString}\t\t{nanFraction[thisMethod][thisScenario]:2.3f}")


    # write times in a nice table for latex =)
    print ("[TABLE]")
    makeTable ('time_table.tex', runTimes)
    print ("[TABLE_DONE]")


    # make some boxplot
    print ("[BOXPLOTS]")
    theBoxPlots (scenarios, allMleS, nanFraction)
    print ("[BOXPLOTS_DONE]")


if __name__ == "__main__":
    main()
