import numpy
import pandas
import seaborn
import scipy
import matplotlib.pyplot as plt
import pathlib
import pickle
import bz2




RESULTS_PICKLE_FILE = pathlib.Path ("final_lls.pickle.bz2")
PLOT_FILE = pathlib.Path ("QQ_ROC_BOX_onset.pdf")




def loadPickledResults (pickleFile):

    # unpickle
    pickleData = pickle.load (bz2.open (pickleFile, "rb"))
    numSimReps = pickleData['numSimReps']
    constMLEs = pickleData['constMLEs']
    onsetMLEs = pickleData['onsetMLEs']

    # then build a dataframes for seaborn
    mles = {
        'onset' : onsetMLEs,
        'constant' : constMLEs,
    }
    frames = {}
    for (thisSimScenario, thisMLEs) in mles.items():
        thisColumns = ['onsetScenario', 'trueS2', 'rep', 's2Hat']
        if (thisSimScenario == 'onset'):
            thisColumns.append ('tHat')
        thisFrame = None
        # fill it
        for (k,v) in thisMLEs.items():
            (thisScenario, thisS2String) = k
            thisDict = {
                'onsetScenario' : numpy.repeat (thisScenario, numSimReps),
                'trueS2' : numpy.repeat (float(thisS2String), numSimReps),
                'rep' : numpy.arange (numSimReps),
                'baseLL' : v[:,0],
                'maxLL' : v[:,1],
                's2Hat' : v[:,2],
            }
            if (thisSimScenario == 'onset'):
                thisDict['tHat'] = v[:,3]

            # and extend the big frame
            thisDF = pandas.DataFrame (thisDict)
            if (thisFrame is None):
                thisFrame = thisDF
            else:
                thisFrame = pandas.concat([thisFrame, thisDF], axis=0)
        
        # store the new frame
        frames[thisSimScenario] = thisFrame

    # then put the frame with the onset results together
    preFrame = frames['onset']
    onsetFrame = preFrame[preFrame['onsetScenario'] == 'onset']

    return onsetFrame


def boxplotOnset (ax, theTitle, onsetFrame):
    ax = seaborn.boxplot(data=onsetFrame, x="trueS2", y="tHat", fliersize=0.5, orient="v", ax=ax)

    ax.grid(visible=True, axis='y', lw=0.3)
    # add aux lines
    trueT = 2000
    ax.axline ((0,trueT), slope=0, color='C3', linestyle='dashed', lw=0.5)
    ax.set_title (theTitle, loc='left', fontsize='medium')
    ax.tick_params(axis='both', labelsize="small")
    ax.set_xlabel(r'True $s_{AA}$', fontsize="small")
    ax.set_ylabel(r'Inferred $\hat{t}_{o}$', fontsize="small")

    return ax


def boxplotSelCoeff (ax, theTitle, onsetFrame):
    ax = seaborn.boxplot(data=onsetFrame, x="trueS2", y="s2Hat", fliersize=0.5, orient="v", ax=ax)
    theYLims = (-0.01, 0.015)
    ax.set_ylim (theYLims[0], theYLims[1])

    # see about outliers =(
    # first get the numbrs
    simS2 = numpy.unique(onsetFrame["trueS2"])
    numAboveYLimit = -1 * numpy.ones(len(simS2))
    numBelowYLimit = -1 * numpy.ones(len(simS2))
    for (sIdx, thisS) in enumerate(simS2):
        thisSFrame = onsetFrame[numpy.isclose (onsetFrame['trueS2'], thisS)]
        numAboveYLimit[sIdx] = numpy.sum(thisSFrame['s2Hat'] > theYLims[1])
        numBelowYLimit[sIdx] = numpy.sum(thisSFrame['s2Hat'] < theYLims[0])

    # put the outliers in the plot
    # PLEASE NOTE: this is horrible, but that's what seaborn wants
    # put in some dummies, so spacing is alway right
    genericParams = dict(arrowprops=dict(arrowstyle="->", color="#66666622"),
                            bbox=dict(boxstyle="round", pad=0.125, edgecolor='#aaaaaa66', fc='#ffffff88'),
                            ha='center', va='center', color='black', annotation_clip=True, fontsize="small")
    belowFactor = 0.35
    aboveFactor = 0.19
    yRange = theYLims[1] - theYLims[0]
    # below
    if (numBelowYLimit.sum() > 0):
        dummy = ax.annotate(text=str(42), xy=(0, theYLims[0]), xytext=(0, theYLims[0] - belowFactor*yRange))
        dummy.set_alpha(0)
    # above
    if (numAboveYLimit.sum() > 0):
        dummy = ax.annotate(text=str(42), xy=(0, theYLims[1]), xytext=(0, theYLims[1] + aboveFactor*yRange))
        dummy.set_alpha(0)
    # and more, if needed      
    for (sIdx, thisS) in enumerate(simS2):
        # that's how we get it =(
        numBelow = numBelowYLimit[sIdx]
        numAbove = numAboveYLimit[sIdx]
        if (numBelow > 0):
            ax.annotate(text=str(int(numBelow)), xy=(sIdx, theYLims[0]), xytext=(sIdx, theYLims[0] - belowFactor*yRange), **genericParams)
        if (numAbove > 0):
            ax.annotate(text=str(int(numAbove)), xy=(sIdx, theYLims[1]), xytext=(sIdx, theYLims[1] + aboveFactor*yRange), **genericParams)

    # back to normal formatting
    ax.grid(visible=True, axis='y', lw=0.3)
    # add aux lines
    left_x, right_x = ax.get_xlim()
    ax.hlines(y=list(map(float, simS2)), xmin=left_x, xmax=right_x, colors=['black'] * len(simS2), linestyle='dashed', lw=0.5)
    ax.tick_params(axis='both', labelsize="small")
    ax.set_title (theTitle, loc='left', fontsize='medium')
    ax.set_xlabel(r'True $s_{AA}$', fontsize="small")
    ax.set_ylabel(r'Inferred $\hat{s}_{AA}$', fontsize="small")
    
    return ax
    
    
def qqPlotOnset (ax, theTitle, onsetFrame):

    neutralFrame = onsetFrame[numpy.isclose(onsetFrame['trueS2'], 0)]
    neutralLLRstats = 2 * (neutralFrame['maxLL'] - neutralFrame['baseLL'])
    pValuesChiOne = sorted(1 - scipy.stats.chi2.cdf (neutralLLRstats, 1))
    logPChiOne = -numpy.log10(pValuesChiOne)
    pValuesChiTwo = sorted(1 - scipy.stats.chi2.cdf (neutralLLRstats, 2))
    logPChiTwo = -numpy.log10(pValuesChiTwo)
    expPValues = numpy.arange (1, len(neutralLLRstats) + 1)
    expPValues = expPValues / (len(neutralLLRstats) + 1)
    logPExp = -numpy.log10(expPValues)

    theMax = max (logPChiOne.max(), logPChiTwo.max(), logPExp.max())

    theLegend = []
    ax.plot (logPExp, logPChiOne, '.')
    theLegend.append (r"$\chi^2(1)$")
    ax.plot (logPExp, logPChiTwo, '.')
    theLegend.append (r"$\chi^2(2)$")
    ax.axline ((0,0), slope=1, color='C3')

    # some formatting
    ax.legend (theLegend, loc='best', fontsize="x-small")
    ax.set_title (theTitle, loc='left', fontsize='medium')
    ax.tick_params(axis='both', labelsize="small")
    ax.set_xlabel(r'$-\log_{10}(\text{quantile})$', fontsize="small")
    ax.set_ylabel(r'$-\log_{10}(p)$', fontsize="small")
    ax.set_xlim(-0.1, theMax + 0.5)
    ax.set_ylim(-0.1, theMax + 0.5)
    ax.axhline(0, 0, 1, ls=':', c='darkgray', lw=1)
    ax.axvline(0, 0, 1, ls=':', c='darkgray', lw=1)

    return ax


def getTFPNstats (neutralLLR, selLLR, numThs, smallestNonZero=0.1):
    # make sure things are sorted, default is ascending
    assert (len(neutralLLR) == len(selLLR))
    neutralLLR = numpy.sort(neutralLLR)
    selLLR = numpy.sort(selLLR)
    # order is: TP, FP, TN, FN
    TFPNstats = -1 * numpy.ones((numThs+1,4))
    # get some ths
    maxLLR = max(neutralLLR.max(), selLLR.max())
    ths = numpy.concatenate ([[0], numpy.geomspace (smallestNonZero, maxLLR + 0.1, numThs)])
    
    neutralBelowTH = numpy.searchsorted(neutralLLR, ths)/len(neutralLLR)
    selBelowTH = numpy.searchsorted(selLLR, ths)/len(selLLR)
    # TP
    TFPNstats[:,0] = 1 - selBelowTH
    # FP
    TFPNstats[:,1] = 1 - neutralBelowTH
    # TN
    TFPNstats[:,2] = neutralBelowTH
    # FN
    TFPNstats[:,3] = selBelowTH

    # give it away
    assert (TFPNstats.min() >= 0)
    return TFPNstats


def rocOnset (ax, theTitle, onsetFrame):

    simS2 = numpy.unique(onsetFrame["trueS2"])
    line_colors = plt.cm.Dark2(range(len(simS2)))

    neutralFrame = onsetFrame[numpy.isclose(onsetFrame['trueS2'], 0)]
    neutralLLR = numpy.sort(2 * (neutralFrame['maxLL'] - neutralFrame['baseLL']))

    theLegend = []
    for (s2Idx, thisS2) in enumerate(simS2):
        thisSelFrame = onsetFrame[numpy.isclose(onsetFrame['trueS2'], thisS2)]
        thisSelLLR = numpy.sort(2 * (thisSelFrame['maxLL'] - thisSelFrame['baseLL']))
        TFPN = getTFPNstats (neutralLLR, thisSelLLR, 100)
        TPR = TFPN[:,0] / (TFPN[:,0] + TFPN[:,3])
        FPR = TFPN[:,1] / (TFPN[:,1] + TFPN[:,2])
        ax.plot(FPR, TPR, color=line_colors[s2Idx], label=str(s2Idx))
        theLegend.append (thisS2)

    # title and legend
    ax.set_title (theTitle, loc='left', fontsize="medium")
    ax.legend (theLegend, loc="best", edgecolor='None', facecolor='None',
               fontsize="x-small", title_fontsize="small", frameon=False, title="True $s_{AA}$")

    # some more cosmetics
    ax.set_xlabel('False Positive Rate', fontsize="small")
    ax.set_ylabel('True Positive Rate', fontsize="small")
    ax.set_xlim(-0.04, 1.04)
    ax.set_ylim(-0.04, 1.04)
    ax.tick_params(labelsize="small")
    ax.axhline(0, 0, 1, ls=':', c='darkgray', lw=1)
    ax.axvline(0, 0, 1, ls=':', c='darkgray', lw=1)
    ax.axhline(1, 0, 1, ls=':', c='darkgray', lw=1)
    ax.axvline(1, 0, 1, ls=':', c='darkgray', lw=1)
    ax.axline((0, 0), (1, 1), ls=":", lw=0.5, c="0.5")

    return ax


def combinedPlot (onsetFrame, plotFile):

    # all together now
    Arrangement = [['A', 'B'],
                ['C', 'D']]
    fig, ax = plt.subplot_mosaic(Arrangement, figsize=(6.5, 5), gridspec_kw={'width_ratios': [2, 4], 'height_ratios': [1, 1]})

    # panel A: QQ-plot
    ax['A'] = qqPlotOnset (ax['A'], r"$\mathbf{A)}$ QQ-plot", onsetFrame)

    # panel B: boxplot selection
    ax['B'] = boxplotSelCoeff (ax['B'], r"$\mathbf{B)}$ Inferred selection coefficient", onsetFrame)

    # panel C: ROC
    ax['C'] = rocOnset (ax['C'], r"$\mathbf{C)}$ ROC curves", onsetFrame)

    # panel B: boxplot onset
    ax['D'] = boxplotOnset (ax['D'], r"$\mathbf{D)}$ Onset of selection", onsetFrame)

    plt.subplots_adjust(left=0.1,
                        bottom=0.1, 
                        right=0.975, 
                        top=0.9, 
                        wspace=0.4, 
                        hspace=0.6)

    fig.savefig(plotFile)


def main():

    print ("[LOAD_RESULTS]")
    onsetFrame = loadPickledResults (RESULTS_PICKLE_FILE)
    print ("[LOAD_RESULTS_DONE]")

    print ("[COMBINED_PLOT]")
    combinedPlot (onsetFrame, PLOT_FILE)
    print ("[COMBINED_PLOT_DONE]")


if __name__ == "__main__":
    main()

