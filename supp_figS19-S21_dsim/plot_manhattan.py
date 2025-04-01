import pathlib
import pickle
import bz2
import matplotlib.pyplot as plt
import numpy




OUTPUT_DIR = 'results/'
FISHER_FILE = pathlib.Path(OUTPUT_DIR, "analyzed_F0-F60_fisher.pickle.bz2")
MANHATTAN_FILE = "dsim_manhattan_fisher.pdf"




def loadFisher():
    # pickle in a jar
    pickleData = pickle.load (bz2.open (FISHER_FILE, "rb"))
    commonPositions = pickleData["commonPositions"]
    fisherPValues = pickleData["fisherPValues"]
    return (commonPositions, fisherPValues)


def plotFisher (commonPositions, fisherPValues):

    # Like this? Seems to be a global "scaling"
    plt.rcParams['font.size'] = 10

    (fig, ax) = plt.subplots (figsize=(6.5, 4))

    # go through all chromosomes we have
    plotOffset = 0
    xTicks = []
    XTickLabels = []
    numFisherTests = 0
    # for now only 1 chromosome
    for thisChr in ['chr2L']:
        print (f"-- {thisChr}")
        plotPositions = commonPositions[thisChr]
        logPValues = -numpy.log10(fisherPValues[thisChr])
        assert (len(plotPositions) == len(logPValues)), (len(plotPositions), len(logPValues))
        ax.plot ((plotOffset + plotPositions)/1e6, logPValues, 'o', ms=0.5, rasterized=True)
        # xTicks.append (plotOffset + 0.5*(plotPositions.min() + plotPositions.max()))
        # XTickLabels.append (thisChr)
        plotOffset = plotOffset + plotPositions.max()
        numFisherTests += len(plotPositions)

    # ax.set_xticks (xTicks, XTickLabels)
    ax.set_xticks ([0, 5, 10, 15, 20])
    ax.axline ((0,-numpy.log10(0.05/numFisherTests)), slope=0, ls='--', color='tab:orange')
    ax.axline ((0,0), slope=0, ls='-', lw=0.5, color='black')
    ax.set_title (f"")
    ax.set_xlabel("Position on chromosome 2L (Mbp)")
    ax.set_ylabel(r'$-\log_{10}p_{\chi^2(1)}$')
    ax.set_ylim ((-0.5,18))

    fig.tight_layout()
    fig.savefig (MANHATTAN_FILE, dpi=600)
    # plt.clf()


def plotManhattan():

    # load the data
    print ("[LOAD_FISHER]")
    (commonPositions, fisherPValues) = loadFisher()
    print ("[LOAD_FISHER_DONE]")

    # plot it
    print ("[PLOT_FISHER]")
    plotFisher (commonPositions, fisherPValues)
    print ("[PLOT_FISHER_DONE]")


def main():
    plotManhattan()


if __name__ == "__main__":
    main()
