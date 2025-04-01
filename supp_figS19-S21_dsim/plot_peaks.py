import pathlib
import pickle
import bz2
import numpy
import pandas
import matplotlib.pyplot as plt




# this file contains the p-values combined across biological replicates using Fisher's mwthod
FISHER_PVALUES_FILE = pathlib.Path ("results/analyzed_F0-F60_fisher.pickle.bz2")
# this file contains the annotations for build M252
M252_ORIGINAL_ANNOTATIONS_FILE = pathlib.Path ("dsim-M252-popgen-ann-r1.1.gtf")
# store processed annotations here
M252_PROCESSED_ANNOTATIONS_FILE = pathlib.Path ("chr2L_processed_annotations.pickle.bz2")
# here are the real names for the relevant genes
FLYBASE_GENE_NAMES_FILE = pathlib.Path ("flybase_gene_names.txt")

# plot regions on this chromosome
PLOT_CHR = 'chr2L'

# these are the two regions to plot
PLOT_REGIONS = [
    (3.165e6,3.2204e6),
    (6.065e6, 6.134e6),
]
# these are the filenames for the plots
PLOT_FILENAMES = [
    "dsim_peak1.pdf",
    "dsim_peak2.pdf",
]




def loadFisherPValues():
    
    # this is the file with the p-values
    pickleData = pickle.load (bz2.open (FISHER_PVALUES_FILE, 'r'))

    # here are the p-values
    pos = pickleData['commonPositions'][PLOT_CHR]
    pValues = pickleData['fisherPValues'][PLOT_CHR]
    assert (len(pos) == len(pValues))
    logPValues = -numpy.log10(pValues)
    bfThres = -numpy.log10(0.05/len(pValues))

    return (pos, logPValues, bfThres)


def loadM252Annotations():

    # which chromosome do we want to work with?
    thisChr = PLOT_CHR.split('chr')[1]

    # only do annotations processing if processed file does not exist
    if (not M252_PROCESSED_ANNOTATIONS_FILE.is_file()):
        print ("[PROCESSED] File does not exist, so process original data.")
        # process the original annotations
        # here are the annotations
        dsimFrame = pandas.read_csv (str(M252_ORIGINAL_ANNOTATIONS_FILE), sep='\t', header=None, dtype={7: str})
        # this is from https://www.ensembl.org/info/website/upload/gff.html
        dsimFrame.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

        # only specified chromosome
        print (dsimFrame.shape)
        dsimFullAnnoFrame = dsimFrame[dsimFrame["seqname"] == thisChr]
        dsimFullAnnoFrame.reset_index (inplace=True, drop=True)
        print (dsimFullAnnoFrame.shape)

        # unravel the attributes
        # this is slow, but I don't think we have a choice
        splitAttributes = numpy.char.split (numpy.array(dsimFullAnnoFrame["attribute"], dtype=str), ";")
        # go through attribute lists
        for (thisIdx, thisList) in enumerate(splitAttributes):
            # go through attributes in this list
            for thisAttribute in thisList:
                attributeString = thisAttribute.strip()
                # only non-empty ones
                if (attributeString != ''):
                    (k, v) = attributeString.split(" ")
                    # if column doesn't exist, add it
                    if (k not in dsimFullAnnoFrame):
                        dsimFullAnnoFrame.loc[:, [k]] = None
                    # add the value
                    dsimFullAnnoFrame.loc[thisIdx,k] = v.strip('" ')
        
        # get rid of the column
        dsimFullAnnoFrame = dsimFullAnnoFrame.drop(['attribute'],axis=1)

        # and save
        pickleData = {
            "dsimFullAnnoFrame" : dsimFullAnnoFrame,
        }
        pickle.dump (pickleData, bz2.open (M252_PROCESSED_ANNOTATIONS_FILE, 'w'))
    else:
        print ("[PROCESSED] File exist, load it.")
        # processed file already exists, load it
        pickleData = pickle.load (bz2.open (M252_PROCESSED_ANNOTATIONS_FILE, 'r'))
        dsimFullAnnoFrame = pickleData["dsimFullAnnoFrame"]

    # however you got it, return it
    return dsimFullAnnoFrame


def plotRegionsGeneTrack (pos, logPValues, bfThres, dsimFullAnnoFrame, flybaseToGene):

    print ("[COLLECT_GENE_COORDS]")
    # it helps to have a reduced frame just for the gene coordinates
    # reduce info to just a list of genes with min and max pos of any of its features
    uniqueGenes = numpy.unique (dsimFullAnnoFrame["gene_id"])
    starts = []
    ends = []
    for thisGene in uniqueGenes:
        geneFrame = dsimFullAnnoFrame[dsimFullAnnoFrame["gene_id"] == thisGene]
        starts.append (geneFrame["start"].min())
        ends.append (geneFrame["end"].max())
    geneCoordFrame = pandas.DataFrame ({
        'gene_id' : uniqueGenes,
        'start' : starts,
        'end' : ends,
    })
    print ("[COLLECT_GENE_COORDS_DONE]")

    # parameters for ticks and labels
    xTicks = [
        [x*1e6 for x in [3.17, 3.18, 3.19, 3.20, 3.21, 3.22]],
        [x*1e6 for x in [6.07, 6.08, 6.09, 6.10, 6.11, 6.12, 6.13]],
    ]

    # set "global" fontsize
    plt.rcParams['font.size'] = 9

    # need to know who exceeds thresholds
    exceedMask = (logPValues > bfThres)

    for (rIdx, plotRegion) in enumerate(PLOT_REGIONS):

        print ("new peak:")

        # use only pValues in the given region
        plotMask = (plotRegion[0] < pos) & (pos < plotRegion[1])
        # whoo does and does not exceed in this region?
        plotExceedMask = plotMask & exceedMask
        plotNotExceedMask = plotMask & ~exceedMask

        # plot p-values
        (fig, ax) = plt.subplots (figsize=(6, 3.5))
        ax.plot (pos[plotNotExceedMask], logPValues[plotNotExceedMask], 'o', ms=1.5, color='C0')
        ax.plot (pos[plotExceedMask], logPValues[plotExceedMask], 'o', ms=1.5, color='C1')
        ax.axhline (y=bfThres, color='C1', ls='--', lw=0.5)
        ax.axhline (y=0, color='black', ls='-')

        # ticks
        ax.set_xticks (xTicks[rIdx], labels=[x/1e6 for x in xTicks[rIdx]])
        logPTicks = [0,5,10,15]
        ax.set_yticks (logPTicks, labels=logPTicks)

        # plot some vertical guides
        plotExceedPos = pos[plotExceedMask]
        minExceed = min(plotExceedPos)
        maxExceed = max(plotExceedPos)
        ax.axvline (x=minExceed, color='black', ls='--', lw=0.5)
        ax.axvline (x=maxExceed, color='black', ls='--', lw=0.5)
        print ((min(plotExceedPos), max(plotExceedPos)))

        # plot gene track
        trackY = -2
        trackText = 0.5
        offset = 1
        for (gIdx, thisGene) in geneCoordFrame.iterrows():
            # any gene that overlaps wiith our plotregion gets plotted
            if (((plotRegion[0] < thisGene["start"]) and (thisGene["start"] < plotRegion[1])) or
                ((plotRegion[0] < thisGene["end"]) and (thisGene["end"] < plotRegion[1]))):

                # but we don't plot it if it has no name
                thisGeneId = thisGene["gene_id"]
                realGeneName = thisGeneId
                if (thisGeneId in flybaseToGene):
                    realGeneName = flybaseToGene[thisGeneId]
                    if (realGeneName == "na"):
                        continue

                yPos = trackY * offset
                offset += 1
                trackStart = max(thisGene["start"], plotRegion[0])
                trackEnd = min(thisGene["end"], plotRegion[1])
                ax.plot ([trackStart,trackEnd], [yPos, yPos], ls='-', color='C2')
                # plot beginning or end marker if whole gene in sight
                if (numpy.isclose (trackStart,thisGene["start"])):
                    ax.plot ([trackStart], [yPos], marker='.', color='C2')
                if (numpy.isclose (trackEnd,thisGene["end"])):
                    ax.plot ([trackEnd], [yPos], marker='.', color='C2')

                # and maybe some exons?
                thisGeneFrame = dsimFullAnnoFrame[dsimFullAnnoFrame['gene_id'] == thisGeneId]
                exonOnlyFrame = thisGeneFrame[thisGeneFrame['feature'] == 'exon']
                for thisExon in list(zip (exonOnlyFrame['start'], exonOnlyFrame['end'])):
                    ax.plot ([thisExon[0],thisExon[1]], [yPos, yPos], ls='-', lw=7, color='C2', solid_capstyle="butt")

                # put label of gene in plot
                ax.text (trackStart, yPos + trackText, realGeneName, clip_on=True)

                # print gene id if within +/- 10kb of significant region
                printRegion = (minExceed  - 1e4, maxExceed + 1e4)
                if (((printRegion[0] < thisGene["start"]) and (thisGene["start"] < printRegion[1])) or
                    ((printRegion[0] < thisGene["end"]) and (thisGene["end"] < printRegion[1]))):
                    # print (thisGeneId, (thisGene["start"], thisGene["end"]))
                    print (f"{thisGeneId} -> {realGeneName}")

        # labels and limits
        ax.set_xlim (plotRegion[0], plotRegion[1])
        ax.set_xlabel ("Position on Chromosome 2L (Mbp)")
        ax.set_ylabel (r"$-\log_{10}(p)$", loc="top")

        # save plot
        fig.tight_layout()
        fig.savefig (PLOT_FILENAMES[rIdx])


def loadFlybaseGeneNames():

    # load the file
    nameFrame = pandas.read_csv (str(FLYBASE_GENE_NAMES_FILE), sep='\t', header=None, index_col=False)

    # build a dictionary
    flybaseToGene = {}
    for nameRow in nameFrame.iterrows():
        flybaseToGene[str(nameRow[1][0])] = str(nameRow[1][1])

    # and return it
    # print (flybaseToGene)
    return flybaseToGene


def main():

    # load p-values
    print ("[LOAD_P_VALUES]")
    (positions, logPValues, bfThres) = loadFisherPValues()
    print ("[LOAD_P_VALUES_DONE]")

    # load annotation
    print ("[LOAD_ANNOTATIONS]")
    dsimFullAnnoFrame = loadM252Annotations()
    print ("[LOAD_ANNOTATIONS_DONE]")

    # load real gene names
    print ("[LOAD_REAL_NAMES]")
    flybaseToGene = loadFlybaseGeneNames()
    print ("[LOAD_REAL_NAMES_DONE]")

    # make the zoomed plots for both peaks
    print ("[PLOT_REGIONS]")
    plotRegionsGeneTrack (positions, logPValues, bfThres, dsimFullAnnoFrame, flybaseToGene)
    print ("[PLOT_REGIONS_DONE]")


if __name__ == "__main__":
    main()

