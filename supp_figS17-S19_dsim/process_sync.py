import pandas
import numpy
import sys
import pickle
import bz2
import time




# accoding to https://code.google.com/archive/p/popoolation2/wikis/Tutorial.wiki
# The allele frequencies are in the format A:T:C:G:N:del, i.e: count of bases 'A', count of bases 'T',... and deletion count in the end (character '*' in the mpileup)
# we do the conversion to integers later
NUC_IDX2CHAR = ['A', 'T', 'C', 'G', 'N', 'DEL']
NUC_CHAR2IDX = {}
for (idx, char) in enumerate(NUC_IDX2CHAR):
    NUC_CHAR2IDX[char] = idx
NUC_RANGE = numpy.arange(NUC_CHAR2IDX['A'],NUC_CHAR2IDX['G'] + 1)
# column indices in sync-file
SYNC_POSITION_COL = 1
SYNC_REFERENCE_COL = 2
SYNC_OFFSET = 3

# number of replicates and timepoints in this specific sync-file
SYNC_NUM_REPLICATES = 10
SYNC_NUM_TIMEPOINTS = 7
SYNC_CMH_PVALUES_COL = SYNC_OFFSET+SYNC_NUM_TIMEPOINTS*SYNC_NUM_REPLICATES




def extractGeno (rawSyncEntry):
    # as per sync-file standard, split stuff by colon
    return rawSyncEntry.strip().split(':')


def loadInputFile (inputFile):

    # read the file
    syncFrame = pandas.read_csv (inputFile, sep='\t', header=None)

    # get some info about the SNPs
    thePositions = syncFrame.iloc[:,SYNC_POSITION_COL].to_numpy()
    refAlleles = syncFrame.iloc[:,SYNC_REFERENCE_COL].to_numpy()
    refAlleleIdxs = numpy.vectorize(lambda x : NUC_CHAR2IDX[x])(refAlleles)
    # if this holds everywhere, that's pretty good
    assert (refAlleleIdxs.min() >= NUC_RANGE.min())
    # I guess we might have some 'N' in the reference
    assert (refAlleleIdxs.max() <= NUC_RANGE.max() + 1), (refAlleleIdxs.max(), NUC_RANGE.max() + 1)

    # transform the data into a formate we like
    # first get a list of numbers for every locus x (rep x time)
    preReadCounts = syncFrame.iloc[:,SYNC_OFFSET:SYNC_OFFSET + SYNC_NUM_TIMEPOINTS*SYNC_NUM_REPLICATES].map(extractGeno).to_numpy()
    # then transform to int and shape into full proper numpy array
    originalShape = preReadCounts.shape
    readCounts = numpy.array(preReadCounts.flatten().tolist(), dtype=int).reshape(originalShape + (-1,))

    # get the p-values just for good measure
    cmhPValues = syncFrame.iloc[:,SYNC_CMH_PVALUES_COL].to_numpy()

    # return some results
    return (thePositions, refAlleleIdxs, readCounts, cmhPValues)


def convertData (thePositions, refAlleleIdxs, readCounts, cmhPValues):

    # get major/minor alleles: sort, but decreasing
    # also, should only be done for the counts of real alleles
    sortedIdxs = numpy.flip (numpy.argsort (readCounts[:,:,NUC_RANGE], axis=2), axis=2)
    # these are undefined if there is a tie
    # we could assert something, but it does not even matter for allele freq
    majorAlleleIdx = sortedIdxs[:,:,0]
    minorAlleleIdx = sortedIdxs[:,:,1]
    # polarize by major allele at first time point
    firstMajorAlleleIdx = majorAlleleIdx[:,:SYNC_NUM_REPLICATES]
    firstMinorAlleleIdx = minorAlleleIdx[:,:SYNC_NUM_REPLICATES]


    # # in this data, some sequencing pools have:
    # print ("[GET_STATS]")
    # # - reads for DEL
    # hasDeletions = readCounts[:,:,NUC_CHAR2IDX['DEL']] > 0
    # print (f"num pools with deletions: {numpy.sum(hasDeletions)}")
    # # - reads for N (seems to be very rare)
    # hasMissing = readCounts[:,:,NUC_CHAR2IDX['N']] > 0
    # print (f"num pools with missing: {numpy.sum(hasMissing)}")
    # # - N in the reference sequence
    # refIsMissing = (refAlleleIdxs == (NUC_RANGE.max() + 1))
    # print (f'num missing in "base" seq: {numpy.sum(refIsMissing)}')
    # # - reads for three alleles (most only a few, but sometimes more)
    # numProperAlleles = numpy.sum(readCounts[:,:,NUC_RANGE] > 0, axis=2)
    # hasNoProperAlleles = (numProperAlleles < 1)
    # isMonoAllelic = (numProperAlleles == 1)
    # isBiAllelic = (numProperAlleles == 1) | (numProperAlleles == 2)
    # moreThanBiAllelic = ~(isBiAllelic | hasNoProperAlleles)
    # print (f'num pools more than 2 allels: {numpy.sum(moreThanBiAllelic)}')
    # print (readCounts[moreThanBiAllelic])
    # # - reads for major and minor allele that are not the reference allele (often only few reads, but sometimes more)
    # trueBiAllelic = isBiAllelic & ~isMonoAllelic
    # majorIsRef = (majorAlleleIdx == refAlleleIdxs[:, numpy.newaxis])
    # minorIsRef = (minorAlleleIdx == refAlleleIdxs[:, numpy.newaxis])
    # refIsNotMajorOrMinor = ~(majorIsRef | minorIsRef)
    # wrongAlleleCount = (trueBiAllelic & refIsNotMajorOrMinor & ~refIsMissing[:, numpy.newaxis])
    # print (f'# pools major/minor not ref: {numpy.sum(wrongAlleleCount)}')
    # print (readCounts[wrongAlleleCount])
    # # THIS SHOULD ALL BE FILTERED PROPERLY
    # print ("[GET_STATS_DONE]")
    # # allowing for just one or two reads to be wrong might already save a lot
    # # but for now, we just take the numbers of reads for the major and minor allele


    # get the major and minor read counts
    polarMajorAlleleIdx = numpy.tile(firstMajorAlleleIdx, (1, SYNC_NUM_TIMEPOINTS))
    polarMinorAlleleIdx = numpy.tile(firstMinorAlleleIdx, (1, SYNC_NUM_TIMEPOINTS))
    (majorRows, majorCols) = numpy.indices (polarMajorAlleleIdx.shape)
    majorReadCounts = readCounts[majorRows, majorCols, polarMajorAlleleIdx]
    minorReadCounts = readCounts[majorRows, majorCols, polarMinorAlleleIdx]

    # return stuff
    return (thePositions, majorReadCounts, minorReadCounts, cmhPValues)


def saveData (outputFile, thePositions, majorReadCounts, minorReadCounts, cmhPValues):

    # pickle it
    pickleData = {
        'thePositions' : thePositions,
        'majorReadCounts' : majorReadCounts,
        'minorReadCounts' : minorReadCounts,
        'cmhPValues' : cmhPValues,
    }
    pickle.dump (pickleData, bz2.open (outputFile, "wb"))


def processInputFile (inputFile, outputFile):

    # load the input
    startTime = time.time()
    print ("[LOAD]")
    (thePositions, refAlleleIdxs, readCounts, cmhPValues) = loadInputFile (inputFile)
    print (readCounts.shape, len(thePositions), len(refAlleleIdxs))
    print (f"[LOAD_DONE] {int(time.time() - startTime)}")
    

    # get the counts
    startTime = time.time()
    print ('[CONVERT]')
    (thePositions, majorReadCounts, minorReadCounts, cmhPValues) = convertData (thePositions, refAlleleIdxs, readCounts, cmhPValues)
    print (f"[CONVERT_DONE] {int(time.time() - startTime)}")


    # and save them
    startTime = time.time()
    print ('[SAVE]')
    saveData (outputFile, thePositions, majorReadCounts, minorReadCounts, cmhPValues)
    print (f"[SAVE_DONE] {int(time.time() - startTime)}")


def main():

    if (len(sys.argv) != 3):
        print ("usage: python <script_name> <input_file.sync.gz> <output_file.pickle.bz2>")
        sys.exit (1)

    inputFile = sys.argv[1]
    outputFile = sys.argv[2]
    print (f"[EXTRACT] {inputFile} -> {outputFile}")
    
    processInputFile (inputFile, outputFile)


if __name__ == "__main__":
    main()

