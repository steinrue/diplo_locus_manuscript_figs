import pathlib
import numpy
import time
# custom modules
import execution_engine




# how to execute
THE_EXECUTOR = execution_engine.FileExecutor ("SLURM_ANALYZE.txt", append=True)
# THE_EXECUTOR = execution_engine.SubprocessExecutor()
RSCRIPT_BINARY = "Rscript"
# where is the simulated data?
RESULTS_DIR = pathlib.Path ("results_lls")
# fixed parameters
DIFF_NE = {
    'g4000' : 10000,
    'g160' : 1000,
}
# selection
FIXED_H = ["0.5", "1.0"]




def llsScriptR (inputFile, outputFile, Ne, fixedH):

    # this could be used to get bootstrap p-values: simulate.p.value = FALSE, N.pval = 1000
    # but not for now, because takes time and unstable

    if (numpy.isclose (fixedH, 0.5)):
        dominanceString = ', h = 0.5, method="LLS"'
    else:
        dominanceString = f', h = {fixedH:.4f}, method="NLS"'

    schkriptSchtring = f'''
#  load required libraries
library("data.table")
library("foreach")
library("stringi")
library("matrixStats")
library("Rcpp")
# and load poolseq
library("poolSeq")

inputFile <- "{str(inputFile)}"

# first read just the beginning of the file to get the sampling times
con <- file (inputFile, "r")
timesLine <- NULL
while (TRUE) {{
  line <- readLines(con, n = 1)
  if (length(line) == 0) {{
    break
  }}
  if (startsWith (line, "## Sampling times:")) {{
    timesLine <- line
    break    
  }}
}}
close(con)
stopifnot (!is.null(timesLine))

# should have the line with times, get times out
timeString <- strsplit(strsplit(timesLine, ':')[[1]][2], ";")[[1]][1]
samplingGens <- lapply (strsplit(timeString, ","), as.numeric)[[1]]

# load the rest of the data
allData <- read.table (inputFile, header = TRUE)
# should have 2 numbers per gen and one ID column
stopifnot (ncol(allData) == 2*length(samplingGens)+1)

# put it into convenient format
timeIdx <- seq(1, length(samplingGens), 1)
derCol <- paste("d", timeIdx, sep='')
nCol <- paste("n", timeIdx, sep='')
samples <- allData[derCol]
sampleSizes <- allData[nCol]

# run through replicates and estimate
sEstimate <- c()
for (repIdx in 1:nrow(samples)) {{
  freqs <- samples[repIdx,]/sampleSizes[repIdx,]
  covs <- sampleSizes[repIdx,]
  thisEst <- estimateSH (freqs, samplingGens, {Ne}, cov = covs{dominanceString})
  sEstimate <- c(sEstimate, thisEst$s)
}}

# output the results nicely formatted
results <- data.frame (allData['ID'], sEstimate)
write.table (results, "{str(outputFile)}", row.names=FALSE, sep="\\t", quote=FALSE)

'''
    
    return schkriptSchtring


def run_lls():

    # create the directory for the output
    # RESULTS_DIR.mkdir()
    # RESULTS_DIR.mkdir (exist_ok=True)
    assert (RESULTS_DIR.is_dir())


    # build commands in each scenario and run them
    for totalGensKey in DIFF_NE.keys():
        # run both hs
        for stringH in FIXED_H:
            # but only run analysis for const coeff
            # for toChange in ["const", "var"]:
            for toChange in ["const"]:

                # name of scenario and files
                thisScenario = f"{totalGensKey}_h{stringH}_{toChange}"
                print (f"+++++ {thisScenario}")

                # input file
                sampleFile = pathlib.Path (RESULTS_DIR, f"{thisScenario}_samples.count")
                # check whether it exists
                assert (sampleFile.is_file())

                # output file
                outputFile = pathlib.Path (RESULTS_DIR, f'{thisScenario}_lls_results.table')

                # get the R script together
                theSchkript = llsScriptR (sampleFile, outputFile, DIFF_NE[totalGensKey], float(stringH))
                schkriptName = pathlib.Path (RESULTS_DIR, f'{thisScenario}_lls.R')

                # write the script
                # print (schkriptName)
                ofs = open (schkriptName, 'w')
                ofs.write (theSchkript)
                ofs.close()

                # run the script
                schkriptCmd = f"{RSCRIPT_BINARY} {str(schkriptName)}"
                print ("[LLS]")
                start_time = time.time()
                THE_EXECUTOR.runCmd (schkriptCmd)
                end_time = time.time()
                print (end_time - start_time)
                print ("[LLS_END]")


def main():

    run_lls()


if __name__ == "__main__":
    main()