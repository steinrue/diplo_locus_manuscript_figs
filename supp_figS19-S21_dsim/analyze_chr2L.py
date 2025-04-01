import execution_engine
import numpy
import pathlib
import os

# what are the input files?
# CHR_NAMES = [f"chr{x}" for x in ['2L', '2R', '3L', '3R', '4', 'X']]
CHR_NAMES = [f"chr{x}" for x in ['2L']]
INPUT_FILES = {}
INPUT_DIR = 'data/'
for thisChr in CHR_NAMES:
    INPUT_FILES[thisChr] = pathlib.Path (INPUT_DIR, f"{thisChr}_F0-F60_all_reps.pickle.bz2")
OUTPUT_DIR = 'results/'
OUTPUT_FILE_SKELETON = str (pathlib.Path (OUTPUT_DIR, "analyzed_%s_F0-F60_rep%d.pickle.bz2"))
# how many replicates in each file
NUM_REPLICATES = 10
# misc
THE_EXECUTOR = execution_engine.FileExecutor ('SLURM_ANALYZE_ALL.txt', append=False)
# THE_EXECUTOR = execution_engine.SubprocessExecutor()
NUM_CPUS = 32
BLOCK_SIZE = 100000
ANALYZE_REPLICATE_SCRIPT = 'analyze_replicate_blocked.py'




def analyzeAllData():

    # create directory for results
    # os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(OUTPUT_DIR)
    
    # go through all files and replicates and put the respective commands together
    for (thisChr, thisInputFile) in INPUT_FILES.items():
        for thisRep in numpy.arange(NUM_REPLICATES):

            # put command together for this one
            #  python analyze_replicate_blocked.py 8 1000 chr4_F0-F60_all_reps.pickle.bz2 0 analyzed_chr4_F0-F60_rep0.pickle.bz2
            thisCmdPieces = [
                f"python {ANALYZE_REPLICATE_SCRIPT} {NUM_CPUS} {BLOCK_SIZE}",
                f"{thisInputFile} {thisRep}",
                f"{OUTPUT_FILE_SKELETON % (thisChr, thisRep)}"
            ]
            thisCmd = " ".join(thisCmdPieces)
            print (thisCmd)
            THE_EXECUTOR.runCmd (thisCmd)


def main():

    # do it
    analyzeAllData()


if __name__ == "__main__":
    main()

