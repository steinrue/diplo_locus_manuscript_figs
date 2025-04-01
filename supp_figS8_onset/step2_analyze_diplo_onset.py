import pathlib
import numpy
import pandas
import time
# custom modules
import execution_engine




# how to execute
THE_EXECUTOR = execution_engine.FileExecutor ("SLURM_ANALYZE.txt", append=False)
# THE_EXECUTOR = execution_engine.SubprocessExecutor()
# where is the simulated data?
RESULTS_DIR = pathlib.Path ("results_diplo_locus")
# fixed parameters
# diffusion
DIFF_M_ALPHA = 0
DIFF_M_BETA = DIFF_M_ALPHA
DIFF_NE = 10000
# selection
FIXED_H = "0.5"
MAX_ABS_S2 = "0.1"
NUM_S2_GRID = 71
# FIRST_EPOCH_S2_FILE = pathlib.Path ("first_epoch_s2.txt")
# onset
MIN_ONSET_GEN = 0
MAX_ONSET_GEN = 4000
NUM_ONSET_GRID = 51
# diplo-locus
DIPLO_LOCUS_EXECUTABLE = 'DiploLocus-likelihood'
NUM_HMM_STATES = 500




def run_diplo_locus():

    # create the directory for the output
    # RESULTS_DIR.mkdir()
    # RESULTS_DIR.mkdir (exist_ok=True)
    assert (RESULTS_DIR.is_dir())

    # get the s2s for the first epoch out of the given file
    # firstEpochS2 = numpy.loadtxt (FIRST_EPOCH_S2_FILE)
    firstEpochS2 = [0]

    # one command for each possible onset time
    # +/- 1 just for good measure
    # for thisOnsetTime in numpy.linspace (MIN_ONSET_GEN+1, MAX_ONSET_GEN-1, NUM_ONSET_GRID):
    onsetTimes = numpy.linspace (MIN_ONSET_GEN+1, MAX_ONSET_GEN-1, NUM_ONSET_GRID)
    onsetTimes = numpy.concatenate ([[float("nan")], onsetTimes])
    for thisOnsetTime in onsetTimes:

        if (numpy.isnan(thisOnsetTime)):
            # thisOnsetTime should never be used again
            firstNan = True

        # input file
        sampleFile = pathlib.Path (RESULTS_DIR, f"onset_samples.count")
        # check whether it exists
        assert (sampleFile.is_file())

        # read a bit to get the sampling times
        samplingTimes = None
        ifs = open (sampleFile, 'r')
        for line in ifs.readlines ():
            line = line.strip()
            if (line.startswith ('## Sampling times:')):
                samplingTimes = numpy.array ([float(x) for x in line.split(':')[-1].split(';')[0].split(',')])
                break
        ifs.close()
        assert (samplingTimes is not None)
        # print (samplingTimes)

        # do this for every s2 in the first epoch
        for thisFirstS2 in firstEpochS2:

            if (numpy.isnan(thisOnsetTime)):
                # thisOnsetTime should never be used again
                thisScenario = f"onset_nan"
            else:
                # should be ok to make it an int
                thisOnsetTime = int(thisOnsetTime)
                thisScenario = f"onset_{thisOnsetTime}_{thisFirstS2:.6e}"

            # name of scenario and files
            print (f"+++++ {thisScenario}")

            # DiploLocus-likelihood  -i example_counts.txt --u01 1e-6 --Ne 2500 --gen_time 1 --sample_times 0,10,20
            # --init uniform --fix_h 0.5 --geom_s2_range="-0.05,0.05,20" -o example_out --get_off_grid_max --get_MLR --get_chi2_pval
            diploCmdList = [
                f"{DIPLO_LOCUS_EXECUTABLE}",
                f"-i {str(sampleFile)}",
                f"--u01 {DIFF_M_ALPHA} --u10 {DIFF_M_BETA}",
                f"--Ne {DIFF_NE}" ,
                "--gen_time 1",
                f"--sample_times {','.join([str(int(x)) for x in samplingTimes])}",
                "--init uniform",
                f"--fix_h {FIXED_H}",
                f'--geom_s2_range="-{MAX_ABS_S2},{MAX_ABS_S2},{NUM_S2_GRID}"',
                f"-o {str(pathlib.Path (RESULTS_DIR, f'{thisScenario}_results'))}",
                f"--numStates {NUM_HMM_STATES}",
                # "--get_off_grid_max --get_MLR --get_chi2_pval",
                # "--get_off_grid_max --get_MLR",
            ]

            # the one withtout onset should just work, right?
            if (not numpy.isnan(thisOnsetTime)):
                # add piecewise stuff for onset
                pieceFrame = pandas.DataFrame()
                pieceFrame['starts'] = numpy.concatenate([[0],[thisOnsetTime]])
                pieceFrame['ends'] = numpy.concatenate([[thisOnsetTime], [4000]])
                # pieceFrame['s1'] = [0,'x']
                # pieceFrame['s2'] = [0,'x']
                pieceFrame['s1'] = [thisFirstS2,'x']
                pieceFrame['s2'] = [thisFirstS2,'x']
                pieceFile = pathlib.Path (RESULTS_DIR, f"{thisScenario}_diplo.pieces")
                pieceFrame.to_csv (pieceFile, sep='\t', index=False, header=False)

                # and add to command
                diploCmdList.append (f"--piecewise --specify_piece {str(pieceFile)}")
            else:
                if (firstNan):
                    firstNan = False
                else:
                    break

            print ("[DIPLO_LOCUS]")
            diploCmd = ' '.join(diploCmdList)

            # and run it, but take the time
            start_time = time.time()
            THE_EXECUTOR.runCmd (diploCmd)
            end_time = time.time()

            print (end_time - start_time)
            print ("[DIPLO_LOCUS_END]")


def main():

    run_diplo_locus()


if __name__ == "__main__":
    main()