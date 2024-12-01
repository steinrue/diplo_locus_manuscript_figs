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
DIFF_NE = {
    'g4000' : 10000,
    'g160' : 1000,
}
# selection
FIXED_H = ["0.5", "1.0"]
MAX_ABS_S2 = "0.75"
NUM_GRID = 31
SELECTION_CHANGE_TIMES = {
    'g4000' : [1000,3000],
    'g160' : [40,120],
}
# diplo-locus
DIPLO_LOCUS_EXECUTABLE = 'DiploLocus-likelihood'
NUM_HMM_STATES = 500




def run_diplo_locus():

    # create the directory for the output
    # RESULTS_DIR.mkdir()
    # RESULTS_DIR.mkdir (exist_ok=True)
    assert (RESULTS_DIR.is_dir())

    # build commands in each scenario and run them
    assert (DIFF_NE.keys() == SELECTION_CHANGE_TIMES.keys())
    for totalGensKey in DIFF_NE.keys():
        for stringH in FIXED_H:
            for toChange in ["const", "var"]:

                # name of scenario and files
                thisScenario = f"{totalGensKey}_h{stringH}_{toChange}"
                print (f"+++++ {thisScenario}")

                # input file
                sampleFile = pathlib.Path (RESULTS_DIR, f"{thisScenario}_samples.count")
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

                # DiploLocus-likelihood  -i example_counts.txt --u01 1e-6 --Ne 2500 --gen_time 1 --sample_times 0,10,20
                # --init uniform --fix_h 0.5 --geom_s2_range="-0.05,0.05,20" -o example_out --get_off_grid_max --get_MLR --get_chi2_pval
                diploCmdList = [
                    f"{DIPLO_LOCUS_EXECUTABLE}",
                    f"-i {str(sampleFile)}",
                    f"--u01 {DIFF_M_ALPHA} --u10 {DIFF_M_BETA}",
                    f"--Ne {DIFF_NE[totalGensKey]}" ,
                    "--gen_time 1",
                    f"--sample_times {','.join([str(int(x)) for x in samplingTimes])}",
                    "--init uniform",
                    f"--fix_h {stringH}",
                    f'--geom_s2_range="-{MAX_ABS_S2},{MAX_ABS_S2},{NUM_GRID}"',
                    f"-o {str(pathlib.Path (RESULTS_DIR, f'{thisScenario}_results'))}",
                    f"--numStates {NUM_HMM_STATES}",
                    # "--get_off_grid_max --get_MLR --get_chi2_pval",
                    "--get_off_grid_max --get_MLR",
                ]
                # might have to add piecewise stuff
                if (toChange == "var"):
                    # weird file to output
                    pieceFrame = pandas.DataFrame()
                    thisChangeTimes = SELECTION_CHANGE_TIMES[totalGensKey]
                    pieceFrame['starts'] = numpy.concatenate([[0],thisChangeTimes])
                    pieceFrame['ends'] = numpy.concatenate([thisChangeTimes, [4000]])
                    pieceFrame['s1'] = [0,'x',0]
                    pieceFrame['s2'] = [0,'x',0]
                    pieceFile = pathlib.Path (RESULTS_DIR, f"{thisScenario}_diplo.pieces")
                    pieceFrame.to_csv (pieceFile, sep='\t', index=False, header=False)
                    # and add to command
                    diploCmdList.append (f"--piecewise --specify_piece {str(pieceFile)}")
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