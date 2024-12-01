import pathlib
import numpy
import time
# custom modules
import execution_engine

metaRNG = numpy.random.default_rng (4715)




# how to execute
THE_EXECUTOR = execution_engine.FileExecutor ("SLURM_ANALYZE.txt", append=True)
# THE_EXECUTOR = execution_engine.SubprocessExecutor()
# WFABC1_BINARY = "./WFABC_v1.1/binaries/Mac/wfabc_1"
# WFABC2_BINARY = "./WFABC_v1.1/binaries/Mac/wfabc_2"
WFABC1_BINARY = "./WFABC_v1.1/binaries/Linux/wfabc_1"
WFABC2_BINARY = "./WFABC_v1.1/binaries/Linux/wfabc_2"
# where is the simulated data?
RESULTS_DIR = pathlib.Path ("results_wfabc")
# fixed parameters
DIFF_NE = {
    'g4000' : 10000,
    'g160' : 1000,
}
# selection
FIXED_H = ["0.5", "1.0"]




def run_wfabc():

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
                sampleFile = pathlib.Path (RESULTS_DIR, f"{thisScenario}_loci.txt")
                # check whether it exists
                assert (sampleFile.is_file())
                # no explicit output needed

                # first wfabc1 for summary stats
                wfabc1Cmd = f"{WFABC1_BINARY} -nthreads 1 -nboots 0 -seed {metaRNG.integers(99999999)} {sampleFile}"

                # then wfabc2 for actual ABC stuff
                cmdList = [
                    f"{WFABC2_BINARY}",
                    "-nthreads 1",
                    f"-fixed_N {2*DIFF_NE[totalGensKey]}",
                    f"-seed {metaRNG.integers(99999999)}",
                    # "-min_s -0.75 -max_s 0.75",
                ]
                if (totalGensKey == 'g160'):
                    # cmdList += ["-mean_s 0 -sd_s 0.15"]
                    cmdList += ["-mean_s 0 -sd_s 0.075"]
                elif (totalGensKey == 'g4000'):
                    # cmdList += ["-mean_s 0 -sd_s 0.15"]
                    cmdList += ["-mean_s 0 -sd_s 0.0075"]
                else:
                    assert (False)
                floatH = float(stringH)
                if (not numpy.isclose(floatH, 0.5)):
                    cmdList += [f"-min_h {floatH - 0.0001:.5f} -max_h {floatH:.5f}"]
                cmdList += [
                    f"{str(sampleFile)}",
                ]
                wfabc2Cmd = ' '.join(cmdList)

                # combine both commands, so they can be run with one call (hopefully)
                combinedCmd = wfabc1Cmd + '; ' + wfabc2Cmd
                # print (combinedCmd)


                # run things, but togethor
                print ("[WFABC]")
                start_time = time.time()
                THE_EXECUTOR.runCmd (combinedCmd)
                end_time = time.time()
                print (end_time - start_time)
                print ("[WFABC_END]")


def main():

    run_wfabc()


if __name__ == "__main__":
    main()