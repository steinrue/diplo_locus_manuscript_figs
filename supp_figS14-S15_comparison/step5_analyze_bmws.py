import pathlib
import numpy
import time
# custom modules
import execution_engine




# how to execute
THE_EXECUTOR = execution_engine.FileExecutor ("SLURM_ANALYZE.txt", append=True)
# THE_EXECUTOR = execution_engine.SubprocessExecutor()
BMWS_BINARY = {
    'const' : 'bmws',
    'var' : 'python bmws/bmws_change.py',
}
# where is the simulated data?
RESULTS_DIR = pathlib.Path ("results_bmws")
# fixed parameters
DIFF_NE = {
    'g4000' : 10000,
    'g160' : 1000,
}
# selection
SELECTION_CHANGE_TIMES = {
    'g4000' : [1000,3000],
    'g160' : [40,120],
}
FIXED_H = ["0.5"]




def run_bmws():

    # create the directory for the output
    # RESULTS_DIR.mkdir()
    # RESULTS_DIR.mkdir (exist_ok=True)
    assert (RESULTS_DIR.is_dir())


    # build commands in each scenario and run them
    assert (DIFF_NE.keys() == SELECTION_CHANGE_TIMES.keys())
    for totalGensKey in DIFF_NE.keys():
        # only additive
        for stringH in FIXED_H:
            # but const and varying
            for toChange in ["const", "var"]:

                # name of scenario and files
                thisScenario = f"{totalGensKey}_h{stringH}_{toChange}"
                print (f"+++++ {thisScenario}")

                # input files
                vcfFile = pathlib.Path (RESULTS_DIR, f"{thisScenario}.vcf")
                assert (vcfFile.is_file())
                metaFile = pathlib.Path (RESULTS_DIR, f"{thisScenario}.meta")
                assert (metaFile.is_file())

                # output file
                outFile = pathlib.Path (RESULTS_DIR, f"{thisScenario}.results")

                # command for bmws
                cmdList = [
                    f"{BMWS_BINARY[toChange]}",
                    "analyze",
                    f"{str(vcfFile)}",
                    f"{str(metaFile)}",
                    "-d pseudohaploid",
                    "-g 1",
                    f"-n {DIFF_NE[totalGensKey]}",
                ]
                # figure out what to do with regularization
                theLambda = -1
                if (toChange == 'const'):
                    theLambda = 7
                elif (toChange == 'var'):
                    theLambda = 4.5
                assert (theLambda > 0)
                cmdList += [f"-l {theLambda:.2f}"]
                # need to add this for var
                if (toChange == 'var'):
                    # selection interval
                    thisSelInterval = SELECTION_CHANGE_TIMES[totalGensKey]
                    cmdList += [f"--selection_interval {int(thisSelInterval[0])},{int(thisSelInterval[1])}"]
                # and regular stuff
                cmdList += [
                    f">{str(outFile)}",
                ]
                bmwsCmd = ' '.join(cmdList)

                # run things
                print ("[BMWS]")
                start_time = time.time()
                THE_EXECUTOR.runCmd (bmwsCmd)
                end_time = time.time()
                print (end_time - start_time)
                print ("[BMWS_END]")


def main():

    run_bmws()


if __name__ == "__main__":
    main()