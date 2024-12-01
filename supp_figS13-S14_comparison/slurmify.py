import sys
import os
import pathlib




# TMP_DIR = 'slurmify_tmp'
NUM_CORES = 1
NUM_GB = 16
NUM_HOURS = 24
# NUM_HOURS = 8

# this command is run to benchmark against
BENCHMARK_CMD = 'DiploLocus-likelihood --u01 0 --u10 0 --Ne 10000 --gen_time 1' \
    + ' -i results_diplo_locus/g4000_h1.0_var_samples.count --sample_times 0,500,1000,1500,2000,2500,3000,3500,4000' \
    + ' --init uniform --fix_h 1.0 --geom_s2_range="-0.75,0.75,31" -o results_diplo_locus/BENCH%s_g4000_h1.0_var_results' \
    + ' --numStates 500 --get_off_grid_max --get_MLR --piecewise' \
    + ' --specify_piece results_diplo_locus/g4000_h1.0_var_diplo.pieces'


def readCommands (commandFile):

    # read file
    ifs = open (commandFile, 'r')
    theLines = ifs.readlines()
    ifs.close()

    # go through and clean
    cleaned = []
    for thisLine in theLines:
        cleaned.append (thisLine.strip())

    return cleaned


def prepareSlurmScripts (commandList, metaSubmitScript, basename, benchmark):

    # set up the tmp dir
    tmpDirName = f"slurmify_{basename}"
    tmpDir = pathlib.Path (tmpDirName)
    # tmpDir.mkdir()
    tmpDir.mkdir (exist_ok=True)


    # set up meta script
    metafs = open (metaSubmitScript, 'w')

    # go through commands
    cmdIdx = 0
    for thisCmd in commandList:

        # no empty commands =)
        if (thisCmd.strip() == ''):
            continue

        # generic name
        thisName = f"slurm_{basename}_{cmdIdx}"
        thisScriptName = pathlib.Path (tmpDir, f"{thisName}.sbatch")
        ofs = open (thisScriptName, "w")

        # needs to be in script
        ofs.write ("#!/bin/bash\n")
        # give it a name
        ofs.write (f"#SBATCH --job-name={thisName}\n")
        # only 1 node
        ofs.write ("#SBATCH --nodes=1\n")
        # only 1 task
        ofs.write ("#SBATCH --ntasks=1\n")
        # we want this many cpus
        ofs.write (f"#SBATCH --cpus-per-task={NUM_CORES}\n")
        # also a bit of time
        ofs.write (f"#SBATCH --time={NUM_HOURS}:00:00\n")
        # some memory
        ofs.write (f"#SBATCH --mem={NUM_GB}gb\n")
        # output
        ofs.write (f"#SBATCH --output={tmpDirName}/%x_j%j.out\n")
        ofs.write (f"#SBATCH --error={tmpDirName}/%x_j%j.err\n")
        ofs.write ("\n")


        # unfortunately, python needs this
        ofs.write ('export PYTHONUNBUFFERED="everything_is_fine"\n')


        # maybe benchmark shenanigans
        if (benchmark):
            ofs.write ('echo "[PRE_BENCH]"\n')
            ofs.write ('date +%s\n')
            ofs.write ('echo "[PRE_BENCH]"\n')

            ofs.write (BENCHMARK_CMD % (f"{basename}{cmdIdx}") + '\n')
            ofs.write ('echo "[PRE_CMD]"\n')
            ofs.write ('date +%s\n')
            ofs.write ('echo "[PRE_CMD]"\n')

        # write the actual command
        ofs.write (thisCmd + '\n')

        # final benchmark shenanigans
        if (benchmark):
            ofs.write ('echo "[POST_CMD]"\n')
            ofs.write ('date +%s\n')
            ofs.write ('echo "[POST_CMD]"\n')


        # be nice
        ofs.close()

        # add something to a meta script
        metafs.write ("sbatch %s\n" % thisScriptName)
    
        cmdIdx += 1
        
    # close and then we done?
    metafs.close()


def main():

    # "parse" coomand line arguments
    if (len(sys.argv) != 4):
        print ("usage: python <script_name> <CMD-file> <base_name> <benchmark or no_benchmark>")
        return
    
    cmdFile = sys.argv[1]
    baseName = sys.argv[2]
    benchmarkString = sys.argv[3]
    # print (cmdFile, baseName)

    # to benchmark or not to benchmark
    if (benchmarkString.strip() == 'benchmark'):
        benchmark = True
    elif (benchmarkString.strip() == 'no_benchmark'):
        benchmark = False
    else:
        assert (False), 'Provide either "benchmark" or "no_benchmark" as option.'

    # get the commands
    commandList = readCommands (cmdFile)

    prepareSlurmScripts (commandList, f'submit_{baseName}.sh', baseName, benchmark)


if __name__ == "__main__":
    main()
