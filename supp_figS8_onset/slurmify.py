import sys
import os
import pathlib




# TMP_DIR = 'slurmify_tmp'
NUM_CORES = 1
# NUM_CORES = 32
NUM_GB = 16
# NUM_GB = 32
# NUM_GB = 64
NUM_HOURS = 4




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


def prepareSlurmScripts (commandList, metaSubmitScript, basename):

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

        # write the actual command
        ofs.write (thisCmd + '\n')

        # be nice
        ofs.close()

        # add something to a meta script
        metafs.write ("sbatch %s\n" % thisScriptName)
    
        cmdIdx += 1
        
    # close and then we done?
    metafs.close()


def main():

    # "parse" coomand line arguments
    if (len(sys.argv) != 3):
        print ("usage: python <script_name> <CMD-file> <base_name>")
        return
    
    cmdFile = sys.argv[1]
    baseName = sys.argv[2]
    # print (cmdFile, baseName)

    # get the commands
    commandList = readCommands (cmdFile)

    prepareSlurmScripts (commandList, f'submit_{baseName}.sh', baseName)


if __name__ == "__main__":
    main()
