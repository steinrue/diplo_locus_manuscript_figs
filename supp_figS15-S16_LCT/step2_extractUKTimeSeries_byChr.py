import os
import sys

import pandas
import pathlib
import collections
import subprocess

ALLEN_VERSION = "v54.1"
ALLEN_ROOT = pathlib.Path(f'aadr_{ALLEN_VERSION}_1240K/aadr_{ALLEN_VERSION}_1240K_public')
INDIVIDUALS_FILE = pathlib.Path('extracted', f'UK_{ALLEN_VERSION}_1240K_noFam_strictPASS_from4500.table')
# assume that directory "extracted" exists (if files exist, they are overwritten)
OUTPUT_VCF_DIR = pathlib.Path('extracted')

# some scripts that we need
EIGENSTRAT_CONVERSION_SCRIPT = pathlib.Path('gdc/eigenstrat2vcf.py')
CORRECT_HAPLO_ENCODING_SCRIPT = pathlib.Path('diploToHaploVcf.py')
# field different in different versions
# v54.1 & v52.2
GENETIC_ID = 'Genetic_ID'


def extractTimeSeries(c):
    # new directory
    # os.makedirs (vcfDir)

    # load the individuals and make an individuals file for conversion
    thisBasename = os.path.splitext(os.path.basename(INDIVIDUALS_FILE))[0]
    eigentstratIndFile = pathlib.Path(OUTPUT_VCF_DIR, f'{thisBasename}.inds')
    individualsFrame = pandas.read_csv(INDIVIDUALS_FILE, sep='\t')
    print(individualsFrame.shape)

    ofs = open(eigentstratIndFile, 'w')
    for thisInd in individualsFrame[GENETIC_ID]:
        ofs.write(f"{thisInd}\n")
    ofs.close()

    # load all the snps
    snpFile = pathlib.Path(f'{ALLEN_ROOT}.snp')
    snpFrame = pandas.read_csv(snpFile, delim_whitespace=True)
    print(snpFrame.shape)

    chromHist = collections.Counter(snpFrame.iloc[:, 1])
    # make sure the requested chromosome is included
    assert (c in chromHist.keys()), (c, chromHist.keys())

    # parse select chromosome
    chromName = str(c)
    # 23 is X
    if c == 23:
        chromName = 'X'
    # 24 is Y
    elif c == 24:
        chromName = 'Y'
    else:
        pass

    # make a snpfile for eigenstrat
    eigenstratSnpFile = pathlib.Path(OUTPUT_VCF_DIR, f'{thisBasename}_c{chromName}.snps')
    ofs = open(eigenstratSnpFile, 'w')
    thisSnpFrame = snpFrame.loc[snpFrame.iloc[:, 1] == c]
    for thisSnp in thisSnpFrame.iloc[:, 0]:
        ofs.write(f"{thisSnp}\n")
    ofs.close()

    # prepare the output file
    outputDiploVCF = os.path.join(OUTPUT_VCF_DIR, f'{thisBasename}_c{chromName}.diplo_vcf')

    # and then put eigentstrat command together
    stratCmd = f"python {EIGENSTRAT_CONVERSION_SCRIPT} -r {ALLEN_ROOT} -i {eigentstratIndFile} -s {eigenstratSnpFile} >{outputDiploVCF}"
    print(stratCmd)

    # run it
    result = subprocess.run(stratCmd, shell=True, capture_output=True)
    print(result.stdout)

    # clean up
    print(f"Removing: {eigenstratSnpFile}")
    os.remove(eigenstratSnpFile)

    # make it a proper vcf with haploid calls encoded correctly
    outputVCF = os.path.join(OUTPUT_VCF_DIR, f'{thisBasename}_c{chromName}.vcf')

    # put command together to convert it to pseudo-haploid
    haploCmd = f"python {CORRECT_HAPLO_ENCODING_SCRIPT} {outputDiploVCF} > {outputVCF}"
    print(haploCmd)

    # run it
    result = subprocess.run(haploCmd, shell=True, capture_output=True)
    print(result.stdout)

    # clean up
    print(f'Removing: {outputDiploVCF}')
    os.remove(outputDiploVCF)


def main():
    chrom = sys.argv[1]
    extractTimeSeries(int(chrom))


if __name__ == "__main__":
    main()
