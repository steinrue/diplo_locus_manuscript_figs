import pandas
import sys
import collections
import numpy




def diploToHaploVcf (inputVcf):
    
    # we need the comments at the beginning of the vcf file though
    ifs = open (inputVcf, "r")
    for line in ifs:
        # this trick allows us to get the vcf header into the dataframe
        if (line.startswith ("##")):
            print (line.strip())
        else:
            assert (line.startswith ('#CHROM'))
            # get the header
            columnNames = line.strip().split()
            break

    # try to do the rest with pandas
    vcfFrame = pandas.read_csv (ifs, sep='\t', header=None, names=columnNames)

    # close the stream here
    ifs.close()

    # which columns have the genotypes?
    # genotypes start after format column
    assert ('FORMAT' in columnNames)
    firstGenoColumn = columnNames.index('FORMAT') + 1
    # until the end

    # iterate over genotype columns
    nrInds = 0
    for cIdx in numpy.arange (firstGenoColumn, vcfFrame.shape[1]):

        # get the name of the individual
        indName = columnNames[cIdx]
        nrInds +=1
        
        # count the genotypes for this individuals
        thisGeno = vcfFrame.iloc[:,cIdx]
        thisCounter = collections.Counter (thisGeno)
        thisGenoSet = set(thisCounter.keys())
        
        # we only want valid genotypes (also, if missing, is missing in both)
        assert (thisGenoSet.issubset (set(['./.', '0/0', '0/1', '1/0', '1/1'])))

        # do we have any heterozygotes?
        if (thisGenoSet.issubset (set(['./.', '0/0', '1/1']))):
            # no, only homozygotes
            # it is very likely that this is a pseudo-haploid individual
            # so we make this assumption
            assert (("DG" not in indName) or ("UDG" in indName)), indName

            # change it to be a proper haplotype
            vcfFrame.iloc[:,cIdx] = thisGeno.replace (to_replace=['./.', '0/0', '1/1'], value=['.', '0', '1'])

        else:
            # we have heterozygotes, so this should be a regular diplotype
            assert (("DG" in indName) and ("UDG" not in indName)), indName
            # leave it untouched
            pass

    # write it to stdout
    # header is alread in stdout, so just write the rest
    # seems like WINDOWS needs us to specify a lineterminator
    vcfFrame.to_csv (sys.stdout, sep='\t', index=False, lineterminator='\n')


def main():
    
    # need one parameter
    if (len(sys.argv) != 2):
        print ("usage: python <script_name> <input_vcf_file>")
        exit (-1)

    inputVcfFile = sys.argv[1]

    diploToHaploVcf (inputVcfFile)


if __name__ == "__main__":
    main()