import pandas
import numpy
import collections
import pathlib

# these have to be changed to modify input and output
ALLEN_VERSION = "v54.1"
ANNO_FILE = pathlib.Path(f"aadr_{ALLEN_VERSION}_1240K/aadr_{ALLEN_VERSION}_1240K_public.anno")
# assume that directory "extracted" exists (if files exist, they are overwritten)
OUTPUT_FILE = pathlib.Path('extracted', f'UK_{ALLEN_VERSION}_1240K_noFam_strictPASS_from4500.table')

# some constants for the field names
# might have to be adjusted based on different versions (and do assume compressed column names)
# v54
POLITICAL_ENTITY = 'Political_Entity'
GENETIC_ID = 'Genetic_ID'
TARGET_COVERAGE = '1240k_coverage'
SNPS_HIT = 'SNPs_hit_1240k'
# NO_RELATIVE_ENTRIES = ['n/a (no relatives detected)', 'n/a(norelativesdetected)', 'n/a (No relatives detected)']


def loadAnnotationFile(annoFile):
    # read it
    anno = pandas.read_csv(annoFile, sep='\t', dtype=str)

    # reduce the column names
    # fuse first two names, or take the only one
    newColumns = []
    for c in anno.columns:
        fields = c.split()
        if len(fields) <= 1:
            # just this one
            newColumns.append(fields[0])
        elif (fields[0] == "SNPs") and (fields[1] == "hit"):
            # this one a bit special
            if any([x == '1240k' for x in fields]):
                newColumns.append(f"{fields[0]}_{fields[1]}_1240k")
            elif any([x == 'HO' for x in fields]):
                newColumns.append(f"{fields[0]}_{fields[1]}_HO")
            elif any([x == '3.2M' for x in fields]):
                newColumns.append(f"{fields[0]}_{fields[1]}_3.2M")
            elif any([x == 'non-padded' for x in fields]):
                newColumns.append(f"{fields[0]}_{fields[1]}_non-padded")
            else:
                # this should be ok for version 50.0?
                newColumns.append(f"{fields[0]}_{fields[1]}")
        elif (fields[0] == "Y") and (fields[1] == "haplogroup"):
            # this one's special too
            if any([x == 'terminal' for x in fields]):
                newColumns.append(f"{fields[0]}_{fields[1]}_terminal")
            elif any([x == 'ISOGG' for x in fields]):
                if any([x == 'curation' for x in fields]):
                    newColumns.append(f"{fields[0]}_{fields[1]}_ISOGG_curation")
                else:
                    newColumns.append(f"{fields[0]}_{fields[1]}_ISOGG")
            else:
                assert False
        elif (fields[0] == "Xcontam") and (fields[1] == "ANGSD"):
            # many special snowflakes
            if any([x == 'SNPs' for x in fields]):
                newColumns.append(f"{fields[0]}_{fields[1]}SNPs")
            elif any([x == 'MOM' for x in fields]):
                if any([x == 'point' for x in fields]):
                    newColumns.append(f"{fields[0]}_{fields[1]}_point")
                elif any([x == 'Z-score' for x in fields]):
                    newColumns.append(f"{fields[0]}_{fields[1]}_Z-score")
                elif any([x == 'CI' for x in fields]):
                    newColumns.append(f"{fields[0]}_{fields[1]}_CI")
                else:
                    assert False
            else:
                assert False
        else:
            # fuse first two
            newColumns.append(f"{fields[0]}_{fields[1]}")

    # set the new names
    anno.columns = newColumns
    # print(anno.columns)

    # make sure we don't have duplicates in the new names
    # print(anno.columns)
    assert (numpy.all([x <= 1 for x in collections.Counter(anno.columns).values()]))

    return anno


def extractUK(anno):

    print("[EXTRACT]")
    # take all that have political entity UK
    entityMask = (anno[POLITICAL_ENTITY] == 'United Kingdom')

    # that are 1240K captured
    sourceMask = (anno['Data_source'].str.contains("1240[k|K]"))

    # and have PASS for the filter column
    # we take strict PASS only here
    passMask = (anno['ASSESSMENT'] == 'PASS')

    # remove relatives
    # noRelativeMask = numpy.isin(anno['Family_ID'], NO_RELATIVE_ENTRIES)
    noRelativeMask = (anno.Family_ID.str.contains(r"[N|n]o\s*[r|R]elative")) | (anno.Family_ID == "..")

    # take these individuals
    finalMask = entityMask & sourceMask & passMask & noRelativeMask

    print(f"{numpy.sum(finalMask)} individuals in UK that pass filter")

    # time mask
    ## convert to numbers
    # presentMask = (anno['Date_mean'] == 0)
    ancientMask = (anno['Date_mean'].astype(int) > 0)
    post4500Mask = (anno['Date_mean'].astype(int) <= 4500)

    print(f"\t{numpy.sum(finalMask & ancientMask)} in past")
    print(f"\t{numpy.sum(finalMask & ancientMask & post4500Mask)} more recent than 4500 YBP")

    finalIndividuals = anno.loc[finalMask & ancientMask & post4500Mask]

    print("[EXTRACT_DONE]")

    return finalIndividuals


def writeAnnotation(annoUK, OUTPUT_FILE):
    # just write it
    annoUK.to_csv(OUTPUT_FILE, sep='\t', index=False)


def main():
    # load the file
    anno = loadAnnotationFile(ANNO_FILE)

    # extract the UK individuals
    annoUK = extractUK(anno)

    # remove dups
    print('Removing duplicates...')
    before_nrow = annoUK.shape[0]
    annoUK = annoUK.drop_duplicates(subset=['Master_ID'])
    after_nrow = annoUK.shape[0]
    print(f'Removing {before_nrow - after_nrow} duplicate samples. {after_nrow} samples left.')

    # and write the new table file
    writeAnnotation(annoUK, OUTPUT_FILE)


if __name__ == "__main__":
    main()
