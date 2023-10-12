"""Read vcf, info, and inds file, output filtered count file, table of sample info, and list of snps

Usage: python %prog <vcf_file> <info_file> <inds_file> <snps_file> <filter1> <filter2> ... <filterN>

`filters` is a string of one or more (comma-delimited) filter abbr.'s below:
* ".anno" based filters:
    - "SG": Only take shot-gun sequenced samples
    - "noSG": Exclude shot-gun sequenced samples
    - "DG": Only take diploid samples
    - "noDG": Exclude diploid samples
    - "1240K": Only consider samples
    - "noModern": Exclude contemporary (t=0) samples
    - "noFam": Exclude individuals with family members in the record
    - "from[5678]to[1234]YBP": Only consider samples dated to be younger than 5678YBP and older than 1234YBP
* VCF-based filters:
    - "obsSamp.05num": total number of observed GT calls is no less than 5% of total number of samples
    - "sampTimeKge2": Only take samples with more than two (`ge` 2) sampling time points that has more than 1 observations
    - "MAF0.01": Minor allele frequency of all (pooled) samples no less than 0.01
"""
import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning)
# actual imports
import sys, os, time
import numpy as np
import pandas as pd
import pandas.api.types as ptypes
# local import
from diplo_locus.utility import parse_vcf_input, parse_ind_arg
# regex preps
import re

int_regex = re.compile(r'[^\.](\d+)')
num_regex = re.compile(r'\d+\.?\d*')
frac_regex = re.compile(r'([0]{0,1}\.\d+)')
float_regex = re.compile(r'\d+\.?\d*e?\d+')


def _split_locus_name(locus_name):
    split_names = locus_name.split("_")
    if len(split_names) > 3:
        Chr, physPos = split_names[0], split_names[1]
        rsID = '_'.join(split_names[2:])
    elif len(split_names) == 3:
        Chr, physPos, rsID = split_names
    else:
        return False
    return int(Chr), int(physPos), rsID


def v54_vcf_to_counts(vcffile, infoTableOrFile, gen_time=1., inds=None):
    """Convert a VCF from v54 parsing pipeline to a pd.DataFrame of allele counts
    :param vcffile: VCF(.gz) file name
    :param infofile: info table file name. Must have one "ID" column and one column named {"Gen_Ago", "Gen", "Time_Ago", "Time"}
    :param inds: list of sample IDs (subset of those in infofile) to consider
    :param gen_time: Number of time units per generation. Default is 1. Needed when info table presents  "Time" or "Time_Ago"
    :return: a `pd.DataFrame` with columns ["locus_name", "Chr", "physPos", "rsID", "d1", "n1", ... , "dK", "nK"]
    """
    # read files
    if isinstance(inds, str):
        ind_list = parse_ind_file(inds)
    else:
        assert isinstance(inds, list), f'`inds` is either a filename or a list-like object, not {inds}.'
        ind_list = inds

    if isinstance(infoTableOrFile, (list, pd.DataFrame, dict)):
        ## steal some codes from utility.py
        # make sure it's df
        infoTable = pd.DataFrame(infoTableOrFile)
        # assert 'Time_Ago' in infoTable.columns, infoTable.columns
        # assert 'ID' in infoTable, infoTable.columns
        # check for missing value
        if np.any(infoTable["Time_Ago"].isnull()):
            infoTable = infoTable[~infoTable["Time_Ago"].isnull()]
        # t0 needs to be the same unit too
        # if t0 is None:
        t0 = max(infoTable["Time_Ago"])
        # infoTable["Gen"] = np.array((t0 - np.array(infoTable["Time_Ago"])) / gen_time, dtype=int)
        # convert Time_Ago
        infoTable['Gen'] = (t0 - infoTable['Time_Ago']) / gen_time
        # just in case
        if np.any(infoTable["Gen"].isnull()):
            infoTable = infoTable[~infoTable["Gen"].isnull()]
        # print(f't0={t0}, {infoTable["Time_Ago"].describe()}, gen_time={gen_time}') #\n{list(infoTable["Time_Ago"])}
        # algorithm only accept integer generations for now
        infoTable['Gen'] = infoTable['Gen'].round().astype(int)
        # group by gen
        Pools = list(infoTable.groupby(by="Gen")["ID"])
        # pandas algorithm ensures that time_points are sorted
        time_points, sample_pools = map(list, zip(*Pools))
        sampleTimes = np.array([0] + time_points)
        sampleIDpools = [list(pool) for pool in sample_pools]
    else:
        assert os.path.isfile(infoTableOrFile)
        sampleTimes, sampleIDpools = read_info_file(infoTableOrFile, inds=ind_list, gen_time=gen_time)

    K = len(sampleTimes) - 1
    # now we can also know the number of samples per timepoint
    n_max = [len(pool) for pool in sampleIDpools]
    # locus_names, samples, sampleSizes = parse_vcf_input(vcffile, sampleIDpools, pseudo=pseudo, sampleTimes=sampleTimes)
    # parse_vcf_input(vcffile: str, samplePools: list, sampleTimes=None, snps=None, inds=None,
    #                     forced_haps: list or set=[], forced_dips: list or set=[], minMAF: float=0.)
    Parsed_vcf = parse_vcf_input(vcffile, sampleIDpools, sampleTimes=sampleTimes, inds=ind_list)
    locus_names, samples, sampleSizes, trimmedSampleTimes, het_flags, trimmedIDpools = Parsed_vcf
    # sanity check
    trimmed = []
    kept = 0
    for k, time_point in enumerate(sampleTimes[1:]):
        if time_point not in trimmedSampleTimes:
            trimmed.append(k)
        else:
            kept += 1
    assert (len(trimmed) + kept) == K, f'len(trimmed)={len(trimmed)}, kept={kept}, num_time_points={K}'
    # matrix shape should match too
    assert samples.shape[1] == len(trimmedSampleTimes) - 1
    assert samples.shape == sampleSizes.shape
    # update other related variables
    sampleTimes = trimmedSampleTimes
    sampleIDpools = trimmedIDpools
    K = len(sampleTimes) - 1
    n_max = [len(pool) for pool in sampleIDpools]

    # sanity check
    assert np.all(samples <= sampleSizes), np.sum(samples > sampleSizes)
    maxObsSampleSizes = np.max(sampleSizes, axis=0)
    size_check = np.where(n_max >= maxObsSampleSizes)
    # assert np.all(size_check), f'{sum(size_check)} locus; n_max={np.array(n_max)[~size_check]}\nobs_size={np.array(maxObsSampleSizes)[~size_check]}'
    size_check_failed = np.where(n_max < maxObsSampleSizes)
    assert np.all(
        n_max >= maxObsSampleSizes), f'{sum(size_check)} locus; n_max={np.array(n_max)[size_check_failed]}\nobs_size={np.array(maxObsSampleSizes)[size_check_failed]}'

    # turn them into df
    samples = pd.DataFrame(samples, columns=[f'd{i}' for i in range(1, K + 1)])
    samples['locus_name'] = locus_names
    # samples['Chr'], samples['physPos'], samples['rsID'] = _split_locus_name(locus_names)
    split_names = samples['locus_name'].apply(_split_locus_name)
    # print(type(split_names), split_names.shape, split_names[:5])
    samples['Chr'], samples['physPos'], samples['rsID'] = zip(*split_names)
    samples.physPos = samples.physPos.astype(int)
    # samples.rsID = samples.rsID.astype(str)
    # print(samples.describe())
    print(samples.shape)
    sampleSizes = pd.DataFrame(sampleSizes, columns=[f'n{i}' for i in range(1, K + 1)])
    print(sampleSizes.shape)
    # merge the two
    Samples = pd.concat([samples, sampleSizes], axis=1)
    # print(Samples.shape)
    # reorder the columns
    allele_count_header = ([''] * int(2 * K))
    allele_count_header[::2] = [f'd{i}' for i in range(1, K + 1)]
    allele_count_header[1::2] = [f'n{i}' for i in range(1, K + 1)]
    Samples = Samples[["locus_name", "Chr", "physPos", "rsID"] + allele_count_header]
    # print(Samples.describe())
    # print(Samples.shape)
    assert np.all(
        np.array(Samples[[f'd{i}' for i in range(1, K + 1)]]) <= np.array(Samples[[f'n{i}' for i in range(1, K + 1)]]))

    return Samples, K, n_max, sampleTimes


# filter_args is a list of strings (already split)
## return: infoTable
def filter_v54_anno_inds(infofile, country: str, ID_col: str, time_col: str, from_ybp: float, to_ybp: float,
                         SG_filter=None, DG_filter=None,
                         fam_filter=None, only_1240k=None, modern_filter=None, access_filter=None):
    # read info table
    infoTable = pd.read_csv(infofile, sep="\t",
                            dtype=str)  # , comment="$" basically a separator that has not appeared in the doc
    # this is v54.1-specific
    infoTable.columns = ['GenID', 'MasterID', 'Skeletal_code', 'Skeletal_element',
                         'Pub_year', 'Publication', 'Dating_method', 'MeanYBP', 'SdYBP', 'Full_Date',
                         'Death_Age', 'GroupID', 'Locality', 'Country', 'Lat', 'Long', 'Pulldown_method', 'Data_source',
                         'Num_library', 'Coverage_1240k', 'SNPhit_1240k', 'SNPhit_HO', 'Mol_Sex', 'FamID',
                         'YhapGroup_term_mut',
                         'YhapGroup_ISOGG', 'mtDNA_Coverage', 'mtDNA_HapGroup', 'mtDNA_match', 'Damage_rate_1stNt',
                         'Sex_ratio', 'Library_type', 'Libraries', 'Assessment', 'Assessment_warnings', 'Unnamed']
    # only keep a subset of them
    infoTable = infoTable[['GenID', 'MasterID', 'Pub_year', 'Publication', 'Dating_method', 'MeanYBP', 'SdYBP',
                           'Death_Age', 'GroupID', 'Locality', 'Country', 'Lat', 'Long', 'Pulldown_method',
                           'Data_source',
                           'Num_library', 'Coverage_1240k', 'SNPhit_1240k', 'SNPhit_HO', 'Mol_Sex', 'FamID',
                           'Sex_ratio', 'Library_type', 'Libraries', 'Assessment', 'Assessment_warnings']]
    tag = ''
    if country is not None:
        # assert isinstance(country, str), f'country {country}, type {type(country)}'
        if isinstance(country, str):
            if country == "UK":
                infoTable = infoTable[infoTable.Country.str.contains('United Kingdom')]
                tag += '_UK'
            elif sum(infoTable.Country.str.contains(country, regex=False)) == 0:
                assert np.any(infoTable.Locality.str.contains(country, regex=True)), \
                    f'{country} not in \"Locality\" or \"Country\" column.\n{infoTable.Country.value_counts()}\n{infoTable.Locality.value_counts()}'
                infoTable = infoTable[infoTable.Locality.str.contains(country, regex=False)]
                tag += f'_{country}'
            else:
                infoTable = infoTable[infoTable.Country.str.contains(country, regex=False)]
                tag += f'_{country}'
        elif isinstance(country, list):
            assert np.all(np.array([isinstance(cntry, str) for cntry in
                                    country])), f'All elements for country arg need to be str. {country}'
            infoTable = infoTable[infoTable.Country.isin(country)]
            tag += f'_{"-".join(country)}'
    else:
        print('Will not filter based on countries or localities.')
        print(infoTable.Country.value_counts().index)
        print(infoTable.Locality.value_counts().index)

    # print(infoTable.head())
    if SG_filter is not None:
        if SG_filter:
            ## some ones might be missing if only filter SG tag, e.g. "S_Adygei-2.DG", "AITI_87_d" <-- resequenced
            # only keep shot gun seqs
            print('Only use samples of shot-gun sequences')
            ## some .DG were also from shotgun sequenced
            ## this step will not remove diploids
            before_nrow = infoTable.shape[0]
            infoTable = infoTable[infoTable.Data_source.str.contains("[S|s]hotgun")]
            after_nrow = infoTable.shape[0]
            print(f'Only keep {after_nrow} shotgun-sequenced samples out of {before_nrow} total.\n')
            tag += '_SG'
        else:
            print('Remove samples of shot-gun sequences')
            before_nrow = infoTable.shape[0]
            infoTable = infoTable[~infoTable.Data_source.str.contains("[S|s]hotgun")]
            # the two criteria aren't completely overlapping. Need to do both
            # infoTable = infoTable[~infoTable.GenID.str.contains('SG')]
            infoTable = infoTable[~infoTable.GenID.str.endswith("SG")]
            after_nrow = infoTable.shape[0]
            print(f'{after_nrow} samples left after removing {before_nrow - after_nrow} shotgun sequences.\n')
            tag += '_noSG'
    else:
        print('Will not filter samples by sequencing method.')
        print(infoTable.Data_source.value_counts())

    if only_1240k is not None:
        if only_1240k:
            print('Only use samples from 1240k pull down')
            ## some .DG were also from shotgun sequenced
            ## this step will not remove diploids
            before_nrow = infoTable.shape[0]
            infoTable = infoTable[infoTable.Data_source.str.contains("1240[k|K]")]
            after_nrow = infoTable.shape[0]
            print(f'Only keep {after_nrow} 1240k samples out of {before_nrow} total.\n')
            tag += '_1240K'

    if DG_filter is not None:
        if DG_filter:
            print('Only use samples with diploid genomes')
            # these ones are all labeled with .DG tag
            before_nrow = infoTable.shape[0]
            infoTable = infoTable[infoTable.Data_source.str.contains("[D|d]iploid")]
            # sanity check
            assert np.all(infoTable.GenID.str.contains("DG")), infoTable[~infoTable.GenID.str.contains("DG")]
            after_nrow = infoTable.shape[0]
            print(f'Only keep {after_nrow} diploid samples out of {before_nrow} total.\n')
            tag += '_DG'
        else:
            print('Remove samples with diploid genomes')
            before_nrow = infoTable.shape[0]
            # everything else should be alright, even if some don't satisfy /.DG$/
            infoTable = infoTable[~infoTable.Data_source.str.contains("[D|d]iploid")]
            after_nrow = infoTable.shape[0]
            print(f'{after_nrow} samples left after removing {before_nrow - after_nrow} diploid sequences.\n')
            tag += '_noDG'
    else:
        print('Will not filter samples by ploidy.\n')
        print(infoTable.Data_source.value_counts())

    if fam_filter is not None:
        if not fam_filter:
            print('Only keep samples with no recorded relatives')
            before_nrow = infoTable.shape[0]
            ## there're 41 samples from Nepal(23), Japan(10), and S.Korea(8) without info for FamID ("..")
            ## All from 2 pub: RobbeetsNingNature2021, LiuJeongNatComm2022
            infoTable = infoTable[
                (infoTable.FamID.str.contains(r"[N|n]o\s*[r|R]elative")) | (infoTable.FamID == "..")]
            after_nrow = infoTable.shape[0]
            print(
                f'{after_nrow} samples left after removing {before_nrow - after_nrow} samples with known relatives.\n')
            tag += '_noFam'
    else:
        print('Will not filter samples based on records of relatives. Existing family IDs:')
        print(infoTable.FamID.value_counts().shape)

    if modern_filter is not None:
        if modern_filter:
            print('Remove contemporary (t=0) samples')
            before_nrow = infoTable.shape[0]
            if type(list(infoTable.MeanYBP)[0]) is str:
                infoTable.MeanYBP = infoTable.MeanYBP.astype(int)
            infoTable = infoTable[infoTable.MeanYBP != 0]
            after_nrow = infoTable.shape[0]
            print(
                f'{after_nrow} samples left after removing {before_nrow - after_nrow} contemporary (t=0) samples.\n')
            tag += '_noModern'

    if access_filter is not None or access_filter != 'none':
        if access_filter == 'loose':
            print(f'Only keep samples that loosely pass Accessment filter (ie, _contain_ PASS)')
            before_nrow = infoTable.shape[0]
            infoTable = infoTable[infoTable.Accessment.str.contains('PASS')]
            after_nrow = infoTable.shape[0]
            print(
                f'{after_nrow} samples left after removing {before_nrow - after_nrow} samples without \"PASS\".\n')
            tag += '_loosePASS'
        elif access_filter == 'strict':
            print(f'Only keep samples that loosely pass Accessment filter (ie, _contain_ PASS)')
            before_nrow = infoTable.shape[0]
            infoTable = infoTable[infoTable.Assessment == 'PASS']
            after_nrow = infoTable.shape[0]
            print(
                f'{after_nrow} samples left after removing {before_nrow - after_nrow} samples without an exact \"PASS\".\n')
            tag += '_strictPASS'
        else:
            print(f'Cannot recognize the argument `--PASS={access_filter}`')
            print(infoTable.Accessment.value_counts())

    # make sure ybps are of type number
    if not ptypes.is_numeric_dtype(infoTable.MeanYBP):
        infoTable.MeanYBP = pd.to_numeric(infoTable.MeanYBP)

    if from_ybp is not None:
        print(f'Only consider samples younger than {from_ybp} YBP.')
        # assert isinstance(from_ybp, (int, float)), type(from_ybp)
        if not isinstance(from_ybp, (int, float)):
            try:
                from_ybp = int(from_ybp)
            except ValueError:
                from_ybp = float(from_ybp)
        before_nrow = infoTable.shape[0]
        infoTable = infoTable[infoTable.MeanYBP <= from_ybp]
        after_nrow = infoTable.shape[0]
        print(
            f'{after_nrow} samples left after removing {before_nrow - after_nrow} samples older than {from_ybp} YBP.\n')
        tag += f'_from{from_ybp:g}'
    else:
        print('No upper cap for sample age.')

    if to_ybp is not None:
        print(f'Only consider samples older than {to_ybp} YBP.')
        # assert isinstance(to_ybp, (int, float)), type(to_ybp)
        if not isinstance(to_ybp, (int, float)):
            try:
                to_ybp = int(to_ybp)
            except ValueError:
                to_ybp = float(to_ybp)
        before_nrow = infoTable.shape[0]
        infoTable = infoTable[infoTable.MeanYBP > to_ybp]
        after_nrow = infoTable.shape[0]
        print(
            f'{after_nrow} samples left after removing {before_nrow - after_nrow} samples younger than {to_ybp} YBP.\n')
        tag += f'_to{to_ybp:g}'
    else:
        print('No lower cap for sample age.')

    print(infoTable.MeanYBP.describe())

    # remove duplicates
    print('Removing duplicates...')
    before_nrow = infoTable.shape[0]
    infoTable = infoTable.drop_duplicates(subset=['MasterID'])
    after_nrow = infoTable.shape[0]
    print(f'Removing {before_nrow - after_nrow} duplicate samples. {after_nrow} samples left.')
    # now we're finally done with everything (for now)
    # rename the column for future convenience <-- no longer necessary
    infoTable = infoTable.rename(columns={'MeanYBP': time_col, 'GenID': ID_col})
    # print(f'{infoTable.shape[0]} individuals/samples')

    ## exclude the "_" at the beginning of the tag
    return infoTable, tag[1:]


def _format_to_dotless_e(number):
    # Format the number in "e" notation and convert it to a float
    formatted_number = "{:e}".format(float(number))
    # exponent = float(formatted_number.split("e")[1])
    e_part1, e_part2 = formatted_number.split("e")
    # remove zeros at the end:
    e_part2 = int(e_part2)
    if number >= 1:
        while not np.isclose(float(e_part1) - int(e_part1.split(".")[0]), 0):
            part1_update = float(e_part1) * 10
            e_part2 -= 1
            e_part1 = str(part1_update)
        # print(e_part1, e_part2)
        # only look at integer part
        while e_part1.endswith("0"):
            e_part1 = e_part1[:-1]
            e_part2 -= 1
            # print(e_part1, e_part2)
    else:
        # only look at decimal part
        while e_part1.startswith("0"):
            e_part1 = e_part1[1:]
            e_part2 += 1
            # print(e_part1, e_part2)
    # Convert the number to an integer to remove the decimal point and trailing zeros
    # formatted_integer = int(float(e_part1))
    return "{:d}e{:d}".format(int(float(e_part1)), int(e_part2) + 1)


def filter_snps(Samples, K, n_max, minK, missingness: float, minMAF: float, chroms=None, pos_arg=None):
    """Remove loci whose data don't pass the filters. Input: allelic count matrix"""
    sample_cols = [f'd{i}' for i in range(1, (K + 1))]
    size_cols = [f'n{j}' for j in range(1, (K + 1))]
    tag = ''

    ## parse chrom
    if chroms is not None:
        Chroms = chroms.split(',')
        Chroms = [int(c) for c in Chroms if len(c) > 0]
        Samples = Samples[Samples.Chr.isin(Chroms)]
        print(f'Keeping {Samples.shape[0]} SNPs on chromosome(s) {Chroms}')
        tag += f'_chr{chroms}'
    else:
        Chroms = None
    ## then position
    if pos_arg is not None:
        assert (len(Samples.Chr.values) == 1) or (Chroms is not None), f'Chroms: {Chroms}, position range {pos_arg}'
        numbers = float_regex.findall(pos_arg)
        if len(numbers) == 2:
            left_pos, right_pos = sorted([float(p) for p in numbers])
            # pos_range = (left_pos, right_pos)
            tag += f'_{_format_to_dotless_e(left_pos)}-{_format_to_dotless_e(right_pos)}'
        elif len(numbers) == 1:
            # need to tell which one
            assert "-" in pos_arg, pos_arg
            pos_range = pos_arg.split("-")
            if int_regex.search(pos_range[0]):
                left_pos = int(pos_range[0])
                right_pos = None
            elif int_regex.search(pos_range[1]):
                left_pos = None
                right_pos = int(pos_range[1])
            else:
                print(f'cannot parse the position flag: {pos_arg}.')
                sys.exit(1)
            # pos_range = (left_pos, right_pos)
        else:
            print(f'cannot parse the position flag: {pos_arg}.')
            sys.exit(1)
        before = Samples.shape[0]
        Samples = Samples[(Samples.physPos >= left_pos) & (Samples.physPos <= right_pos)]
        after = Samples.shape[0]
        print(f'Keeping {after} SNPs (from {before}) positioned between {int(left_pos)} to {int(right_pos)}')

    assert np.all(np.array(Samples[sample_cols]) <= np.array(Samples[size_cols]))
    Samples.loc[:, 'pooled_d'] = Samples[sample_cols].sum(axis=1)
    Samples.loc[:, 'pooled_n'] = Samples[size_cols].sum(axis=1)

    ## now check out number of times observed
    if minK > 0:
        Samples['times_obs'] = (Samples[size_cols] > 0).sum(axis=1)
        before = Samples.shape[0]
        Samples = Samples[Samples.times_obs >= minK]
        after = Samples.shape[0]
        print(f'Remove {before - after} SNPs with <{minK} times/pools of observations. {after} SNPs left.')
        tag += f'_minK{minK:d}'

    # filter for missingness <--- float in (0, 1)
    before = Samples.shape[0]
    # print(Samples.pooled_n.dtype)
    Samples = Samples[(Samples.pooled_n > missingness * sum(n_max))]
    tag += f'_minObs{missingness:g}'
    after = Samples.shape[0]
    print(f'{after} SNPs (out of {before}) passed missingness filter (:= pooled_n > sum(n_max)*{missingness:g} ).')

    ## lastly: MAF
    assert 0 <= minMAF < 0.5, f'Invalid minor allele frequency cutoff: {minMAF}'
    tag += f'_MAF{str(minMAF)[1:]}'
    # get maf
    Samples.loc[:, 'pooled_d'] = Samples.pooled_d.where(
        Samples.pooled_d.copy() <= 0.5 * Samples.pooled_n.copy(),
        (Samples.pooled_n.copy() - Samples.pooled_d.copy()))
    before = Samples.shape[0]
    Samples = Samples[Samples.pooled_n > 0]
    after = Samples.shape[0]
    print(f'Remove {before - after} sites with zero samples')
    ## actual MAF filter
    before = Samples.shape[0]
    Samples['pooledMAF'] = Samples.pooled_d / Samples.pooled_n  # .replace(0, np.nan) <-- there shouldn't be n=0 anymore
    Samples = Samples[Samples.pooledMAF > minMAF]
    after = Samples.shape[0]
    print(f' Remove {before - after} SNPs whose pooled MAF <= {minMAF}. {after} left.')
    # remove aux column?
    Samples = Samples.drop(columns=['pooled_d', 'pooled_n'])

    # done going through all filter args, return df
    ## exclude the "_" at the beginning of the tag
    return Samples, tag[1:]


def subset_VCF(inVCFname, outVCFname, inds=None, snps=None, notes=''):
    if inVCFname.lower().endswith('.gz'):
        input_vcf = gzip.open(inVCFname, 'rt')
    elif inVCFname.lower().endswith('.vcf'):
        input_vcf = open(inVCFname, 'r')
    else:
        raise ValueError(f'Unrecognized filetype: .{inVCFname.split(".")[-1]}')

    output_vcf = open(outVCFname, 'w')
    # setting up
    from diplo_locus.utility import _get_ids_from_vcf_header
    twohash_regex = re.compile(r'^##')
    onehash_regex = re.compile(r'^#[^#](.*)')

    in_l = input_vcf.readline()
    header=[]
    sampleIndexPool = []
    while in_l != "":
        if twohash_regex.search(in_l):
            assert len(header) == 0, len(header)
            output_vcf.write(in_l)
        elif onehash_regex.search(in_l):
            # it's header line!
            header = in_l.strip().split("\t")
            if inds is None:
                sampleIndexPool = range(9, len(header))
                output_header = header
            else:
                assert isinstance(inds, list), type(inds)
                sampleIndexPool = _get_ids_from_vcf_header(header, inds)
                output_header = header[:9] + list(inds)
            # write notes
            if notes != '':
                output_vcf.write(notes)
            # write header
            output_vcf.write("\t".join(output_header) + "\n")
        else:
            assert len(header) > 9, header
            l = in_l.strip().split("\t")
            assert len(l) > 9, ValueError('Please check the vcf file format.')
            assert sampleIndexPool != [], f'samplePools = {inds},\nheader = {header}'
            # pass if not in the snp list
            if snps is not None:
                assert isinstance(snps, list), f'type(snps)={type(snps)}'
                if l[2] not in snps:
                    in_l = input_vcf.readline()
                    continue
            # assume all snps given are pre-filtered & are thus bi-allelic and all
            # directly write out everything
            outline = l[:9] + [l[samp_id] for samp_id in sampleIndexPool]
            output_vcf.write("\t".join(outline) + "\n")
        # that's it!
        in_l = input_vcf.readline()
    print('Done.')


import argparse


def args_parser(parser: argparse.ArgumentParser):
    inputs = parser.add_argument_group('Input')
    inputs.add_argument('--vcf', dest='vcffile', required=True, help='VCF files to be parsed.')
    inputs.add_argument('--anno', dest='annofile', required=True,
                        help='\".anno\" file as presented in the v54.1 1240K database.')
    inputs.add_argument('--gen_time', dest='gen_time', default=26.9, type=float,
                        help='Average generation time, in years, to be used for converting years to generations.'
                             'Default value is 26.9 years (Wang et al. 2023; https://doi.org/10.1126/sciadv.abm7047)')

    outputs = parser.add_argument_group('Output Options')
    outputs.add_argument('--out', '-o', '--outname', dest='outprefix', default='v54.1_1240K', #required=True,
                         help='path and prefix to output files.')
    outputs.add_argument('--writeVCF', '--outputVCF', dest='writeVCF', default=False,
                         action='store_true', help='Option to write out the filtered VCF file. Default is False.')
    outputs.add_argument('--writeCounts', '--outputCounts', dest='writeCounts', action=argparse.BooleanOptionalAction, default=True,
                         help='Choose whether the allele count file will be written.')

    sFilters = parser.add_argument_group('Sample-based Filters')
    # for file-reading
    sFilters.add_argument('--ID_col', dest='IDcolName', default='GenID',
                          help='Name of the column in info table for ID names of the sample (as shown in VCF). Default is \"ID\".')
    sFilters.add_argument('--time_ago_col', dest='TimeAgoColName', default='MeanYBP',
                          help='Name of the column in info table for each sample\'s sampling times (backward; ascending'
                               ' from present to past). Default name is \"MeanYBP\", unit is years before present (:=1950).')

    # for detail criteria
    sFilters.add_argument('--inds', dest='inds', default=None,
                          help='List (or file name) of sample IDs to include. '
                               'Default is to include all samples shared between vcf & info files.')
    sFilters.add_argument('--country', '--place', '--locality', '--region', dest='where', default=None,
                          help='Key word in \"Political Entity\" or \"Locality\" column in the .anno file.')
    sFilters.add_argument('--SG', '--Shotgun', dest='SG', default=None, action=argparse.BooleanOptionalAction,
                          help='Option to include (\"True\") or exclude (\"False\") shotgun-sequenced samples,'
                               ' as indicated in the \"Data source\" column. Default is to ignore such info.')
    sFilters.add_argument('--DG', '--Diploid', dest='DG', default=None, action=argparse.BooleanOptionalAction,
                          help='Option to include (\"True\") or exclude (\"False\") samples with diploid variant calls,'
                               ' as indicated in the \"Data source\" column. Default is to ignore such info.')
    sFilters.add_argument('--1240K_only', '--1240k_only', dest='only1240K', default=None, action=argparse.BooleanOptionalAction,
                          help='Option to only include samples sequenced via 1240K captures. That is, \"Data source\"'
                               '_is exactly_ \"1240[K|k]\".')
    sFilters.add_argument('--PASS', dest='PASS', choices=['none', 'loose', 'strict'], default='loose',
                          help='Options on the requirement for quality filters.\"none\": ignore such info;\n'
                               '\"loose\" (default): samples whose \"Assessment\" column _contains_ \"PASS\";\n'
                               '\"strict\": samples whose \"Assessment\" column _is exactly_ \"PASS\".')
    sFilters.add_argument('--Fam', dest='Fam', default=None, action=argparse.BooleanOptionalAction,
                          help='Option to include (\"True\") or exclude (\"False\") samples with known relatives recorded, '
                               'as indicated in the \"Family ID\" column. Default is to ignore such info.')
    sFilters.add_argument('--Modern', default=None, dest='Mod', action=argparse.BooleanOptionalAction,
                          help='Option to include (\"True\") or exclude (\"False\") samples with known relatives recorded, '
                               'as indicated in the \"Family ID\" column. Default is to ignore such info.')
    sFilters.add_argument('--fromYBP', '--from', '--sinceYBP', '--since', dest='from_year', default=None,
                          help='The oldest time (including itself), in years before present, based on the sample\'s mean'
                               ' calibrated YBP, to consider. Default is to take all.')
    sFilters.add_argument('--toYBP', '--to', dest='to_year', default=None,
                          help='The latest time (excluding itself), in years before present, based on the sample\'s mean'
                               ' calibrated YBP, to consider. Default is to take all.')

    vFilters = parser.add_argument_group('VCF-based Filters')
    vFilters.add_argument('--force_hap', dest='forced_haps', default='none',
                          help='Comma-separated string or a plain-text file that lists IDs (matching VCF column names) to be considered as haploids even though some may present diploid genotype calls. If \"all\", all samples will be considered as haploids, whose genotypes will be counted as matching haplotype (i.e. half the alleles). Force quit if any specified individual has heterozygote variant or any unmentioned diploid-formatted samples lack heterozygote variant calls. Default is \"none\".')
    vFilters.add_argument('--force_dip', dest='forced_dips', default='none',
                          help='Comma-separated string or a plain-text file that lists samples IDs (separated by comma) to be considered as diploids even though some may have loci with only haplotype calls. If \"all\", all samples will be considered as diploids, whose haploid GTs will be counted as matching homozygote diploid GTs (i.e. double-count the allele). Default is \"none\".')
    vFilters.add_argument('--snps', dest='snps', default=None,
                          help='List (or file name) of SNPs to include. Default is to include all qualifying SNPs in the VCF.')
    vFilters.add_argument('--chr', '-c', dest='chrom',
                          help='Name of the chromosome to extract SNPs from. Must match values in the VCF. '
                               'Default is to include all qualifying SNPs in the VCF.')
    vFilters.add_argument('--pos', dest='posRange',
                          help='Range of physical positions (\"<from>-<to>\") to consider on the proposed'
                               ' chromosome. Default is to include all qualifying SNPs in the VCF.')
    missingness = vFilters.add_mutually_exclusive_group()
    missingness.add_argument('--max_missing', dest='maxMissing', type=float, default=1.,
                             help='For each SNP, the max allowed proportion of samples to have missing data.')
    missingness.add_argument('--min_observed', dest='minObs', type=float, default=0.,
                             help='For each SNP, the min allowed proportion of samples to have valid observations (not missing).')
    vFilters.add_argument('--minK', dest='minK', default=0, type=int,
                          help='For each SNP, the minimal number of time points to have non-zero observations in.')
    vFilters.add_argument('--minMAF', dest='minMAF', type=float, default=0.,
                          help='Minimum threshold (non-inclusive) for minor allele frequencies in the pooled samples. '
                               'SNPs with sub-threshold frequencies will not be considered for analyses. Default is 0.')
    return parser


def main(CLargs):
    parser = argparse.ArgumentParser()
    parser = args_parser(parser)
    args = parser.parse_args(CLargs)
    # print(args)

    # read anno
    notes = f'## Read sample info from {args.annofile}\n'
    infoTable, samp_tag = filter_v54_anno_inds(args.annofile, country=args.where, ID_col=args.IDcolName,
                                               time_col=args.TimeAgoColName, SG_filter=args.SG, DG_filter=args.DG,
                                               fam_filter=args.Fam, only_1240k=args.only1240K,
                                               modern_filter=args.Mod, from_ybp=args.from_year,
                                               to_ybp=args.to_year, access_filter=args.PASS)
    # some sanity check
    assert (args.IDcolName in infoTable.columns) and (args.TimeAgoColName in infoTable.columns), infoTable.columns
    IDcol = args.IDcolName
    TimeAgoCol = args.TimeAgoColName
    # "translate" column names into ones the subsequent codes will recognize
    # infoTable = infoTable.rename(columns={args.IDcolname: "ID", args.time_ago_col: "Time_Ago"})
    print(f'{infoTable.shape[0]} individuals/samples')
    # write them out
    print(f'Writing info table to {args.outprefix}_{samp_tag}.info')
    infoTable.to_csv(f'{args.outprefix}_{samp_tag}.info', index=False, sep="\t")

    # get ind_list
    if args.inds is not None:
        # or args.snp_subset is not None:
        ind_list = parse_ind_arg(args.inds, "inds")
        notes += f'## Only consider {len(ind_list)} samples from `--inds`\n'
        infoTable = infoTable[infoTable[IDcol].isin(ind_list)]
    else:
        ind_list = list(infoTable[IDcol])

    # now we get gen_time and group them by generations
    ## default t0 is the first/oldest sampling time
    t0 = max(infoTable[TimeAgoCol])
    print(
        f'Setting the earliest sampling time t0={t0} as generation 0, adopting average generation time of {args.gen_time}.')
    infoTable['Gen'] = (t0 - infoTable[TimeAgoCol]) / args.gen_time
    # in case there's nan
    if np.any(infoTable.Gen.isnull()):
        print(f'{sum(infoTable.Gen.isnull())} samples with gen.isnull()')
        print(infoTable.loc[infoTable.Gen.isnull(), [IDcol, TimeAgoCol, 'Gen']])

    # group by gen
    Pools = list(infoTable.groupby(by="Gen")[IDcol])
    # pandas algorithm ensures that time_points are sorted
    time_points, sample_pools = map(list, zip(*Pools))
    sampleTimes = np.array([0] + time_points)
    samplePools = [list(pool) for pool in sample_pools]

    # do vcf
    ## parse snps
    if args.snps is not None:
        snp_list = parse_ind_arg(args.snps, "snps")
        notes += f'## {len(snp_list)} loci from `--snps`\n'
    else:
        snp_list = None

    # deal with ploidy
    ## read all GTs as-is unless flags specify
    all_IDs = {ID for pool in samplePools for ID in pool}
    if args.forced_haps != 'none' and args.forced_haps != "all":
        hap_IDs = parse_ind_arg(args.forced_haps)
        if len(hap_IDs) == 0:
            print('Invalid input for `--force_hap`:', args.forced_haps)
        else:
            hap_IDs = set(hap_IDs)
            print(
                f'{len(hap_IDs)} samples will be deliberately counted as haploids: {", ".join(list(hap_IDs))}')
    elif args.forced_haps == "all":
        hap_IDs = all_IDs
        print(f'All samples will be deliberately counted as haploids.')
    else:
        hap_IDs = {}

    if args.forced_dips != "none" and args.forced_dips != "all":
        dip_IDs = parse_ind_arg(args.forced_dips)
        if len(dip_IDs) == 0:
            print('Invalid input for `--force_dip`:', args.forced_dips)
        else:
            dip_IDs = set(dip_IDs)
            print(
                f'{len(dip_IDs)} samples will be deliberately counted as diploids: {", ".join(list(dip_IDs))}')
    elif args.forced_dips == "all":
        dip_IDs = all_IDs
        print(f'All samples will be deliberately counted as diploids.')
    else:
        dip_IDs = {}

    # sanity check
    if len(hap_IDs) > 0 and len(dip_IDs) > 0:
        assert not hap_IDs.isdisjoint(
            dip_IDs), f'Sample(s) {hap_IDs.intersection(dip_IDs)} cannot be forced to be diploid and haploid at the same time.'
    # all IDs should be subset of all_IDs
    # print(hap_IDs, dip_IDs)
    assert set(hap_IDs).issubset(all_IDs) and set(dip_IDs).issubset(
        all_IDs), f'Sample IDs specified in `--force_hap` and `--force_dip` must already be included in the VCF. {hap_IDs - all_IDs} {dip_IDs - all_IDs} are not included.'

    print(f'Loading variant counts from {args.vcffile}...')  # <-- leave MAF filter to next step
    ## note-to-self:
    ## parse_vcf_input(vcffile: str, samplePools: list, sampleTimes=None, snps=None, inds=None,
    ##                     forced_haps: list or set=[], forced_dips: list or set=[], minMAF: float=0.)
    ## return locus_names, samples, sampleSizes, sampleTimes, het_flags, sampleIndexPools
    Parsed_vcf = parse_vcf_input(args.vcffile, samplePools, sampleTimes, snp_list, ind_list,
                                 forced_haps=hap_IDs, forced_dips=dip_IDs)
                                 # missingness=(args.maxMissing, args.minObs),
    IDs, samples, sampleSizes, trimmedSampleTimes, het_flags, trimmedIndexPools = Parsed_vcf

    # in case it's trimmed
    if len(trimmedSampleTimes) != len(sampleTimes):
        assert len(trimmedSampleTimes) < len(
            sampleTimes), f'Before SNP-based trimming, sample times: {sampleTimes}\nAfter: {trimmedSampleTimes}'
        # sanity check
        trimmed = []
        kept = 0
        for k, time_point in enumerate(sampleTimes[1:]):
            if time_point not in trimmedSampleTimes:
                trimmed.append(k)
            else:
                kept += 1
        assert (len(trimmed) + kept) == samples.shape[1], f'len(trimmed)={len(trimmed)}, kept={kept}'
        # matrix shape should match too
        assert samples.shape[1] == len(trimmedSampleTimes) - 1
        assert samples.shape == sampleSizes.shape
        # update other related variables
        sampleTimes = trimmedSampleTimes

    # sampleIDpools = trimmedIndexPools sampleIDpools
    K = len(sampleTimes) - 1
    n_max = [len(pool) for pool in trimmedIndexPools]

    # sampleIDPools = [header[idx] for pool in trimmedSampleIndexPools for idx in pool] # <-- no use
    # make Samples dataFrame
    # print(K, samples.shape)
    samples = pd.DataFrame(samples, columns=[f'd{i}' for i in range(1, K + 1)], dtype=int)
    samples['locus_name'] = IDs
    # samples['Chr'], samples['physPos'], samples['rsID'] = _split_locus_name(locus_names)
    split_names = samples['locus_name'].apply(_split_locus_name)
    # print(type(split_names), split_names.shape, split_names[:5])
    samples['Chr'], samples['physPos'], samples['rsID'] = zip(*split_names)
    samples.physPos = samples.physPos.astype(int)
    # samples.rsID = samples.rsID.astype(str)
    # print(samples.describe())
    print(samples.shape)
    sampleSizes = pd.DataFrame(sampleSizes, columns=[f'n{i}' for i in range(1, K + 1)])
    print(sampleSizes.shape)
    # merge the two
    Samples = pd.concat([samples, sampleSizes], axis=1)
    # print(Samples.shape)
    # reorder the columns
    allele_count_header = ([''] * int(2 * K))
    allele_count_header[::2] = [f'd{i}' for i in range(1, K + 1)]
    allele_count_header[1::2] = [f'n{i}' for i in range(1, K + 1)]
    Samples = Samples[["locus_name", "Chr", "physPos", "rsID"] + allele_count_header]

    # ^ also made sure to have ID & Time_Ago columns
    print(f'{infoTable.shape[0]} samples after filtering')

    # filter snps
    # define missingness as the proportion left / not missing
    missingness = max([1 - args.maxMissing, args.minObs])
    print('maxMissing=', args.maxMissing, '\nminObserved=', args.minObs,
          '\nmissingness=', missingness)
    # assert 0 <= args.minMissing <= args.maxObs, f'Invalid (lower={args.minMissing}, upper={args.maxObs}) bound for missingness'
    # `missingness` be a tupple of 2 float
    filteredSamples, snp_tag = filter_snps(Samples, K, n_max, minK=args.minK, missingness=float(missingness),
                                           minMAF=args.minMAF, chroms=args.chrom, pos_arg=args.posRange)
    print(f'{filteredSamples} variants after filtering.')

    # writes out snps
    outname = f'{args.outprefix}_{samp_tag}_{snp_tag}'
    print(f'Writing filtered SNP IDs to {outname}.snps ...')
    filteredSamples.rsID.to_csv(outname + '.snps', index=False, header=False)

    # assemble preamble notes
    notes = '## Samp.Gen.Ago: ' + ', '.join(list(map(str, sampleTimes)))
    notes += '\n## Max.Samp.Size: ' + ', '.join(list(map(str, n_max))) + '\n'
    if args.writeVCF:
        outVCFname = outname + '.vcf'
        print(f'Writing filtered VCF to {outVCFname} ...')
        subset_VCF(args.vcffile, outVCFname, inds=ind_list, snps=snp_list, notes=notes)

    # write notes
    print(f'Writing temporal allele counts to {outname}.count ...')
    with open(outname + '.count', 'w') as count:
        count.write(notes)
    count.close()
    filteredSamples.to_csv(outname + '.count', sep="\t", index=False, na_rep='n/a', mode='a')


if __name__ == "__main__":
    main(sys.argv[1:])
    # temp_main_annoOnly()
