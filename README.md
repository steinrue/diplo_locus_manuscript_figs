
# Figures for _"`DiploLocus`: A lightweight toolkit for inference and simulation of time-series genetic data under general diploid selection"_

This repository contains replicable scripts for generating simulations and replicating figures as shown in the manuscript.

------------------
## Set up

In order for scripts in this repository to run, the Python version must be at least 3.8. The user must also have the [`diplo_locus` package installed](https://github.com/steinrue/diplo_locus/). Both the API and CLI will be used.

Download this repository with
```shell
# download
git clone https://github.com/steinrue/diplo_locus_manuscript_figs
# enter the dir
cd diplo_locus_manuscript_figs
```

All scripts in the instructions below are used only in Unix command-line interface.


<a id="toc"> </a>

## Table of Contents


 * [Fig. S1 - S7: Data simulation and likelihood-based inference](#DLsim)
 * [Fig. S8 - S12: `DiploLocus` Likelihood-based inference on `SLiM` simulations](#Slim)
 * [Fig. S13 - S14: Examining genotype data around *LCT* gene in UK population from AADR database](#LCT)
 * [Fig. 1 in main text: Combination of previous figures](#main_fig)


## Data simulation and likelihood-based inference using `diplo_locus` API
<a id="DLsim"> </a>

Navigate to the folder and create a folder to keep all the simulated data:
```shell
# navigate to the dir
cd figS1-S7_DLsims/
# create folder for simulated data
mkdir simulations/
```

### Data generation
To simulate replicates and compute their likelihoods, run
```shell
# Run analysis with variable s_AA: 
## Simulations initiated with a Watterson's standard neutral
## Usage: 
## python DL_sim2mlr_stdVar.py <outprefix> <nuet_samples_pkl> <init_distn> <minMAF> [seed]
python DL_sim2mlr_stdVar.py simulations/varS2 "none_yet" initDistn_watterson-N1e4.txt 0.05 87235

## Likewise, to simulate data from initial frequency of 0.01:
## Usage:
## python DL_sim2mlr_initFreq.py <outprefix> <neut_samples_pkl> <minMAF> [seed]
python DL_sim2mlr_initFreq.py simulations/varS2 "none_yet" 0.05 10827

# Run analyses with variable S_aA:
## Usage:
## python () <outprefix> <neut_samples_pkl> <init_distn> <minMAF> [seed]
python DL_sim2mlr_s1-var.py simulations/varS1 "none_yet" initDistn_watterson-1Ne4.txt 0.05 9567 
```
For scenarios with variable true $s_\text{AA}$ (i.e. `s2`), with either type of simulation initial condition, 4 values for the dominance coefficient `h`, $h \in \{0, 0.5, 1, 5\}$, and 7 values for the selection coefficient `s2`, $s_\text{AA} \in \{0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.01, 0.02\}$ are analyzed.

For scenarios with variable true $s_\text{Aa}$ (i.e. `s1`), with either type of simulation initial condition, 2 values for the dominance coefficient `h`, $h \in \{5, \infty\}$ (in the $h=\infty$, it's equivalent to setting `s2=0` with the given `s1` value), and 7 values for the selection coefficient `s1`, $s_\text{AA} \in \{0.001, 0.002, 0.003, 0.004, 0.005, 0.01, 0.02\}$ are analyzed.

 For each combination of parameters, 500 replicates are generated, conditioning on each replicate's...
* population allele frequency never falls below $1/(4N_e)$ (equivalent to allele not being lost), and
* pooled minor sample allele frequency be greater than 0.05.

When computing likelihood surface, a symmetrical geometric grid (including zero) of $s$ values is adopted. The user can obtain the same grid with the following codes in an interactive python session:
```python
>>> import numpy as np
>>> from diplo_locus.utility import _get_geom_grid
>>> s_grid = _get_geom_grid(-0.75, 0.75, 50, Ne=1e4)
>>> s_grid
# array([-7.50000000e-01, -5.10969052e-01, -3.48119163e-01, -2.37170825e-01,
#        -1.61582602e-01, -1.10084945e-01, -7.50000000e-02, -5.10969052e-02,
#        -3.48119163e-02, -2.37170825e-02, -1.61582602e-02, -1.10084945e-02,
#        -7.50000000e-03, -5.10969052e-03, -3.48119163e-03, -2.37170825e-03,
#        -1.61582602e-03, -1.10084945e-03, -7.50000000e-04, -5.10969052e-04,
#        -3.48119163e-04, -2.37170825e-04, -1.61582602e-04, -1.10084945e-04,
#        -7.50000000e-05,  0.00000000e+00,  7.50000000e-05,  1.10084945e-04,
#         1.61582602e-04,  2.37170825e-04,  3.48119163e-04,  5.10969052e-04,
#         7.50000000e-04,  1.10084945e-03,  1.61582602e-03,  2.37170825e-03,
#         3.48119163e-03,  5.10969052e-03,  7.50000000e-03,  1.10084945e-02,
#         1.61582602e-02,  2.37170825e-02,  3.48119163e-02,  5.10969052e-02,
#         7.50000000e-02,  1.10084945e-01,  1.61582602e-01,  2.37170825e-01,
#         3.48119163e-01,  5.10969052e-01,  7.50000000e-01])
```

The parameter settings described above are hard-coded in the python scripts. To consider other alternatives, the user should modify the codes at their own discretion.

With the hard-coded parameters, it takes ~4.5 core hours (on HPC) to generate either of the `varS2` simulations, and ~2.5 hours for `varS1` simulations. The `.pkl` files generated from these pipelines take about ~9 Gb space.


### Plotting

#### Receiver Operating Characteristic (ROC) curves (Fig. S1--3)
To visualize the simulated data as shown in the manuscript, run
```shell
# plot ROC from the off-grid max-likelihood likelihood ratios
## Usage:
## python plot_DLsim_ROCs.py <generic_s2_pkl> <generic_s1_pkl> <sim_init> <fig_prefix>
### generic names of .pkl files that contain the likelihood computation results. Use the "\*" wildcard character(s) to match all files
s2_pkl_name="simulations/varS2_fromSeed87235_simInitSV_h\*LLinit\*_offGridMLE.pkl" 
s1_pkl_name="simulations/varS1_fromSeed9567_simInitSV_simH\*LLinit\*_offGridMLE.pkl" 
### prefix of the ROC figure
roc_prefix="simInitSV_s2-87235_s1-9527"
# now plot
python plot_DLsim_ROCs.py ${s2_pkl_name} ${s1_pkl_name} "Watterson neutral" ${roc_prefix} 
```
The commands above will generate multi-panel ROC curves for simulations initiated with the given distribution.

(placeholder for pics)

Likewise, ROCs for simulations with the fixed initial frequency of 0.01 can be generated using the same script:
```shell
# plot ROC from the off-grid max-likelihood likelihood ratios
### generic name of .pkl files that contain the likelihood computation results. Use the "\*" wildcard character(s) to match all files
s2_pkl_name="simulations/varS2_fromSeed10827_simInitFreq.01_h\*LLinit\*_offGridMLE.pkl" 
s1_pkl_name="simulations/varS1_fromSeed9567_simInitFreq.01_simH\*LLinit\*_offGridMLE.pkl" 
### prefix of the ROC figure
roc_prefix="simInitFreq.01_s2-10827_s1-9527"
# now plot
python plot_DLsim_ROCs.py ${s2_pkl_name} ${s1_pkl_name} "initial freq f=0.01" ${roc_prefix} 
```

(placeholder for pics)

([Back to TOC](#toc))

#### Boxplots for max-likelihood estimates (MLEs) (Fig. S4--6)

The plotting script `plot_DLsim_contrastBoxes.py` has the same usage as `plot_DLsiims_ROCs.py`.

In addition to two plots---one for either initial condition adopted in likelihood computation---it will also produce two tables with outlier replicates whose MLEs fall outside the range covered by the y-axes in the plot (that is, $[-0.01, 0.015]$).

```shell
# plot boxplots for MLEs
## simulations initiated with standing variation (SV):
s2_pkl_name="simulations/varS2_fromSeed87235_simInitSV_h\*LLinit\*_offGridMLE.pkl" 
s1_pkl_name="simulations/varS1_fromSeed9567_simInitSV_simH\*LLinit\*_offGridMLE.pkl" 
### prefix of the figure
prefix="simInitSV_s2-87235_s1-9527"
# now plot
python plot_DLsim_contrastBoxes.py ${s2_pkl_name} ${s1_pkl_name} "Watterson neutral" ${prefix} 

## simulations initiated with fixed frequency (Freq.01):
s2_pkl_name="simulations/varS2_fromSeed10827_simInitFreq.01_h\*LLinit\*_offGridMLE.pkl" 
### $s1_pkl_name stays the same
prefix="simInitFreq.01_s2-10827_s1-9527"
# now plot
python plot_DLsim_contrastBoxes.py ${s2_pkl_name} ${s1_pkl_name} "initial freq f=0.01" ${prefix} 
```

(placeholder for pics)


([Back to TOC](#toc))
## `DiploLocus` likelihood-based inference on `SLiM` simulations
<a id="Slim"> </a>





([Back to TOC](#toc))

## Examining genotype data around *LCT* gene in UK population from AADR database
<a id="LCT"> </a>

Data for this example is extracted from [Allen Acient DNA Resources databse v54.1](https://reich.hms.harvard.edu/allen-ancient-dna-resource-aadr-downloadable-genotypes-present-day-and-ancient-dna-data).
```shell
# Move to the corresponding folder (assume one starts at SLiM folders)
cd ../figS13-S14_LCT/

# Download dataset
## This step might take a while depending on the internet bandwidth
##  the tarball `v54.1_1240K_public.tar` takes ~5 Gb disk space.
mkdir v54.1_1240K
cd v54.1_1240K
wget https://reichdata.hms.harvard.edu/pub/datasets/amh_repo/curated_releases/V54/V54.1/SHARE/public.dir/v54.1_1240K_public.tar 
# Decompress the tarball
tar -xvf v54.1_1240K_public.tar
# check files
ls ./
#v54.1_1240K_public.anno  v54.1_1240K_public.ind  v54.1_1240K_public.tar
#v54.1_1240K_public.geno  v54.1_1240K_public.snp
## Return to the main folder
cd ../ 
```
The extracted files from this tarball are in EIGENSTRAT format and contain data for all genomes from published ancient DNA studies. To make subsequent file-handling easier, we first generate a subset of all UK samples for VCF conversion:
```shell
# create the folder to store vcfs
mkdir extracted/
# run script
python extractIDsFromAnnoUK.py
```
This step will generate `extracted/UK_v54.1_all_with_duplicates.table`, which will be used to make VCFs later.

To this end, we use a tool from the [gdc](https://github.com/mathii/gdc) package to convert EIGENSTRAT-formatted data to VCF format. This step could be done with the following commands:

```shell
# make sure the submodule folder isn't empty:
git submodule update --recursive

# convert eigenstrat to vcf:
## here we're only interested in chromosome 2, where LCT gene sits
python extractUKTimeSeries_byChr.py 2
```
The converted VCF file will be `extracted/UK_v54.1_all_with_duplicates_c2.vcf`.

Next, the script `filter_v54_vcfs.py` can be used to further narrow down the individuals and SNPs to be included in the `DiploLocus` input.
```shell
python filter_v54_vcfs.py
#usage: filter_v54_vcfs.py [-h] --vcf VCFFILE --anno ANNOFILE [--gen_time GEN_TIME]
#                          [--out OUTPREFIX] [--ID_col IDCOLNAME] [--time_ago_col TIMEAGOCOLNAME]
#                          [--inds INDS] [--country WHERE]
#                          [--SG | --no-SG | --Shotgun | --no-Shotgun]
#                          [--DG | --no-DG | --Diploid | --no-Diploid]
#                          [--1240K_only | --no-1240K_only | --1240k_only | --no-1240k_only]
#                          [--PASS {none,loose,strict}] [--Fam | --no-Fam] [--Modern | --no-Modern]
#                          [--fromYBP FROM_YEAR] [--toYBP TO_YEAR] [--force_hap FORCED_HAPS]
#                          [--force_dip FORCED_DIPS] [--snps SNPS] [--chr CHROM] [--pos POSRANGE]
#                          [--max_missing MAXMISSING | --min_observed MINOBS] [--minK MINK]
#                          [--minMAF MINMAF]
#filter_v54_vcfs.py: error: the following arguments are required: --vcf, --anno
```

One can also see the full help page with `-h` or `--help` flag.

<details>
<summary>Click here to read the full help page.</summary>

```shell
python filter_v54_vcfs.py -h
#usage: filter_v54_vcfs.py [-h] --vcf VCFFILE --anno ANNOFILE [--gen_time GEN_TIME]
#                          [--out OUTPREFIX] [--ID_col IDCOLNAME] [--time_ago_col TIMEAGOCOLNAME]
#                          [--inds INDS] [--country WHERE]
#                          [--SG | --no-SG | --Shotgun | --no-Shotgun]
#                          [--DG | --no-DG | --Diploid | --no-Diploid]
#                          [--1240K_only | --no-1240K_only | --1240k_only | --no-1240k_only]
#                          [--PASS {none,loose,strict}] [--Fam | --no-Fam] [--Modern | --no-Modern]
#                          [--fromYBP FROM_YEAR] [--toYBP TO_YEAR] [--force_hap FORCED_HAPS]
#                          [--force_dip FORCED_DIPS] [--snps SNPS] [--chr CHROM] [--pos POSRANGE]
#                          [--max_missing MAXMISSING | --min_observed MINOBS] [--minK MINK]
#                          [--minMAF MINMAF]
#
#optional arguments:
#  -h, --help            show this help message and exit
#
#Input:
#  --vcf VCFFILE         VCF files to be parsed.
#  --anno ANNOFILE       ".anno" file as presented in the v54.1 1240K database.
#  --gen_time GEN_TIME   Average generation time, in years, to be used for converting years to
#                        generations.Default value is 26.9 years (Wang et al. 2023;
#                        https://doi.org/10.1126/sciadv.abm7047)
#
#Output Options:
#  --out OUTPREFIX, -o OUTPREFIX, --outname OUTPREFIX
#                        path and prefix to output files.
#
#Sample-based Filters:
#  --ID_col IDCOLNAME    Name of the column in info table for ID names of the sample (as shown in
#                        VCF). Default is "ID".
#  --time_ago_col TIMEAGOCOLNAME
#                        Name of the column in info table for each sample's sampling times
#                        (backward; ascending from present to past). Default name is "MeanYBP", unit
#                        is years before present (:=1950).
#  --inds INDS           List (or file name) of sample IDs to include. Default is to include all
#                        samples shared between vcf & info files.
#  --country WHERE, --place WHERE, --locality WHERE, --region WHERE
#                        Key word in "Political Entity" or "Locality" column in the .anno file.
#  --SG, --no-SG, --Shotgun, --no-Shotgun
#                        Option to include ("True") or exclude ("False") shotgun-sequenced samples,
#                        as indicated in the "Data source" column. Default is to ignore such info.
#  --DG, --no-DG, --Diploid, --no-Diploid
#                        Option to include ("True") or exclude ("False") samples with diploid
#                        variant calls, as indicated in the "Data source" column. Default is to
#                        ignore such info.
#  --1240K_only, --no-1240K_only, --1240k_only, --no-1240k_only
#                        Option to only include samples sequenced via 1240K captures. That is, "Data
#                        source"_is exactly_ "1240[K|k]".
#  --PASS {none,loose,strict}
#                        Options on the requirement for quality filters."none": ignore such info;
#                        "loose" (default): samples whose "Assessment" column _contains_ "PASS";
#                        "strict": samples whose "Assessment" column _is exactly_ "PASS".
#  --Fam, --no-Fam       Option to include ("True") or exclude ("False") samples with known
#                        relatives recorded, as indicated in the "Family ID" column. Default is to
#                        ignore such info.
#  --Modern, --no-Modern
#                        Option to include ("True") or exclude ("False") samples with known
#                        relatives recorded, as indicated in the "Family ID" column. Default is to
#                        ignore such info.
#  --fromYBP FROM_YEAR, --from FROM_YEAR, --sinceYBP FROM_YEAR, --since FROM_YEAR
#                        The oldest time (including itself), in years before present, based on the
#                        sample's mean calibrated YBP, to consider. Default is to take all.
#  --toYBP TO_YEAR, --to TO_YEAR
#                        The latest time (excluding itself), in years before present, based on the
#                        sample's mean calibrated YBP, to consider. Default is to take all.
#
#VCF-based Filters:
#  --force_hap FORCED_HAPS
#                        Comma-separated string or a plain-text file that lists IDs (matching VCF
#                        column names) to be considered as haploids even though some may present
#                        diploid genotype calls. If "all", all samples will be considered as
#                        haploids, whose genotypes will be counted as matching haplotype (i.e. half
#                        the alleles). Force quit if any specified individual has heterozygote
#                        variant or any unmentioned diploid-formatted samples lack heterozygote
#                        variant calls. Default is "none".
#  --force_dip FORCED_DIPS
#                        Comma-separated string or a plain-text file that lists samples IDs
#                        (separated by comma) to be considered as diploids even though some may have
#                        loci with only haplotype calls. If "all", all samples will be considered as
#                        diploids, whose haploid GTs will be counted as matching homozygote diploid
#                        GTs (i.e. double-count the allele). Default is "none".
#  --snps SNPS           List (or file name) of SNPs to include. Default is to include all
#                        qualifying SNPs in the VCF.
#  --chr CHROM, -c CHROM
#                        Name of the chromosome to extract SNPs from. Must match values in the VCF.
#                        Default is to include all qualifying SNPs in the VCF.
#  --pos POSRANGE        Range of physical positions ("<from>-<to>") to consider on the proposed
#                        chromosome. Default is to include all qualifying SNPs in the VCF.
#  --max_missing MAXMISSING
#                        For each SNP, the max allowed proportion of samples to have missing data.
#  --min_observed MINOBS
#                        For each SNP, the min allowed proportion of samples to have valid
#                        observations (not missing).
#  --minK MINK           For each SNP, the minimal number of time points to have non-zero
#                        observations in.
#  --minMAF MINMAF       Minimum threshold (non-inclusive) for minor allele frequencies in the
#                        pooled samples. SNPs with sub-threshold frequencies will not be considered
#                        for analyses. Default is 0.
```

</details>

Using this script, we aim to select sampled individuals that are
* sampled from the United Kingdom,
* dated to be no older than 4500 years before present (1950 AC; YBP),
* sequenced from 1240K SNP-set capture,
* un-related to any other samples in the dataset 

To do this, the user can use command
```shell
## save filenames for convenience
extracted_vcf_all="extracted/UK_v54.1_all_with_duplicates_c2.vcf"
anno_file="v54.1_1240K/v54.1_1240K_public.anno"
### checkout column names
python filter_v54_vcfs.py --vcf $extracted_vcf_all --anno $anno_file --gen_time 28.1 \
                        --ID_col "GenID" --time_ago_col "MeanYBP" \
                        --country UK --no-SG --1240K_only --PASS strict --no-Fam \
                        --fromYBP 4500 --toYBP 0 --force_haps all
```

## Combination of previous figures
<a id="main_fig"> </a>

([Back to TOC](#toc))
