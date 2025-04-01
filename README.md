
## Figures for _"`diplo-locus`: A lightweight toolkit for inference and simulation of time-series genetic data under general diploid selection"_

This repository contains repoducible scripts for generating the simulations and figures for the manuscript and supplemental material.

To cite our application, use the reference:

>   Cheng, X.† & Steinrücken, M.† (2023) `diplo-locus`: A lightweight toolkit for inference and simulation of time-series genetic data under general diploid selection. _biorxiv_: https://doi.org/10.1101%2F2023.10.12.562101.

------------------
### Set up

In order for scripts in this repository to run, the Python version must be at least 3.8. The user must also have the [`diplo_locus`](https://github.com/steinrue/diplo_locus/) package installed. Both the API and CLI will be used.

To install __diplo-locus__ so that the module `diplo_locus` can be loaded in python, run
```shell
pip install "git+https://github.com/steinrue/diplo_locus@v1.2.0"
```
Download this repository with
```shell
# download
git clone https://github.com/steinrue/diplo_locus_manuscript_figs
# enter the directory
cd diplo_locus_manuscript_figs
```

All scripts in the instructions below are desgined (and tested) for a Unix-style command-line interface.

<a id="toc"> </a>

### Figures and tables

The figures in the supplemental material are listed first, since the figures and table in the main text depend on some of these:

#### Supplemental Material:
* [Figure S1 - S7: Data simulation and likelihood-based inference](supp_figS1-S7_DLsims)
* [Figure S8: Estimate onset of selection](supp_figS8_onset)
* [Figure S9 - S13: `DiploLocus` likelihood-based inference on `SLiM` simulations](supp_figS9-S13_SLiM)
* [Figure S14 - S15: Benchmarking `DiploLocus`](supp_figS14-S15_comparison)
* [Figure S16: Likelihood surfaces for ASIP in horses](supp_figS16_ASIP)
* [Figure S17 - S18: Examining genotype data around *LCT* gene in GB population from AADR database](supp_figS17-S18_LCT)
* [Figure S19 - S21: Analyzing Evolve & Resequence data for Drosophila simulans](supp_figS19-S21_dsim)

#### Main text:
* [Figure 1: ROC curves an accuracy of MLEs](main_fig1)
* [Figure 2: Manhattan plots for human and D. sim data; Likelihood surfaces for rs4988235 in humans and ASIP in horses.](main_fig2)
* [Table 1: Runtimes for benchmarck of `DiploLocus`](main_tab1)
