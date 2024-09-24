Multivariate Bayesian Analyses for Nursing Research
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

# 1 Copyright information

Copyright 2023, University of Pittsburgh. All Rights Reserved. License:
[GPL-2](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

# 2 Overview

This repository was created to support a presentation at the
International Society of Nurses in Genetics:

Heinsberg LW. Multivariate Bayesian Approaches for Analyzing Correlated
Phenotypes in Nursing Research. (Expert Lecturer Abstract, Podium).
Presented at the International Society of Nurses in Genetics, November
2023, Providence, Rhode Island.

which has been adapted for (hopeful) publication as a manuscript
entitled:

Heinsberg LW, Davis TS, Maher D, Bender CM, Conley YP, Weeks DE.
Multivariate Bayesian Analyses in Nursing Research: An Introductory
Guide. Submitted to Biological Research for Nurses.

The goal of this repository is to provide detailed code and a fully
synthetic example data set to guide nurse scientists and other
researchers in conducting multivariate Bayesian analyses to examine
associations between correlated phenotypes and single nucleotide
polymorphisms (SNPs, i.e., genetic variants). This guide will focus on
the application of two multivariate Bayesian software programs: (1)
`bnlearn` and (2) `mvBIMBAM`.

**Approach 1:** `bnlearn` - bnlearn is an R package for learning the
graphical structure of Bayesian networks,estimating their parameters,
and performing some useful inference as described in:

Scutari M, Howell P, Balding DJ, Mackay I. Multiple Quantitative Trait
Analysis Using Bayesian Networks. Genetics. 2014 Apr 11;198(1):129–137.
PMID: 25236454 PMCID: PMC4174925 DOI:
<https://doi.org/10.1534/genetics.114.165704>

Zhang JZ\*, Heinsberg LW\*, Krishnan M, Hawley NL, Major TJ, Carlson JC,
Harré Hindmarsh J, Watson H, Qasim M, Stamp LK, Dalbeth N, Murphy R, Sun
G, Cheng H, Naseri T, Reupena MS, Kershaw EE, Deka R, McGarvey ST,
Minster RL, Merriman TR, Weeks DE. Multivariate analysis of a missense
variant in CREBRF reveals associations with measures of adiposity in
people of Polynesian ancestries. Genet Epidemiol. 2023
Feb;47(1):105-118. doi: 10.1002/gepi.22508. Epub 2022 Nov 9. PMID:
36352773; PMCID: PMC9892232. \*First authors.

and at

[www.bnlearn.com](https://www.bnlearn.com/)

Additional details about the method and its interpretation can be found
in `01_mvNUR_bnlearn.Rmd`.

**Approach 2:** `mvBIMBAM` - mvBIMBAM implements a terminal-based
Bayesian approach for genetic association analysis of multiple related
phenotypes, as described in:

Shim H, Chasman DI, Smith JD, Mora S, Ridker PM, Nickerson DA, Krauss
RM, Stephens M. A multivariate genome-wide association analysis of 10
LDL subfractions, and their response to statin treatment, in 1868
Caucasians. PLoS One. 2015 Apr 21;10(4):e0120758. doi:
<https://doi.org/10.1371/journal.pone.0120758>. PMID: 25898129; PMCID:
PMC4405269.

Stephens M. A unified framework for association analysis with multiple
related phenotypes. PLoS One. 2013 Jul 5;8(7):e65245. doi:
<https://doi.org/10.1371/journal.pone.0065245>. Erratum in: PLoS One.
2019 Mar 19;14(3):e0213951. PMID: 23861737; PMCID: PMC3702528.

Additional details about the method and its interpretation can be found
in `02_mvNUR_mvBIMBAM.Rmd`.

Throughout all markdowns, we have flagged lines that will need to be
modified by the user if adapting this code for your own data using
`###CUSTOMIZE**` annotation.

Also note that while `bnlearn` is run directly in R, `mvBIMBAM` is a
terminal-based program with no point/click desktop app. As such, please
note that this example code, particularly the mvBIMBAM approach,
requires at least an introductory understanding of R and Unix. Some
great introductory R/Unix resources are listed at the end of this
document.

# 3 Installation

## 3.1 bnlearn

**bnlearn** is available on [CRAN](https://cran.r-project.org/) and can
be downloaded from its web page in the Packages section
([here](https://cran.r-project.org/web/packages/bnlearn/index.html)).

For example, you can install bnlearn from
[CRAN](https://cran.r-project.org/web/packages/bnlearn/index.html)
using:

    install.packages("bnlearn")

or

    install.packages("https://www.bnlearn.com/releases/bnlearn_latest.tar.gz", repos = NULL, type = "source")

Development snapshots, which include bugfixes that will be incorporated
in the CRAN release as well as new features, can be downloaded from the
links above or installed with a simple:

    install.packages("http://www.bnlearn.com/releases/bnlearn_latest.tar.gz")

The only suggested packages not hosted on
[CRAN](https://cran.r-project.org/) are **graph** and **Rgraphviz**,
which can be installed from
[BioConductor](https://www.bioconductor.org/):

    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install()
    BiocManager::install(c("graph", "Rgraphviz"))

following the instructions present on
[this](https://www.bioconductor.org/packages/release/bioc/html/Rgraphviz.html)
and
[this](https://www.bioconductor.org/packages/release/bioc/html/graph.html)
webpage. Please also note that the **gRain** package, while on CRAN,
depends on packages that are on Bioconductor both directly and through
the **gRbase** package, which depends on **RBGL**:

    BiocManager::install()
    BiocManager::install(c("graph", "Rgraphviz", "RBGL"))
    install.packages("gRain")

Check your R version if you are having any issues installing. As of
August 2024, bnlearn requires R version 4.4.0 or higher.

Additional installation and troubleshooting information can be found at
the bnlearn
[website](https://cran.r-project.org/web/packages/bnlearn/index.html).

## 3.2 mvBIMBAM

mvBIMBAM is designed to be installed from
[GitHub](https://github.com/heejungshim/mvBIMBAM). The binary executable
file for Linux is in the `https://github.com/heejungshim/mvBIMBAM/bin/`
directory while a zip file of the binary executable file, source code,
and example input files for Mac is in
`https://github.com/heejungshim/mvBIMBAM/forMAC` directory. NOTE: The
software executable is called ‘bimbam’, but we refer to the program as
‘mvBIMBAM’ throughout.

HOWEVER, when developing/testing this tutorial, several of us ran into
installation issues for Mac. Ultimately, we had to modify the program’s
makefile to enable successful installation. We have alerted the mvBIMABM
software maintainer of this issue, but have not yet seen it resolved. As
such, we provide a modified Mac source code download to make this
process easier for users of this tutorial (allowable under GNU General
Public License (GPLv3+)).

To install mvBIMBAM, please follow the instructions below, which were
written with the assumption that users have an introductory
understanding of navigating your computer using the terminal/Unix
commands. If you are new to using Unix commands in the terminal, please
visit this [link](https://macpaw.com/how-to/use-terminal-on-mac) for a
brief introduction.

### 3.2.1 Mac

1)  Download
    [`BIMBAM_multi_pheno_MAC_M1_modified.zip`](BIMBAM_multi_pheno_MAC_M1_modified.zip)
    or
    [`BIMBAM_multi_pheno_MAC_Intel_modified.zip`](BIMBAM_multi_pheno_MAC_Intel_modified.zip)
    from this repository. The former is for newer Apple silicon/M1 macs
    and the latter is for older macs. If you are not sure which Mac you
    have, you can click on the apple symbol in the upper left corner of
    your machine and choose the “About This Mac” option. If you have an
    Apple silicon machine, the chip will be listed as Apple M1 or M2. If
    you have the Apple Intel, you will see the processor listed as Intel
    Core i5, i7, or similar.

2)  Unzip the file by double clicking the zipped folder.

3)  Using the terminal, navigate to the newly unzipped mvBIMBAM folder

<!-- -->

    cd Desktop/BIMBAM_multi_pheno_MAC_M1_modified
    pwd 

or

    cd Desktop/BIMBAM_multi_pheno_MAC_Intel_modified
    pwd 

For example, on my machine, when I use the `pwd` command (print working
directory), my folder is located at:
/Users/username/Desktop/BIMBAM_multi_pheno_MAC

4)  Run `make clean` and `make all` terminal commands which will
    activate the software build automation tool to compile the C++
    source code into an executable `bimbam` program.

<!-- -->

    make clean
    make all

Note that there is, unfortunately, no obvious sign that install was
successful. Installation can be confirmed in step 5.

5)  Move the executable file to system’s PATH environment variable, such
    as /usr/local/bin, so the program can be called from any location
    (vs. only the folder where the executable file resides). To do this,
    navigate to the location of the Unix executable file via the
    terminal and then use the move (`mv`) command.

<!-- -->

    sudo mv bimbam /usr/local/bin/

Confirm that the install/move worked via this command:

    which bimbam

For example, on my machine, `bimbam` is located at
`/usr/local/bin/bimbam` now.

Now if you type

    bimbam

at the Unix prompt and hit return, it should print out this message:

     BIMBAM version 0.99a, visit http://stephenslab.uchicago.edu for possible update.
     Developed by Yongtao Guan ytguan.at.gmail.com, all rights to be lefted to GNU GPL.
     References: Guan and Stephens (2008), Servin and Stephens (2007), Scheet and Stephens (2006).
     Updated March 14 2009.

#### 3.2.1.1 Modification details

If you are interested in more details related to the installation issue,
and how we modified the makefile to correct this issue, please see
[`InstallationHack.md`](InstallationHack.md).

#### 3.2.1.2 Troubleshooting

Note that you need to have command line tools and GSL (GNU Scientific
Library) installed on your machine before you can install mvBIMBAM.
Please see [`Troubleshooting.md`](Troubleshooting.md) for tips on
getting GSL installed on your machine.

# 4 Synthetic data set overview

To facilitate hands-on learning, a synthetic data set was created for
use with this example code ([data](data)). The synthetic data set was
created from a study of women with breast cancer (BrCa). This data set
contains no real data, but mimics the correlation structure of the true
data set. The original data set contained data for 110 participants and
included data for 6 symptoms, 3 demographic variables, and 6 SNPs. We
created a synthetic version of this data set with fake data for 770
participants for use in this tutorial (for reference, see
[`SynthpopCode.md`](SynthpopCode.md).

Data set name: [`BrCa_synthetic.csv`](BrCa_synthetic.csv)

Variable descriptions:

Symptoms/Outcomes  
-`ID`: Participant ID  
-`EMO_tscore`: Anxiety (PROMIS Anxiety T Score)  
-`bdito`: Depression (Beck’s Depression Inventory)  
-`FAT_tscore`: Fatigue (PROMIS Fatigue Score)  
-`paohcif`: Cognitive Function (Patient’s Assessment of Own Functioning
Inventory - Cognitive Function)  
-`EPSscore`: Sleepiness (Epworth Daytime Sleepiness)  
-`worst_pain`: Pain (Self-reported numeric rating of “average pain”)

Demographics/Covariates  
-`age`: Age in years  
-`education`: Education in years  
-`race`: Race (Self-identified, 0=White, 1=Black)

SNPs/Predictors  
-rs \#s (SNPs)

``` r
# Read in a mapping of variable name to informative label
dict <- read.csv("data/BrCa_synthetic_SimpleDict.csv")
pander::pander(dict, "Data dictionary") 
```

|  Variable  |       Label        |                                                                                                                                                                                                                                        Description                                                                                                                                                                                                                                        |
|:----------:|:------------------:|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|
| Emo_tscore |      Anxiety       | PROMIS Emotional Distress Anxiety (Short Form 8a) - 8 items measuring emotional distress and anxiety (panic, fearfulness, worry, dread, tension,nervousness, restlessness, and somatic symptoms including racing heart and dizziness) in the last 7 days. This instrument generates T-scores, which are standard scores with a mean of 50 and standard deviation of 10 in the U.S. general population. Higher T-scores indicate worse distress and anxiety. T-scores range from 10 to 90. |
|   bdito    |     Depression     |                                                                                                                  The Beck Depression Inventory-II - 21-items measuring severity of depression (i.e., grief, loss, flat affect, withdrawal, mania) over a two-week period. Scores range from 0 to 63, where high scores indicate higher levels or more severe depression.                                                                                                                  |
| FAT_tscore |      Fatigue       |                                                         PROMIS Fatigue (Short Form 8a) - 8 items measuring the experience of fatigue and the interference of fatigue on daily activities over the past 7 days. This instrument generates T-scores, which are standard scores with a mean of 50 and standard deviation of 10 in the U.S. general population. Higher T-scores indicate worse fatigue. T-scores range from 10 to 90.                                                         |
|  paohcif   | Cognitive function |                                                              Patient assessment of own cognitive functioning - 33 items measuring self-reported ‘recent’ cognitive function consisting of 5 dimensions: memory, language and communication, use of hands, sensory-perceptual, higher level cognitive and intellectual functioning. Higher scores indicate worse perceived cognitive functioning. Scores range from 0 to 155.                                                              |
|  EPSscore  |     Sleepiness     |                                                                                                                                                         Epworth Sleepiness Scale - 8 item scale measuring self-reported ‘recent’ daytime sleepiness. Higher score indicate greater daytime sleepiness. Scores range from 0 to 24.                                                                                                                                                         |
|    pain    |        Pain        |                                             Brief Pain Inventory (Short Form) - 9 items a widely used clinical tool for assessing pain (worst pain in the last 24 hours $$0-10$$, least pain in the last 24 hours $$0-10$$, pain on average $$0-10$$, paing at the time of the interview $$0-10$$, pain relief from treatment(s) in the last 24 hours $$0%-100%$$). Higher scores indicate higher levels of pain, each item ranges from 0-10.                                             |
|   rs4880   |     Variant 1      |                                                                                                                                                                                                                          SOD2 gene, hg19 postion chr6: 160113872                                                                                                                                                                                                                          |
| rs5746136  |     Variant 2      |                                                                                                                                                                                                                         SOD2 gene, hg19 postiion chr6: 160103084                                                                                                                                                                                                                          |
| rs1041740  |     Variant 3      |                                                                                                                                                                                                                         SOD1 gene, hg19 postiion chr21: 33040162                                                                                                                                                                                                                          |
| rs10432782 |     Variant 4      |                                                                                                                                                                                                                         SOD1 gene, hg19 postiion chr21: 33036391                                                                                                                                                                                                                          |
| rs4135225  |     Variant 5      |                                                                                                                                                                                                                          TXN gene, hg19 postiion chr9: 113006691                                                                                                                                                                                                                          |
| rs7522705  |     Variant 6      |                                                                                                                                                                                                                         PRDX1 gene, hg19 postiion chr1: 45992300                                                                                                                                                                                                                          |
|    race    |        Race        |                                                                                                                                                                                                                          Self-identified race, 0=White, 1=Black                                                                                                                                                                                                                           |
|    age     |        Age         |                                                                                                                                                                                                                                       Age in years                                                                                                                                                                                                                                        |
| education  |     Education      |                                                                                                                                                                                                                             Self-reported years of education                                                                                                                                                                                                                              |

Data dictionary (continued below)

|          Type          |                                                                                                     Citation                                                                                                      |
|:----------------------:|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|
|        Numeric         | Pilkonis, P.A. et al. (2011) Item banks for measureing emotional distress from the patient-reported outcomes measurement information system (PROMIS): Depression, Anxiety, and Anger, Assesssment, 18(3), 263-283 |
|        Numeric         |                                                         Beck, A.T. et al. (1996) Beck Depression Inventory-II. San Antonio: The Psychological Corporation                                                         |
|        Numeric         |                              Cell, D. et al. (2016) PROMIS fatigue item bank had clinical validity across diverse chronic conditions. Journal of Clinical Epidemiology, 73, 128-134                               |
|        Numeric         |                        Chelune, G. J. et al. (1986) Neuropsychological and personality correlates of patients’ complaints of disability. Advances in clinical Neuropsychology, 1986:95-126                        |
|        Numeric         |                                                  Johns MW. A new method for measuring daytime sleepiness: The Epworth Sleepiness Scale. Sleep 1991; 14(6):540-5                                                   |
|        Numeric         |                                   Cleeland, C.S. et al. (2009). Pain assessment: Global use of the Brief Pain Inventory. Annals, Academy of Medicine, Sigapore, 23(2), 129-138.                                   |
|        Numeric         |                                                                                                                                                                                                                   |
|        Numeric         |                                                                                                                                                                                                                   |
|        Numeric         |                                                                                                                                                                                                                   |
|        Numeric         |                                                                                                                                                                                                                   |
|        Numeric         |                                                                                                                                                                                                                   |
|        Numeric         |                                                                                                                                                                                                                   |
| Numeric, encoded value |                                                                                                                                                                                                                   |
|        Numeric         |                                                                                                                                                                                                                   |
|        Numeric         |                                                                                                                                                                                                                   |

# 5 Examples

To prepare the data before running the examples, see
[`00_mvNUR_DataPrep.Rmd`](00_mvNUR_DataPrep.Rmd) (web version:
[`00_mvNUR_DataPrep.md`](00_mvNUR_DataPrep.md)). After running that
markdown, examples of each method can be found in
[`01_mvNUR_bnlearn.Rmd`](01_mvNUR_bnlearn.Rmd) (web version:
[`01_mvNUR_bnlearn.md`](01_mvNUR_bnlearn.md)) and
[`02_mvNUR_mvBIMBAM.Rmd`](02_mvNUR_mvBIMBAM.Rmd) (web version:
[`02_mvNUR_mvBIMBAM.md`](02_mvNUR_mvBIMBAM.md)).

# 6 R/Unix basics

- <https://danieleweeks.github.io/HuGen2071/preparation.html>
- <https://github.com/ajwills72/rminr>
- <https://www.andywills.info/rminr/#beginners>
- <https://www.youtube.com/watch?v=IrDUcdpPmdI>

# 7 Power calculation resources

- <https://www.andywills.info/rminr/power-bayesian.html>
- <https://solomonkurz.netlify.app/blog/bayesian-power-analysis-part-i/>
- <https://cran.r-project.org/web/packages/BayesianPower/vignettes/bayesianpower.html>

# 8 Other Bayesian software

- [BayesSUR](https://github.com/mbant/BayesSUR)

- [tan](https://pgmpy.org/examples/Structure%20Learning%20with%20TAN.html)

- [jags](https://mcmc-jags.sourceforge.io/)

- [brms](https://paul-buerkner.github.io/brms/)

- [blavaan](https://ecmerkle.github.io/blavaan/)

- [bayesian_first_aid](https://github.com/rasmusab/bayesian_first_aid)

# 9 Contact information

If you have any questions or comments, please feel free to contact me!

Lacey W. Heinsberg, PhD, RN: <law145@pitt.edu>

# 10 Acknowledgments

I’d like to express my gratitude to the following for their support and
contributions to this repository:

- Support from the National Institutes of Health under award numbers
  K99HD107030 and R00HD107030 made this project possible, and for that,
  I’m truly grateful.
- Special thanks to Dr. Tara Davis, Mr. Dylan Maher, and Dr. Daniel
  Weeks for being “test users” and providing their invaluable feedback
  on this guide.
- My deepest gratitude to Dr. Daniel Weeks. Without his guidance and
  inspiration, I would never have ventured into the world of Bayesian
  statistics, or pursued many other endeavors I once thought were beyond
  my capabilities.

*I have a Bayesian inference joke but the first three people I told it
to didn’t laugh and now I’m not so sure it’s funny.* - @JSEllenberg
