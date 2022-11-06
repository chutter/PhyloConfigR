# PhyloConfigR: R Package For phylogenomic configuration, alignment formatting, alignment summary statistics, and dataset filtering

With this R package you can prepare your alignments for various types of phylogenomic analyses: 
1) Summarize alignment statistics (length, samples, parsimony informative sites, missing data, GC content) 
2) Renaming (alignment names, taxa in an alignment)
3) Alignment trimming
4) Filter alignments and create sub datasets from alignment stats
5) Filter gene trees based on alignments
6) Collapse nodes with specific support into polytomies
7) Remove taxa from alignment

In addition, there are numerous functions to aid with specific phylogenetic software
1) Astral-III
  - Create input files
  - Run and import ASTRAL-III results into R
  - Batch analyses multiple datasets
  - Plot Astral-III pie charts (proportion of gene trees that support that relationship) at nodes
2) IQTREE2
  - Plot concordance factors on tree
  - Gene jackknife genomic dataset
3) BPP
  - Set up alignments for immediate input
  - Generate control file
4) PhyloNet
  - Set up alignments for input


## Citation

Publication is in review. Hutter and Duellman, in review. 

For now, you can cite the R package by linking to this GitHub if you use it. 


## 1) Installation

To use all the functions, the following programs are needed if you wish to use them:
  1) ASTRAL-III is available on GitHub here: https://github.com/smirarab/ASTRAL
  2) IQTREE 2, specifically version 2 or higher: http://www.iqtree.org
  3) BPP, latest version 
  
Instructions for installation and testing ASTRAL-III, IQTREE 2, and BPP are included in the respective documentation.

PhyloConfigR depends on three other R packages:
  1) ape (>=3.0)
  2) data.table (>=1.12)
  3) stringr (>= 1.4)
  
To install PhyloConfigR, you can use the R package devtools. Here are step-by-step instructions for installation:

1) Install devtools by typing in your R console: install.packages("devtools", dependencies = T)

2) Install AstralPlane by typing in your R console: devtools::install_github("chutter/PhyloConfigR")

3) Devtools will ask you to install the package dependecies (ape and stringr), select "Yes". If devtools asks you to update packages, you may choose to do so. I would recommend not to install packages from source if devtools asks you. Ape is problemic from source and I could not get it to install on my machine. If devtools is still giving you trouble, you can install the dependencies with "install.packages(c("ape", "stringr"))". Then rerun Step 2 and skip package updates. 

4) Install PhyloConfigR by typing in your R console: devtools::install_github("chutter/PhyloConfigR"), and repeat Step 3. 

5) Devtools should finish and say the packages loaded properly. Load the packages with library(PhyloConfigR). 

And installation should be complete. 


## 2) Setting up R environment

I have included an R script in the main repository with some examples. It is also described here in detail. 

1) first install and load the R package. Its a good idea to install new (or check) every time as this package is being updated frequently. Functions may also be modified and stop working, so check back here for updated tutorial instructions. 

```r
devtools::install_github("chutter/PhyloConfigR")
library(PhyloConfigR)

```

2) You will want a character variable that includes your full path to the astral and iqtree jar files. NOTE: if you move the astral jar file, you will need to move the lib/ directory along with it, as astral depends on it. 


```r
astral.path = "/usr/local/bin/Astral-5-14/astral.5.14.2.jar"
iqtree.path = "/usr/local/bin/IQTREE/bin/iqtree2"
bpp.path = "/usr/local/bin/BPP/bin/bpp"

```

3)Setup your working directory and create if necessary

```r
work.dir = "/Test_FilterZone"
dir.create(work.dir)
setwd(work.dir)
```


