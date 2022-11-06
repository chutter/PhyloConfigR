## PhyloConfigR: R Package For phylogenomic configuration, alignment formatting, alignment summary statistics, and dataset filtering

In extreme cases of incomplete lineage sorting (ILS), it is possible that the most common gene tree topology will not match the true species tree, a phenomenon that has been termed “the anomaly zone”. For species trees in the anomaly zone, concatenation methods can provide strong support for the most common anomalous topology (i.e. anomalous gene trees: “AGTs”) while species tree methods can recover the correct species tree as ILS is into account. However, erroneous gene trees ("EGTs") has been shown to lead to erroneous species tree topologies when gene tree estimation error is high. EGTs resulting from non-biological properties of alignments (e.g. missing data, informative sites, alignment length) produces discordant EGTs from the true species tree, and filtering based on the informativeness of the alignments can lead to more robust species tree estimation (Hutter and Duellman, in review). 

Main features of the package:
  1) Testing for anomaly zones
  2) Calculates comprehensive alignment statistics for a single or folder of alignments 
  3) Creates filtered datasets from alignment statistics (alignment length, informativeness, etc). Generates filtered datasets for:
      - Alignments (concatenated or a folder)
      - Gene trees (using the alignment statistics)
  4) Assess your filtration results through concordance factors (Minh et al. 2020) directly from R (iqtree2 program required)
  5) Plot your filtration concordance factor and anomaly zone results

This package is still in beta testing phase, and more features and expanded functionality will be added in the future. If you find any issues with something not working, or you would like features to be added, go to issues in the top menu bar and submit them. 


# Citation

Publication is in review. Hutter and Duellman, in review. 

For now, you can cite the R package by linking to this GitHub if you use it. 


# Contents

1) Installation
2) Setting up R environment
3) Detecting the anomaly zone
4) Alignment and genetree dataset filtration 
5) Testing effectiveness of filtering on anomaly zone
6) Plot filtering anomaly zone results


## 1) Installation

For the full analysis pipeline, the following programs are needed:
  1) ASTRAL-III is available on GitHub here: https://github.com/smirarab/ASTRAL
  2) IQTREE 2, specifically version 2 or higher: http://www.iqtree.org
  
Instructions for installation and testing ASTRAL-III and IQTREE 2 are included in the respective documentation.

FilterZone depends on two other R packages:
  1) AstralPlane (>=0.1): https://github.com/chutter/AstralPlane
  2) data.table (>=1.12)

The package has two additional R package dependencies, which are treated as imports (i.e. you need them installed, but library(ape) and library(stringr) not needed: 
  ape (>= 5.0)
  stringr (>= 1.4)
  
To install FilterZone, you can use the R package devtools. Here are step-by-step instructions for installation:

1) Install devtools by typing in your R console: install.packages("devtools", dependencies = T)

2) Install AstralPlane by typing in your R console: devtools::install_github("chutter/AstralPlane")

3) Devtools will ask you to install the package dependecies (ape and stringr), select "Yes". If devtools asks you to update packages, you may choose to do so. I would recommend not to install packages from source if devtools asks you. Ape is problemic from source and I could not get it to install on my machine. If devtools is still giving you trouble, you can install the dependencies with "install.packages(c("ape", "stringr"))". Then rerun Step 2 and skip package updates. 

4) Install FilterZone by typing in your R console: devtools::install_github("chutter/FilterZone"), and repeat Step 3. 

5) Devtools should finish and say the packages loaded properly. Load the packages with library(AstralPlane), library(FilterZone), and library(data.table). 

And installation should be complete. 



## 2) Setting up R environment

I have included an R script in the main repository with some examples. It is also described here in detail. 

1) first install and load the R package. Its a good idea to install new (or check) every time as this package is being updated frequently. Functions may also be modified and stop working, so check back here for updated tutorial instructions. 

```r
devtools::install_github("chutter/AstralPlane")
library(AstralPlane)

devtools::install_github("chutter/FilterZone")
library(FilterZone)
library(data.table)

```

2) You will want a character variable that includes your full path to the astral and iqtree jar files. NOTE: if you move the astral jar file, you will need to move the lib/ directory along with it, as astral depends on it. 


```r
astral.path = "/usr/local/bin/Astral-5-14/astral.5.14.2.jar"
iqtree.path = "/usr/local/bin/IQTREE/bin/iqtree2"
```

3)Setup your working directory and create if necessary

```r
work.dir = "/Test_FilterZone"
dir.create(work.dir)
setwd(work.dir)
```


## 3) Detecting the anomaly zone

1) To detect the anomaly zone, you need a species tree estimated with coalescent branch lengths, where ASTRAL-III will provide this for you. As input into ASTRAL-III, you will need gene trees estimated separately for each alignment marker in your dataset. The R package AstralPlane provides some R functions that will streamline your data analysis pipeline:
  a. alignment and gene concatenation
  b. Within gene tree filtering (collapsing nodes with low support, taxa removal)
  c. Prepare gene trees for input into ASTRAL-III
  d. A wrapper to run ASTRAL-III and import results into R
  e. Plotting and results viewing

Instructions can be found here: https://github.com/chutter/AstralPlane

Once you have a species tree through AstralPlane or on your own, you can import this species tree into R. 

2) Create a set of character variables with the path to your tree file. Also indicate your outgroups for rooting the tree. Finally, the save.name is the desired output save name. 

```r
tree.file.path = "/Trees/test-tree.tre"
outgroups = c("Species_A", "Species_B")
save.name = "test-dataset"
```

3) Next, read the tree file into R, where the read.tree function from ape works to read in trees from ASTRAL-III. Alternatively, the file path to the tree file can be input directly into the "tree" parameter in the anomalyZone function and the function will read the tree. 


```r
test.tree = ape::read.tree(tree.file.path)
anom.data = anomalyZone(tree = test.tree,
                        outgroups = outgroup.taxa)
```

Alternatively, to read the file in directly from a file path:

```r
anom.data = anomalyZone(tree = tree.file.path,
                        outgroups = outgroup.taxa)
```

Parameter explanations: 

```r
tree = tree file from ASTRAL-III read into R as a phylo object or a file path to this tree. 
outgroups = your outgroup taxa for rooting the tree
```

4) When you have collected the data for the anomaly zone across the tree, you can view the data.frame that contains the nodes and branches where the anomaly zone was detected. Additionally, you can plot the results on the phylogenetic tree: 

```r
#Estimated run time: 1 second
plot.anomalyZone(tree = uce.tree,
                 data = anom.data,
                 outgroups = outgroup.taxa,
                 save.file = NULL,
                 tip.label.size = 0.5,
                 node.label.size = 1,
                 edge.width = 3)
```

Parameter explanations: 

```r
tree = tree file from ASTRAL-III read into R as a phylo object or a file path to a tree file
data = the output data.frame from the anomalyZone function
outgroups = your outgroup taxa for rooting the tree
save.file = NULL or blank to not save a file; otherwise file name to save PDF
tip.label.size = size of the tip labels, passed to the cex function of ape::plot
node.label.size = size of the node labels, passed to the cex function of ape::nodelabels
```

![](/pics/az-example-plot.svg)


## 4) Alignment and genetree dataset filtration 

The anomaly zone occurs when there are extreme cases of ILS and the most common gene tree topology does not match the true species tree. Species tree methods are designed to take into account ILS, however, they were not designed to take gene tree estimation error into account. The publication that introduces this R package shows that filtering of gene trees based on features of the alignment can result in better supported species trees and less gene tree species discordance. Importantly, filtering can aid in removing anomaly zones from species trees. 

1) To begin, you will first need a folder of alignments in phylip format and a folder of gene trees from IQTREE (other programs will probably work; if not, let me know and I can add them in). Create your working directory first (or use an existing directory). tree.files and align.files link to the gene tree files and alignments that estimated them. The names must match between the genes and alignments (except for the file extension). 


```r
work.dir = "WorkingDirectory"
tree.files = "WorkingDirectory/gene-trees"
align.files = "WorkingDirectory/alignments"
```

2) Next, you will want to select your filters to use, adjusted to the features of your dataset. Here is an example: 

```r
filter.length = c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
                  1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000,
                  2100, 2200, 2300, 2400, 2500) #number of base pairs
filter.sample = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1) #proportion samples
filter.prop.pis = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1) #proportion pis
filter.count.pis = c(10, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500) #count of pis
```

3) To obtain a table of statistics for each alignment, run the summarizeAlignments function. The inputs are the alignment directory path and the file export name. This function is general and can be used for other purposes. 


```r
#Estimated run time: 10 minutes
align.summary = summarizeAlignments(alignment.path = align.files,
                                    file.export = "alignment_stats",
                                    alignment.type = "phylip",
                                    dataset.name = "exons")
```

Parameter explanations: 

```r
alignment.path: path to a folder of multiple sequence alignments in phylip format
file.export: if a name is provided, the table is saved to file
overwrite: if TRUE overwrites file if it exists; FALSE the dataset is skipped
dataset.name: a unique name for your dataset. i.e. exons, introns, UCEs
alignment.type: select the format of the alignment. Phylip is avaialble, will be expanded in the future.
```

The summary table created by the function has the following columns: 

```
alignment name = the file name of the alignment
number_samples = number of samples or taxa in the alignment (number of rows) 
proportion_samples = the proportion of samples in the alignment (number_samples / max samples)
alignment_length = total number of base pairs long for the alignment
count_pis = number of parsimony informative sites in the alignment
proportion_pis = proportion of sites that are informative (count_pis / alignment_length)
count_missing_bp = total number of bases missing from the alignment matrix
proportion_missing_bp = proportion of bases missing from the alignment matrix (count_missing_bp / total bp)
```


4) Now that alignment statistics have been calculated, the filterSummary function can be used to obtain a quick summary of the datasets that will be generated using your selected filters (from above) and the alignment statistics. 

```r
#Estimated run time: 10 seconds
filt.summary = filterSummary(alignment.data = align.summary,
                             alignment.folder = align.dir,
                             dataset.name  = "exons",
                             file.out = "filter_summary",
                             length.filters = filter.length,
                             sample.filters = filter.sample,
                             prop.pis.filters = filter.prop.pis,
                             count.pis.filters = filter.count.pis)
```


Parameter explanations: 

```r
alignment.data: Alignment summary stats calculcated from summarizeAlignments
alignment.folder: The alignment folder from which the stats were calculated from in alignment.data
dataset.name: The name of your dataset, where all filtered datasets will be placed in this folder
file.out: if you wish to save to file, provide a file name for the summary
overwrite: whether to overwrite an existing dataset
length.filters: Your selected length filters as a vector of values for the alignment length
sample.filters: Your selected sampling fraction filters as a vector of values between 0-1
prop.pis.filters: Your selected parsimony informatives sites filter as a vector of values between 0-1
count.pis.filters: Your selected parsimony informatives sites filter as a vector of base pair counts
```

5) After alignment and filtered datasets summary statistics have been calculated, these two can be combined together to create concatenated alignments and parittion files for each filtered dataset. These concatenated alignments can be used in concatenation-based phylogenetics software (e.g. IQTREE, RAXML) to test out the impact of filtering on concatenation analyses. 

```r
#Estimated run time: 10 minutes
filterAlignments(filter.summary = filt.summary,
                 alignment.data = align.summary,
                 alignment.folder = align.dir,
                 format = "concatenated",
                 min.alignments = 5,
                 overwrite = TRUE)
```


Parameter explanations: 

```r
filter.summary: summary data file from filterSummary
alignment.data: summary data file from alignmentSummary
alignment.folder: folder of alignments to be filtered
format: save format for alignments
min.alignments: minimum number of alignments for a filtered set of data
overwrite: if TRUE overwrites file if it exists; FALSE the dataset is skipped
```

6) Prior to filtered summary species tree analysis in ASTRAL-III, the next function creates gene tree datasets for each filtered dataset, prepared for input in ASTRAL-III. Each filtered gene tree dataset is concatenated together (or placed in a folder with format = "folder") and saved to a folder called "filtered-genetrees-concatenated" for the concatenated gene trees or "filtered-genetrees-folders" for directories of gene trees for each filtered dataset. 

```r
#Estimated run time: 10 minutes
filterGeneTrees(filter.summary = filt.summary,
                alignment.data = align.summary,
                genetree.folder = tree.dir,
                format = "concatenated",
                overwrite = TRUE,
                min.trees = 5,
                min.n.samples = 4,
                make.polytomy = TRUE,
                polytomy.limit = 10,
                remove.node.labels = FALSE)
```

Parameter explanations: 

```r
filter.summary: summary data file from filterSummary
alignment.data: summary data file from alignmentSummary
genetree.folder: your target folder of gene trees that correspond to the alignments being filtered
format: save format for genetrees
overwrite: if TRUE overwrites file if it exists; FALSE the dataset is skipped
taxa.remove: species that you would like removed from each gene tree
min.trees: mimimum number of trees to keep filtered set
min.n.samples: the minimum number of samples to keep a gene tree
min.sample.prop: the minimum proportion of samples to keep a gene tree
make.polytomy: collapses polytomies in the gene trees
polytomy.limit: the value at which to collapse a node into a polytomy
remove.node.labels: strips trees of node labels if downstream analyses give you trouble (not recommended)
```

7) To run the new filtered sets of gene trees through ASTRAL-III efficiently, the AstralPlane package function astralRunner can be used on a folder of concateanted gene trees. the concatenated gene tree dataset can be provided to the astralRunner function, which runs ASTRAL-III for each concatenated set of gene trees from each filtered dataset in the "filtered-genetrees-concatenated" folder. Each summary species tree is saved in the output.dir. 

```r
#Estimated run time: hours, depending on how many filters selected
AstralPlane::astralRunner(concat.genetree.folder = "filtered-genetrees-concatenated",
                          output.dir = "filtered-astral",
                          overwrite = TRUE,
                          astral.path = astral.path,
                          astral.t = 2,
                          quiet = TRUE,
                          multi.thread = TRUE,
                          memory = "8g")               
```

Parameter explanations: 

```r
concat.genetree.folder: a folder of genetree files that are concatenated.
output.dir: the output directory name for the astral file
overwrite: overwrite = TRUE to overwrite existing files
astral.path: the absolute file path to the ASTRAL-III jar file
astral.t: the ASTRAL-III "t" parameter for different annotations, t = 2 is all annotation
quiet: TRUE hides the screen output from astral
multi.thread: TRUE to use Astral-MP multithreading 
memory: memory value to be passed to java. Should be in "Xg" format, X = an integer
```

8) Finally, concordance factors from IQTREE 2 can be computed on the species tree using the filtered gene trees, filtered alignments, and their corresponding species trees. The concordance factors for sites and genes will computed for every node in every filtered replicate, and these results can be plotted in the next section to display how filtering affects concordance factors across the entire tree or a focal node. 

```r
AstralPlane::concordanceRunner(alignment.dir = "filtered-alignments-concatenated",
                               species.tree.dir = "filtered-Astral",
                               genetree.dir = "filtered-genetrees-concatenated",
                               output.dir = "concordance-factors",
                               iqtree.path = "iqtree2",
                               overwrite = TRUE,
                               quiet = TRUE,
                               threads  = 6)         
```

Parameter explanations: 

```r
alignment.dir: The alignment folder from which the stats were calculated from in alignment.data
species.tree.dir = output directory from previous step that contains filtered ASTRAL-III species trees
genetree.dir: a folder of genetree files that are concatenated
output.dir: the output directory name to save the concordance factors results from filtration replicates
iqtree.path: path to IQTREE 2 if R cannot find it
overwrite: overwrite = TRUE to overwrite existing files
quiet: TRUE hides the screen output from astral
threads: the number of threads to used, passed to IQTREE
```


## 5) Testing effectiveness of filtering on anomaly zone

1) Now that alignment and filtration statistics have been calculated and filtered ASTRAL-III trees have been estimated, this collection of data can be analyzed together. First the necessary directory paths can be put into character vectors:

```r
work.dir = "WorkingDirectory"
astral.dir = "filtered-astral"
setwd(work.dir)
```

The working directory "work.dir" is the main project folder. The "astral.dir" is the directory of filtered astral datasets saved in the "output.dir" in Step 6 above ("filtered-astral"). 


2) The goal of the next step will be to target and obtain results for specific nodes or clades. To pull results from specific clades, first create an R list object. Include all taxa in the tree from each clade in each list object position, and name the list to correspond to the desired clade names. An example is show below: 

```r
#outgroups
outgroup.taxa = c("Outgroup_genus_1", "Outgroup_genus_2")

#Set up a list of your clades of interest
taxa.set = list()
taxa.set[[1]] = c("Genus_species_1", "Genus_species_2", "Genus_species_3")
taxa.set[[2]] = c("Genus_species_1", "Genus_species_2")
taxa.set[[3]] = c("Genus_species_1", "Genus_species_2", "Genus_species_3", "Genus_species_4")
names(taxa.set) = c("node1", "node2", "node3")
```

3) Load in the alignment and filtered summary data calculated in the previous section using "read.csv". 


```r
align.summary = read.csv("alignment_summary.csv")
filt.summary = read.csv("filter_summary.csv")
```

4) Once the input data is ready, run the filterAnomalies function to collect anomaly and erroneous zone data from all the filtered datasets. This data is calculated across all filtration replicates across all branches and nodes in the tree. 

```r
#Estimated run time: 1 minute
anomaly.data = filterAnomalies(astral.directory = astral.dir,
                               outgroups = outgroup.taxa,
                               filter.data = filt.summary)
```

Parameter explanations: 

```r
astral.directory: directory of filtered astral results
outgroups: outgroups to root your tree
filter.data: your master filtered dataset summary stats

```

5) Next, to obtain concordance factors data from the filtered datasets, run the filterConcordance function. The resulting table will contain the site and gene concordance factors calculated for each node across all the filtration replicates. These data.frames from filterAnomalies and filterConcordance can be used for other analyses. 

```r
#Estimated run time: 1 minute
concord.data = filterConcordance(input.dir = "concordance-factors",
                                 clade.list = taxa.set,
                                 outgroups  = outgroup.taxa)
```

Parameter explanations: 

```r
input.dir: directory of concordance factor data generated from the filtered datasets
clade.list: a named list of clades of interest to test for concordance factors
outgroups: outgroups to root the tree
```

6) Finally, with the "bestFilterTrees" function, the optimal single or multiple best trees from the filtration replications can be saved separately to file for use in publications and for other analyses. 

```r
#### Pull out best tree
bestFilterTrees(anomaly.zone.data = anomaly.data,
                concordance.factors.data = concord.data,
                output.dir = "best-trees",
                min.trees = 20,
                fewest.anomaly.zones = TRUE,
                highest.gene.cf = TRUE,
                highest.post.prob = TRUE,
                all.datasets = TRUE,
                top.best = 1)
```

Parameter explanations: 

```r
anomaly.zone.data: table output from the filterAnomalies function
concordance.factors.data:  table output from the filterConcordance function
output.dir: the name of the output directory to save the plots if TRUE above
single.best: saves a pdf of the single best tree
top.best: saves the top five best trees
all.datasets: TRUE for top five for each datasets FALSE for top five across all datasets
fewest.anomaly.zones: TRUE to return the trees with the fewest anomaly zones
highest.post.prob: TRUE to return the trees with the highest mean posterior probability
highest.gene.cf: TRUE to return the trees with the highest mean gene concordance factors
min.trees: minimum number of trees to keep a filtration replicate. Default: 10
```


This R package is meant to facilitate ASTRAL-III analyses and provide easy R plotting. The package helps prepare analyses from a folder of gene trees, runs astral from R, and creates a new S4 object type "AstralPlane" for easily analyzing the output from ASTRAL-III. The packageprovides several different types of plots, from pie charts on phylogenetic trees representing the quartet frequencies to plotting the gene tree frequencies as far plots. 

![](/pics/header_plot.svg)

This package is still in beta testing phase, and more features and expanded functionality will be added in the future. Now, the package can run a standard ASTRAL-III analysis, read in the analysis results, and create pie-chart tree plots. The goal is for this R package to have the same functionality as DiscoVista (https://github.com/esayyari/DiscoVista) plus additional tools, plotting, and quality of life improvements, but within the R environment. 

If you find any issues with something not working, or you would like features to be added, go to issues in the top menu bar and submit them. 

For now, you can cite the R package by linking to this GitHub if you use it. 

Features coming soon:
  1) Branch barplots (alternative to pie charts)
  2) Occupancy plots
  3) Discordance analyses
  4) Relative gene tree frequency analysis
  5) IQTree concordance factor interface and import


# Installation

The major dependency for AstralPlane is of course the program ASTRAL-III. This program is Java-based and can be run on any machine that can run Java. 

ASTRAL-III is available on GitHub here: https://github.com/smirarab/ASTRAL

Instructions for installation and testing ASTRAL-III are included therein, and once it is up and running AstralPlane should be functional! 

The package has two R package dependencies, which are treated as imports (i.e. you need them installed, but library(ape) and library(stringr) not needed: 
  ape (>= 5.0)
  stringr (>= 1.4)
  
And to install AstralPlane, you can use the R package devtools. Here are step-by-step instructions for installation:

1) Install devtools by typing "install.packages(devtools)" in your R console. 

2) Install AstralPlane by typing in your R console: "devtools::install_github("chutter/AstralPlane")"

3) Devtools will ask you to install the package dependecies (ape and stringr), select "Yes". If devtools asks you to update packages, you may choose to do so. I would recommend not to install packages from source if devtools asks you. Ape is problemic from source and I could not get it to install on my machine. If devtools is still giving you trouble, you can install the dependencies with "install.packages(c("ape", "stringr"))". Then rerun Step 2 and skip package updates. 

4) Devtools should finish and say the package loaded properly. Load the package with library(AstralPlane). 

And you should be done! 


# Mini-Vignette: Usage and single dataset example 

I have included an R script in the main repository with some examples. It is also described here in detail. 

1) first install and load the R package. Its a good idea to install new (or check) every time as this package is being updated frequently. 

```r
devtools::install_github("chutter/AstralPlane")
library(AstralPlane)

```

2) You will want a character variable that includes your full path to the astral jar file. NOTE: if you move this file, you will need to move the lib/ directory along with it, as astral depends on it. 


```r
astral.path = "/usr/local/bin/Astral-5-14/astral.5.14.2.jar"
```

3)Setup your working directory and create if necessary

```r
work.dir = "/Test_Astral"
dir.create(work.dir)
setwd(work.dir)
```

4) Next, this step will walk through a single dataset example. 

First, you will want to add a character variable with the path to your gene tree directory. The gene trees here were generated using IQTree and the script should work with other tree types as well, as long as they have branch lengths and support values. Also indicate your outgroups for rooting the tree later. Finally, the output name can be put in a variable or directly entered, your choice. 

```r
tree.dir = "/Trees/Gene_Trees/trimmed_exons"
outgroups = c("Species_A", "Species_B")
save.name = "test-dataset"

```

5) the setupAstral function is used to take your folder of gene trees, apply some filters to the gene trees, and then save them in a single file that can be read by ASTRAL-III. This should take about a minute to run. 


```r
setupAstral(genetree.folder = tree.dir,
            output.name = save.name,
            min.n.samples = 4,
            min.sample.prop = 0,
            make.polytomy = TRUE,
            polytomy.limit = 10)
```

Parameter explanations: 

```
genetree.folder: a folder of genetrees to prepare for astral analyses
output.name: the save name for your concatenated gene tree file
overwrite: whether to overwrite an existing dataset
taxa.remove: species that you would like removed from each gene tree
min.n.samples: the minimum number of samples to keep a gene tree
min.sample.prop: the minimum proportion of samples to keep a gene tree
make.polytomy: whether to collapse poorly supported nodes into polytomies
polytomy.limit: if make.polytomy = TRUE, the threshold value for node collapsing
```

6) When the setup function finishes running, you can now run ASTRAL-III using the runAstral function. This uses the astral jar directly, and should a minute or two depending on your number of gene trees using multi-threading and around 10 without the multi-threading option. 

```r
runAstral(input.genetrees = save.name,
          output.name = save.name,
          astral.path = astral.path,
          astral.t = 2,
          quiet = FALSE,
          load.tree = FALSE,
          multi.thread = TRUE,
          memory = "8g")
```

Parameter explanations: 

```
input.genetrees: a file of genetrees from setupAstral
output.name: the save name for the ASTRAL-III file
astral.path: the absolute file path to the ASTRAL-III jar file
astral.t: the ASTRAL-III "t" parameter for different annotations, t = 2 is all annotation
quiet: hides the screen output from astral if desired
load.tree: TRUE to laod the tree into R
overwrite: TRUE to overwrite an existing dataset
multi.thread: TRUE to use Astral-MP multithreading 
memory: memory value to be passed to java. Should be in "Xg" format, X = an integer
```

7) Next, you can read in the astral data using the astralPlane S4 Object class, which organizes all the analysis data into different slots in the object that can be accessed using the @ symbol. 


```r
astral.data = createAstralPlane(astral.tree = save.name,
                                outgroups = outgroups,
                                tip.length = 1)
```

Parameter explanations: 

```
astral.tree: phylogenetic tree from ape read.tree or a file path to a tree file
outgroups: a vector of outgroups to root the tree
tip.length: arbitrary value for the terminal tip lengths for plotting purposes
```

8) Finally, you can plot your results using the astralProjection function. You give the function the astralPlane object from the previous step, and select your settings for plotting, and what you would like to plot. An example plot is provided in the main Github repository. 

```r
astralProjection(astral.plane = astral.data,
                 local.posterior = TRUE,
                 pie.plot = TRUE,
                 pie.data = "genetree",
                 save.file = "example_plot.pdf",
                 pie.colors = c("purple", "blue", "green"),
                 node.color.text = c("white"),
                 node.color.bg = c("black"),
                 node.label.size = 0.5,
                 tip.label.size = 0.75,
                 pie.chart.size = 1)
```


Parameter explanations: 

```
astral.plane: AstralPlane S4 object of data generated from AstralPlane function
local.posterior: plot the local posterior support = TRUE
pie.plot: TRUE to plot pie charts on branches. FALSE ignores whatever is selected for "pie.data"
pie.data: 'qscore' the quartet support or 'genetree' proportion of gene trees that support a branch
save.file: if you wish to save to file, put file name. Saves as PDF
pie.colors: select three colors to plot your pie.plot
node.color.text: if local.posterior = TRUE, select the color of posterior support text
node.color.bg: if local.posterior = TRUE, select the color of posterior support background
node.label.size: size of the node labels, passed to cex in plotting function
tip.label.size: size of the tip labels, passed to cex in plotting function
pie.chart.size: size of pie chart, passed to edgelabel plotting function

```

![](/pics/header_plot.svg)


# Mini-Vignette: Many dataset example 

This section is a summary of using the AstralPlane package for batch ASTRAL-III data analysis, which incorporates the previous functions in a wrapper to analyze across multiple datasets. By multiple datasets, this means any sets of data you would like to analyze separately, such as exons, introns or UCEs. 

1). To begin, the working directory and paths can be setup like in the previous section

```r
library(AstralPlane)
astral.path = "/usr/local/bin/Astral-5-14/astral.5.14.2.jar"

work.dir = "/Test_Astral"
genetree.folder = "Genetree_Directory"

dir.create(work.dir)
setwd(work.dir)
```

2). Next, all of you datasets should have gene trees estimated from the alignments, and saved to their own directory (i.e. folder called "exons" or "uces"). There should only be the gene trees in this folder, as other files could cause the script to crash. These gene tree data types should be placed in an over-arching "genetree.folder" directory in order to be analyzed together. 


```r

batchAstral(genetree.datasets = genetree.folder,
            astral.t = 2,
            output.dir = "test-dataset",
            min.n.samples = 4,
            min.sample.prop = 0.1,
            taxa.remove = NULL,
            overwrite = TRUE,
            quiet = F,
            astral.path = astral.path,
            make.polytomy = TRUE,
            polytomy.limit = 10,
            multi.thread = TRUE,
            memory = "8g")

```

Parameter explanations: 

```
genetree.datasets: a folder of genetrees to prepare for astral analyses
astral.t: the ASTRAL-III "t" parameter for different annotations, t = 2 is all annotation
output.dir: the save name for your concatenated gene tree file
min.n.samples: the minimum number of samples to keep a gene tree
min.sample.prop: the minimum proportion of samples to keep a gene tree
taxa.remove: species that you would like removed from each gene tree
overwrite: whether to overwrite an existing dataset
quiet: hides the screen output from astral if desired
astral.path: the absolute file path to the ASTRAL-III jar file
make.polytomy: whether to collapse poorly supported nodes into polytomies
polytomy.limit: if make.polytomy = TRUE, the threshold value for node collapsing
multi.thread: TRUE to use Astral-MP multithreading 
memory: memory value to be passed to java. Should be in "Xg" format, X = an integer
```

After this single function is run, the "test-dataset" output directory will be created with all of your astral results in them! You can create figures and plot the results as above. An example below shows a loop saving the plots for all the datasets: 

```r
#Obtains dataset names
datasets = list.dirs(genetree.folder, full.names = F, recursive = F)

for (i in 1:length(datasets)){

  #Read in the astral data and tree and organize it into different slots
  dataset.name = paste0("test-dataset/", datasets[i], "_astral.tre")
  astral.data = createAstralPlane(astral.tree = dataset.name,
                                  outgroups = outgroups,
                                  tip.length = 1)

  #Plots the astral data
  astralProjection(astral.plane = astral.data,
                   local.posterior = TRUE,
                   pie.plot = TRUE,
                   pie.data = "qscore",
                   save.file = paste0("test-dataset/", datasets[i], ".pdf"),
                   pie.colors = c("purple", "blue", "green"),
                   node.color.text = c("white"),
                   node.color.bg = c("black"),
                   tip.label.size = 0.75,
                   pie.chart.size = 1)

}#end i loop
```


More coming soon!

