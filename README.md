# PhyloConfigR

An R package for preparing, summarizing, and filtering phylogenomic datasets, and for configuring inputs and outputs for phylogenetic software including ASTRAL-III, IQ-TREE 2, BPP, STRUCTURE, and HyPhy.

---

## Overview

PhyloConfigR provides a unified R interface for the core tasks of a modern phylogenomic pipeline:

- Summarizing and filtering alignment datasets
- Estimating and filtering gene trees
- Running coalescent species tree analyses (ASTRAL-III)
- Computing and visualizing IQ-TREE 2 concordance factors
- Detecting anomaly zones
- Preparing inputs for BPP species delimitation and STRUCTURE population analysis
- Running and parsing HyPhy molecular evolution analyses (SLAC, BUSTED, aBSREL, FitMG94)

---

## Citation

Publication in preparation. If you use this package, please cite the GitHub repository: https://github.com/chutter/PhyloConfigR

---

## Installation

### Option 1: Install from GitHub (R only)

Install the R package directly via devtools:

```r
install.packages("devtools", dependencies = TRUE)
devtools::install_github("chutter/PhyloConfigR")
library(PhyloConfigR)
```

R package dependencies (installed automatically): `ape`, `data.table`, `stringr`, `stringi`, `ggplot2`, `jsonlite`, `Biostrings`, `phytools`, `seqinr`

### Option 2: Conda environment (full pipeline)

A conda environment file is provided in `setup-files/environment.yml` that installs all R package dependencies plus the required external programs (IQ-TREE 2, ASTRAL, HyPhy).

```bash
conda env create -f setup-files/environment.yml
conda activate PhyloConfigR
```

Then install the R package inside that environment:

```r
devtools::install_github("chutter/PhyloConfigR")
```

### External software dependencies

The following programs must be installed and accessible on your `PATH` to use the relevant functions:

| Program | Purpose | URL |
|---------|---------|-----|
| IQ-TREE 2 (>=2.0) | Gene tree estimation, concordance factors, concatenation tree | http://www.iqtree.org |
| ASTRAL-III | Coalescent species tree inference | https://github.com/smirarab/ASTRAL |
| Java (>=11) | Required to run ASTRAL | https://adoptium.net |
| HyPhy | Selection analyses (SLAC, BUSTED, aBSREL) | https://hyphy.org |
| BPP | Species delimitation | https://github.com/bpp/bpp |

---

## Function Reference

### Alignment utilities

| Function | Description |
|----------|-------------|
| `summarizeAlignments()` | Calculate per-alignment summary statistics (length, samples, PIS, missing data) |
| `filterAlignments()` | Filter alignments by length, PIS count/proportion, or sample proportion |
| `filterStats()` | Summarize alignment filter results across datasets |
| `filterSummary()` | Generate a combined summary table of alignment filter results |
| `formatAlignmentFolder()` | Rename or reformat alignment files in a folder |
| `alignmentConversion()` | Convert alignments between phylip and nexus formats |
| `concatenateAlignments()` | Concatenate multiple alignments into a single supermatrix |
| `writePhylip()` | Write a DNAbin alignment to phylip format |
| `informativeSites()` | Count or identify parsimony informative sites |
| `heterozygousSites()` | Count heterozygous sites per sample in an alignment |
| `makePolytomy()` | Collapse nodes below a support threshold into polytomies |

### Gene tree estimation and filtering

| Function | Description |
|----------|-------------|
| `estimateGeneTrees()` | Estimate gene trees from a folder of alignments using IQ-TREE 2 |
| `analysis.concatenationTree()` | Estimate a concatenation tree with IQ-TREE 2 |
| `filterGeneTrees()` | Filter gene trees based on alignment statistics |
| `bestFilterTrees()` | Select best gene trees from a filtered set based on posterior probability |

### ASTRAL-III coalescent species trees

| Function | Description |
|----------|-------------|
| `setupAstral()` | Prepare gene tree input files for ASTRAL |
| `runAstral()` | Run ASTRAL-III on a gene tree file |
| `batchAstral()` | Run ASTRAL across multiple filtered datasets |
| `readAstral()` | Read and parse ASTRAL output trees into R |
| `createAstralPlane()` | Create an `AstralPlane` S4 object from an ASTRAL output tree |
| `createAstralPlaneCF()` | Create an `AstralPlane` object incorporating concordance factors |
| `astralProjection()` | Plot an ASTRAL tree with pie charts and posterior support |
| `edgeLengthTable()` | Extract edge length data from an AstralPlane object |

### IQ-TREE 2 concordance factors

| Function | Description |
|----------|-------------|
| `concordanceFactors()` | Compute gene (gCF) and site (sCF) concordance factors with IQ-TREE 2 |
| `readConcordance()` | Read concordance factor output files into R |
| `filterConcordance()` | Filter concordance factor results by node or dataset |

### Anomaly zone detection and visualization

| Function | Description |
|----------|-------------|
| `anomalyZone()` | Detect anomaly zone branches on a tree (Degnan & Rosenberg 2006) |
| `filterAnomalies()` | Summarize anomaly zone counts across filter replicates |
| `plot.anomalyZone()` | Plot anomaly zone status on a phylogeny |
| `plot.filterNode()` | Plot concordance factors vs. filter threshold for a focal node |
| `plot.filterZone()` | Plot concordance factors vs. filter threshold, colored by anomaly zone |

### Species delimitation and population genetics

| Function | Description |
|----------|-------------|
| `generateBPP()` | Generate alignment and IMAP input files for BPP |
| `convert.alignmentsToStructure()` | Convert alignments to STRUCTURE input format |

### HyPhy selection analyses

| Function | Description |
|----------|-------------|
| `hyphy.SLAC()` | Run SLAC (single-likelihood ancestor counting) selection analysis |
| `hyphy.BUSTED()` | Run BUSTED branch-site unrestricted selection test |
| `hyphy.aBSREL()` | Run aBSREL adaptive branch-site random effects test |
| `hyphy.Omega()` | Run FitMG94 to estimate dN/dS (omega) per branch |
| `hyphy.Parser()` | Parse HyPhy JSON output files into R data frames |
| `hyphy.plotSLAC()` | Plot SLAC dN/dS results on a phylogeny |

---

## Typical workflow

```r
library(PhyloConfigR)

# 1. Summarize alignment statistics
align.stats <- summarizeAlignments(
  alignment.path = "alignments/",
  dataset.name   = "exons",
  file.export    = "results/alignment_summary",
  alignment.format = "phylip"
)

# 2. Filter alignments
filterAlignments(
  alignment.path   = "alignments/",
  output.dir       = "filtered/",
  min.length       = 300,
  min.pis          = 3,
  min.taxa         = 0.75,
  alignment.format = "phylip"
)

# 3. Estimate gene trees
estimateGeneTrees(
  alignment.directory = "filtered/",
  output.directory    = "gene-trees/",
  threads             = 8
)

# 4. Run ASTRAL
setupAstral(gene.tree.path = "gene-trees/", output.file = "astral_input.tre")
runAstral(
  astral.path = "/path/to/astral.jar",
  gene.trees  = "astral_input.tre",
  output.file = "astral_output.tre"
)

# 5. Load and plot ASTRAL results
astral.tree <- readAstral("astral_output.tre")
astral.data <- createAstralPlane(astral.tree, outgroups = c("outgroup_sp"), tip.length = 1)
astralProjection(
  astral.plane     = astral.data,
  local.posterior  = TRUE,
  pie.data         = "qscore",
  pie.colors       = c("purple", "blue", "green"),
  save.file        = "astral_plot.pdf"
)

# 6. Compute concordance factors
concordanceFactors(
  species.tree   = "astral_output.tre",
  gene.tree.path = "gene-trees/",
  output.dir     = "concordance/"
)
concord.data <- readConcordance("concordance/")

# 7. Detect anomaly zones
anomaly.data <- anomalyZone(astral.tree, outgroups = c("outgroup_sp"))
plot.anomalyZone(anomaly.data, save.file = "anomaly_zones.pdf")
```

---

## The AstralPlane S4 class

ASTRAL results are stored in an `AstralPlane` S4 object with the following slots:

| Slot | Contents |
|------|----------|
| `@phylo` | The rooted phylo tree object |
| `@samples` | Character vector of tip labels |
| `@outgroups` | Character vector of outgroup tip labels |
| `@nodeData` | Data frame of per-node statistics (pp1, pp2, pp3) |
| `@edgeData` | Data frame of per-edge quartet scores (q1/q2/q3, f1/f2/f3) |

---

## License

GPL-3. See [LICENSE](LICENSE) for details.
