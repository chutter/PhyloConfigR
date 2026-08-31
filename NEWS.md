# PhyloConfigR 0.3.1

## Fixes

* `analysis.concatenationTree()` no longer passes `-msub` on unpartitioned runs
  (`partition.scheme = "none"`). `-msub` only restricts ModelFinder's amino-acid
  model set, so it has no effect when an explicit model is given, and some
  IQ-TREE builds reject the flag in that context with `ERROR: Unknown sequence
  type` (exit status 2), which made every replicate of a gene jackknife fail at
  IQ-TREE launch. It is now added only for the ModelFinder schemes ("merge",
  "file").

* `analysis.concatenationTree()` gained a `seq.type` argument that forces the
  IQ-TREE data type with `-st`, and `analysis.geneJackknife()` now sets it to
  `"DNA"` by default. Some IQ-TREE builds fail auto-detection with `ERROR:
  Unknown sequence type` on replicate matrices that carry a lot of missing data;
  forcing the type is deterministic and, if a locus were genuinely not DNA,
  reports the offending site instead of the opaque error. Set `seq.type = NULL`
  to auto-detect.

* `analysis.geneJackknife()` now removes each replicate's concatenated matrix
  even when the replicate fails, not only on success. The matrix is built just
  before its tree and deleted in a `finally` block, so only the replicates in
  flight occupy disk and an aborted run no longer leaves a directory full of
  matrices behind.

# PhyloConfigR 0.3.0

Adds a gene jackknife.

## New

* `analysis.geneJackknife()` runs a gene jackknife on a folder of locus
  alignments. Each replicate draws loci at random without replacement to a fixed
  size (a number of base pairs or a number of genes), concatenates them with
  `concatenateAlignments()`, and estimates a tree with
  `analysis.concatenationTree()` with UFBoot off, because the resampling is the
  replication. It builds a majority-rule consensus from the replicate trees, or
  can be pointed at a single index with `replicate.subset` so a cluster array
  runs one replicate per task. Replicate `i` is seeded from `seed + i`, so the
  result does not depend on the order the tasks run and any replicate reproduces
  on its own. An optional `locus.lengths` table skips re-reading the alignments
  on every task.

# PhyloConfigR 0.2.1

Fixes to IQ-TREE handling. All four affected functions previously assumed the
executable is named `iqtree2`, which is true of IQ-TREE 2 and not of IQ-TREE 3.

## New

* `findIQTREE()` locates an IQ-TREE executable and reports its version. Tries
  `iqtree2` then `iqtree`, and reads the version back from the executable rather
  than inferring it from the file name. A version that cannot be read is
  reported as `NA` rather than being treated as an error, so wrapper scripts and
  unusual builds still run.

## Fixes

* `analysis.concatenationTree()` could not be called on its default arguments.
  `partition.scheme` and `msub.type` defaulted to their whole choice vectors,
  and comparing a length-3 vector with `==` inside `if()` is an error in R 4.2
  and later. Both now go through `match.arg()`, and `-msub` no longer expands to
  two commands.

* `analysis.concatenationTree()` always passed `-bb`, and IQ-TREE rejects `-bb`
  below 1000, so bootstrapping could not be turned off. `uf.bootstrap = 0` now
  omits the flag, matching the convention `estimateGeneTrees()` already used.
  This makes the function usable inside resampling procedures that do their own
  replication, such as a gene jackknife.

* `analysis.concatenationTree()` ignored the exit status of IQ-TREE and could
  report success after a failed run. It now checks the status and that a
  treefile was produced, and returns the treefile path invisibly.

* `analysis.concatenationTree()` gained a `model` argument for
  `partition.scheme = "none"`, which previously hard-coded `GTR`, along with
  `seed`, `-mem`, a working `resume`, and recursive directory creation.

* `estimateGeneTrees()`, `concordanceFactors()` and `concordanceRunner()` now
  resolve the executable with `findIQTREE()` instead of assuming `iqtree2`.
  `concordanceFactors()` and `concordanceRunner()` take `iqtree.path = NULL` by
  default.

* `NAMESPACE` now imports data.table. Listing it under `Imports` in DESCRIPTION
  is not enough: data.table's `:=` checks that the calling namespace is
  data.table-aware, so `concatenateAlignments()` failed at run time with "[ was
  called on a data.table in an environment that is not data.table-aware" for
  every installed copy of the package. Sourcing the R files hid the problem,
  because the calls then happened in the global environment.

# PhyloConfigR 0.2.0

* Initial packaged release.
