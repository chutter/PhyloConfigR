makeFakeIQTree = function(path){
  script = c(
    "#!/bin/sh",
    "if [ \"$1\" = \"--version\" ]; then",
    "  echo 'IQ-TREE multicore version 2.3.6 for MacOS ARM 64-bit built Jul 30 2024'",
    "  exit 0",
    "fi",
    "while [ \"$#\" -gt 0 ]; do",
    "  case \"$1\" in",
    "    -pre) prefix=\"$2\"; shift 2 ;;",
    "    *) shift ;;",
    "  esac",
    "done",
    "printf '(a:0.1,b:0.1,c:0.1,d:0.1);\\n' > \"${prefix}.treefile\"",
    "printf 'test output\\n' > \"${prefix}.iqtree\""
  )
  writeLines(script, path)
  Sys.chmod(path, mode = "0755")
}

makeAlignment = function(path, n.taxa = 4){
  sequences = c("a ACGT", "b ACGA", "c ATGT", "d TCGT")
  writeLines(c(paste(n.taxa, 4), sequences[seq_len(n.taxa)]), path)
}

test_that("gene trees are estimated and resumed with paths containing spaces", {
  skip_on_os("windows")
  work.dir = tempfile("gene-tree-test-")
  alignment.dir = file.path(work.dir, "alignments with spaces")
  output.dir = file.path(work.dir, "gene trees")
  bin.dir = file.path(work.dir, "bin with spaces")
  dir.create(alignment.dir, recursive = TRUE)
  dir.create(bin.dir)
  makeAlignment(file.path(alignment.dir, "locus.1.phy"))
  makeFakeIQTree(file.path(bin.dir, "iqtree2"))

  result = estimateGeneTrees(alignment.dir, output.dir,
                             iqtree.path = bin.dir,
                             cleanup.files = FALSE)
  expect_equal(result$status, "success")
  expect_true(file.exists(file.path(output.dir, "locus.1.phy.treefile")))
  expect_true(file.exists(file.path(output.dir, "locus.1.phy.iqtree")))

  resumed = estimateGeneTrees(alignment.dir, output.dir,
                              iqtree.path = bin.dir)
  expect_equal(resumed$status, "skipped_complete")

  clean.dir = file.path(work.dir, "clean trees")
  cleaned = estimateGeneTrees(alignment.dir, clean.dir,
                              iqtree.path = bin.dir,
                              cleanup.files = TRUE)
  expect_equal(cleaned$status, "success")
  expect_true(file.exists(file.path(clean.dir, "locus.1.phy.treefile")))
  expect_false(file.exists(file.path(clean.dir, "locus.1.phy.iqtree")))
})

test_that("proportional chunks do not overlap", {
  skip_on_os("windows")
  work.dir = tempfile("gene-tree-test-")
  alignment.dir = file.path(work.dir, "alignments")
  bin.dir = file.path(work.dir, "bin")
  dir.create(alignment.dir, recursive = TRUE)
  dir.create(bin.dir)
  makeFakeIQTree(file.path(bin.dir, "iqtree2"))
  for (i in 1:4){ makeAlignment(file.path(alignment.dir, paste0("locus", i, ".phy"))) }

  first = estimateGeneTrees(alignment.dir, file.path(work.dir, "first"),
                            subset.start = 0, subset.end = 0.5,
                            iqtree.path = bin.dir)
  second = estimateGeneTrees(alignment.dir, file.path(work.dir, "second"),
                             subset.start = 0.5, subset.end = 1,
                             iqtree.path = bin.dir)

  expect_equal(first$locus, c("locus1.phy", "locus2.phy"))
  expect_equal(second$locus, c("locus3.phy", "locus4.phy"))
  expect_length(intersect(first$locus, second$locus), 0)
})

test_that("minimum taxa and destructive output paths are checked", {
  skip_on_os("windows")
  work.dir = tempfile("gene-tree-test-")
  alignment.dir = file.path(work.dir, "alignments")
  output.dir = file.path(work.dir, "trees")
  bin.dir = file.path(work.dir, "bin")
  dir.create(alignment.dir, recursive = TRUE)
  dir.create(bin.dir)
  makeAlignment(file.path(alignment.dir, "small.phy"), n.taxa = 3)
  makeFakeIQTree(file.path(bin.dir, "iqtree2"))

  result = estimateGeneTrees(alignment.dir, output.dir,
                             iqtree.path = bin.dir)
  expect_equal(result$status, "skipped_min_taxa")
  expect_error(
    estimateGeneTrees(alignment.dir, alignment.dir,
                      overwrite = TRUE, resume = FALSE,
                      iqtree.path = bin.dir),
    "cannot be the alignment directory"
  )
})
