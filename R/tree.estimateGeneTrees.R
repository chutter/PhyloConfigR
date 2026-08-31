#' @title estimateGeneTrees
#'
#' @description Batch estimates gene trees for a folder of phylip alignments
#'   using IQ-TREE 2. Runs ModelFinder and UFBoot for each alignment, supports
#'   deterministic job subsets and resuming, and records a run manifest.
#'
#' @param alignment.directory path to a folder of alignment files in phylip format
#'
#' @param output.directory name of the output directory for tree files (default: "gene-trees")
#'
#' @param min.taxa minimum number of taxa required to attempt tree estimation (default: 4)
#'
#' @param subset.start proportion (0-1) of the complete alignment list at which
#'   this job starts (default: 0)
#'
#' @param subset.end proportion (0-1) of the complete alignment list at which
#'   this job ends (default: 1)
#'
#' @param threads number of CPU threads passed to IQ-TREE (-nt flag), or "AUTO"
#'
#' @param memory maximum memory in GB passed to IQ-TREE (-mem flag)
#'
#' @param overwrite if TRUE deletes and recreates the output directory
#'
#' @param resume if TRUE skips alignments whose valid treefile already exists
#'
#' @param quiet if TRUE writes IQ-TREE screen output to a per-locus runner log
#'
#' @param cleanup.files if TRUE removes IQ-TREE auxiliary files after a valid
#'   tree is produced; defaults to FALSE so model and bootstrap information is kept
#'
#' @param iqtree.path path to the iqtree2 executable or the directory containing
#'   it; use NULL if iqtree2 is on the system PATH
#'
#' @param model model or ModelFinder command passed to IQ-TREE (default: "MFP")
#'
#' @param uf.bootstrap number of UFBoot replicates; use 0 to disable (default: 1000)
#'
#' @param seed optional IQ-TREE random seed for reproducible runs
#'
#' @return invisibly returns a data.frame summarizing each locus. The same table
#'   is written to a chunk-specific CSV file in output.directory.
#'
#' @examples
#'
#' estimateGeneTrees(alignment.directory = "alignments",
#'                   output.directory = "gene-trees",
#'                   min.taxa = 4,
#'                   threads = 4,
#'                   overwrite = FALSE,
#'                   resume = TRUE,
#'                   quiet = TRUE)
#'
#' @export

estimateGeneTrees = function(alignment.directory = NULL,
                             output.directory = "gene-trees",
                             min.taxa = 4,
                             subset.start = 0,
                             subset.end = 1,
                             threads = 1,
                             memory = 1,
                             overwrite = FALSE,
                             resume = TRUE,
                             quiet = TRUE,
                             cleanup.files = FALSE,
                             iqtree.path = NULL,
                             model = "MFP",
                             uf.bootstrap = 1000,
                             seed = NULL) {

  #Initial checks
  if (is.null(alignment.directory)){ stop("A folder of alignments is needed.") }
  if (is.null(output.directory)){ stop("An output directory is needed.") }
  if (dir.exists(alignment.directory) == FALSE){ stop("Alignment directory could not be found.") }
  if (length(min.taxa) != 1 || is.numeric(min.taxa) == FALSE || min.taxa < 4){ stop("min.taxa must be at least 4.") }
  if (length(subset.start) != 1 || length(subset.end) != 1 ||
      is.numeric(subset.start) == FALSE || is.numeric(subset.end) == FALSE ||
      subset.start < 0 || subset.end > 1 || subset.start >= subset.end){
    stop("subset.start and subset.end must define an increasing interval between 0 and 1.")
  }
  if (length(memory) != 1 || is.numeric(memory) == FALSE || memory <= 0){ stop("memory must be greater than 0 GB.") }
  if (length(model) != 1 || is.character(model) == FALSE || nzchar(model) == FALSE){ stop("A valid IQ-TREE model is needed.") }
  if (grepl("MERGE", model, ignore.case = TRUE)){ stop("MERGE models require a partition file; use model = 'MFP' for individual alignments.") }
  if (length(uf.bootstrap) != 1 || is.numeric(uf.bootstrap) == FALSE || uf.bootstrap < 0){ stop("uf.bootstrap must be 0 or greater.") }
  if (is.null(seed) == FALSE && (length(seed) != 1 || is.numeric(seed) == FALSE || seed < 1)){ stop("seed must be NULL or a positive number.") }

  if (is.numeric(threads) && length(threads) == 1 && threads >= 1){
    threads = as.character(as.integer(threads))
  } else if (length(threads) == 1 && toupper(as.character(threads)) == "AUTO"){
    threads = "AUTO"
  } else { stop("threads must be a positive number or 'AUTO'.") }

  #Finds IQ-TREE
  if (is.null(iqtree.path)){
    iqtree.exe = unname(Sys.which("iqtree2"))
  } else if (dir.exists(iqtree.path)){
    iqtree.exe = file.path(iqtree.path, "iqtree2")
  } else { iqtree.exe = iqtree.path }
  if (length(iqtree.exe) == 0 || nzchar(iqtree.exe) == FALSE || file.exists(iqtree.exe) == FALSE){
    stop("The iqtree2 executable could not be found.")
  }

  #Protects the alignments from overwrite
  alignment.path = normalizePath(alignment.directory, winslash = "/", mustWork = TRUE)
  output.path = normalizePath(output.directory, winslash = "/", mustWork = FALSE)
  output.prefix = paste0(sub("/+$", "", output.path), "/")
  if (identical(alignment.path, output.path) || startsWith(alignment.path, output.prefix)){
    stop("output.directory cannot be the alignment directory or one of its parent directories.")
  }
  if (resume == TRUE && overwrite == TRUE){
    stop("resume = TRUE and overwrite = TRUE cannot be used together.")
  }

  if (dir.exists(output.path)){
    if (overwrite == TRUE){
      unlink(output.path, recursive = TRUE, force = TRUE)
    } else if (resume == FALSE && length(list.files(output.path, all.files = TRUE, no.. = TRUE)) > 0){
      stop("The output directory is not empty; use resume = TRUE or overwrite = TRUE.")
    }
  }
  if (dir.exists(output.path) == FALSE){ dir.create(output.path, recursive = TRUE) }
  if (dir.exists(output.path) == FALSE){ stop("The output directory could not be created.") }

  #Gathers a deterministic list of phylip alignments
  locus.files = list.files(alignment.path,
                           pattern = "\\.(phy|phylip)$",
                           full.names = TRUE,
                           ignore.case = TRUE)
  file.info.data = file.info(locus.files)
  locus.files = locus.files[is.na(file.info.data$isdir) == FALSE & file.info.data$isdir == FALSE]
  locus.files = locus.files[order(basename(locus.files))]
  if (length(locus.files) == 0){ stop("No phylip alignments were found.") }

  #Assigns the subset before checking completed trees so jobs do not shift
  sub.start = floor(subset.start * length(locus.files)) + 1
  sub.end = if (subset.end == 1) length(locus.files) else floor(subset.end * length(locus.files))
  if (sub.start <= sub.end){
    locus.files = locus.files[seq.int(sub.start, sub.end)]
  } else { locus.files = character() }

  validTree = function(tree.file, expected.taxa){
    if (file.exists(tree.file) == FALSE || file.info(tree.file)$size == 0){ return(FALSE) }
    test.tree = tryCatch(ape::read.tree(tree.file), error = function(e) NULL)
    is.null(test.tree) == FALSE && length(test.tree$tip.label) >= expected.taxa
  }

  makeResult = function(locus, taxa = NA_integer_, sites = NA_integer_, status,
                        exit.status = NA_integer_, attempts = 0L, elapsed = 0,
                        tree.file, message = "", command = ""){
    data.frame(locus = locus, taxa = taxa, sites = sites, status = status,
               exit_status = exit.status, attempts = attempts,
               elapsed_seconds = round(elapsed, 3), tree_file = tree.file,
               message = message, command = command,
               stringsAsFactors = FALSE)
  }

  results = vector("list", length(locus.files))

  #Runs each alignment
  for (i in seq_along(locus.files)) {

    locus.file = locus.files[i]
    locus.name = basename(locus.file)
    output.base = file.path(output.path, locus.name)
    tree.file = paste0(output.base, ".treefile")

    header = tryCatch(scan(locus.file, what = integer(), nmax = 2, quiet = TRUE),
                      error = function(e) integer())
    if (length(header) < 2){
      results[[i]] = makeResult(locus.name, status = "invalid_alignment",
                                tree.file = tree.file,
                                message = "Could not read the phylip header.")
      next
    }
    n.taxa = header[1]
    n.sites = header[2]

    if (n.taxa < min.taxa){
      results[[i]] = makeResult(locus.name, n.taxa, n.sites, "skipped_min_taxa",
                                tree.file = tree.file)
      next
    }
    if (resume == TRUE && validTree(tree.file, n.taxa)){
      results[[i]] = makeResult(locus.name, n.taxa, n.sites, "skipped_complete",
                                tree.file = tree.file)
      next
    }

    iqtree.args = c("-s", locus.file,
                    "-pre", output.base,
                    "-nt", threads,
                    "-mem", paste0(memory, "G"),
                    "-m", model)
    if (uf.bootstrap > 0){ iqtree.args = c(iqtree.args, "-bb", as.integer(uf.bootstrap)) }
    if (is.null(seed) == FALSE){ iqtree.args = c(iqtree.args, "-seed", as.integer(seed)) }
    if (file.exists(tree.file) && validTree(tree.file, n.taxa) == FALSE){ iqtree.args = c(iqtree.args, "-redo") }

    command = paste(shQuote(c(iqtree.exe, iqtree.args)), collapse = " ")
    runner.log = paste0(output.base, ".runner.log")
    error.message = ""
    start.time = proc.time()[3]

    run.status = tryCatch({
      if (quiet == TRUE){
        system2(iqtree.exe, args = vapply(iqtree.args, shQuote, character(1)),
                stdout = runner.log, stderr = runner.log)
      } else {
        system2(iqtree.exe, args = vapply(iqtree.args, shQuote, character(1)))
      }
    }, error = function(e){
      error.message <<- conditionMessage(e)
      127L
    })
    elapsed = proc.time()[3] - start.time

    if (run.status == 0 && validTree(tree.file, n.taxa)){
      results[[i]] = makeResult(locus.name, n.taxa, n.sites, "success",
                                run.status, 1L, elapsed, tree.file,
                                command = command)

      if (cleanup.files == TRUE){
        auxiliary.suffixes = c(".bionj", ".ckp.gz", ".contree", ".iqtree",
                               ".log", ".mldist", ".model.gz", ".runner.log",
                               ".splits.nex", ".ufboot", ".uniqueseq.phy")
        auxiliary.files = file.path(output.path,
                                    paste0(locus.name, auxiliary.suffixes))
        auxiliary.files = auxiliary.files[file.exists(auxiliary.files)]
        if (length(auxiliary.files) > 0){ unlink(auxiliary.files) }
      }
    } else {
      if (nzchar(error.message) == FALSE){
        if (run.status != 0){
          error.message = paste0("IQ-TREE exited with status ", run.status, ".")
        } else { error.message = "IQ-TREE did not produce a valid treefile." }
      }
      results[[i]] = makeResult(locus.name, n.taxa, n.sites, "failed",
                                run.status, 1L, elapsed, tree.file,
                                error.message, command)
    }
  }

  if (length(results) == 0){
    result.table = data.frame(locus = character(), taxa = integer(), sites = integer(),
                              status = character(), exit_status = integer(), attempts = integer(),
                              elapsed_seconds = numeric(), tree_file = character(),
                              message = character(), command = character(),
                              stringsAsFactors = FALSE)
    manifest.name = "gene-tree-manifest-empty.csv"
  } else {
    result.table = do.call(rbind, results)
    manifest.name = paste0("gene-tree-manifest_", sub.start, "-", sub.end, ".csv")
  }
  manifest.file = file.path(output.path, manifest.name)
  utils::write.csv(result.table, manifest.file, row.names = FALSE)
  attr(result.table, "manifest.file") = manifest.file

  invisible(result.table)
}#end function

### END SCRIPT
