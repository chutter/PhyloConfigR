#' @title estimateGeneTrees
#'
#' @description Function for batch trimming a folder of alignments, with the various trimming functions available to select from
#'
#' @param assembly.directory path to a folder of sequence alignments in phylip format.
#'
#' @param target.file available input alignment formats: fasta or phylip
#'
#' @param alignment.contig.name contigs are added into existing alignment if algorithm is "add"
#'
#' @param output.directory available output formats: phylip
#'
#' @param min.percent.id algorithm to use: "add" add sequences with "add.contigs"; "localpair" for local pair align. All others available
#'
#' @param min.match.length TRUE applies the adjust sequence direction function of MAFFT
#'
#' @param min.match.coverage if a file name is provided, save.name will be used to save aligment to file as a fasta
#'
#' @param threads path to a folder of sequence alignments in phylip format.
#'
#' @param memory give a save name if you wnat to save the summary to file.
#'
#' @param trim.target TRUE to supress mafft screen output
#'
#' @param overwrite path to a folder of sequence alignments in phylip format.
#'
#' @param resume contigs are added into existing alignment if algorithm is "add"
#'
#' @param quiet algorithm to use: "add" add sequences with "add.contigs"; "localpair" for local pair align. All others available
#'
#' @param blast.path algorithm to use: "add" add sequences with "add.contigs"; "localpair" for local pair align. All others available
#'
#' @param bbmap.path algorithm to use: "add" add sequences with "add.contigs"; "localpair" for local pair align. All others available
#'
#' @return an alignment of provided sequences in DNAStringSet format. Also can save alignment as a file with save.name
#'
#' @examples
#'
#' your.tree = ape::read.tree(file = "file-path-to-tree.tre")
#' astral.data = astralPlane(astral.tree = your.tree,
#'                           outgroups = c("species_one", "species_two"),
#'                           tip.length = 1)
#'
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
                             cleanup.files = TRUE,
                             iqtree.path = NULL) {

  # #Debug setup
  # library(PhyloCap)
  # setwd("/Users/chutter/Dropbox/Research/1_Main-Projects/1_Collaborative-Projects/Shrew_Genome")
  # alignment.directory = "alignments"
  # output.directory = "gene-trees"
  #
  # # #Main settings
  # subset.start = 0
  # subset.end = 1
  # min.taxa = 4
  # threads = 4
  # memory = 8
  # overwrite = FALSE
  # resume = TRUE
  # quiet = TRUE
  # cleanup.files = TRUE
  #
  # # #program paths
  # iqtree.path = "/usr/local/bin"

  #Same adds to bbmap path
  if (is.null(iqtree.path) == FALSE){
    b.string = unlist(strsplit(iqtree.path, ""))
    if (b.string[length(b.string)] != "/") {
      iqtree.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { iqtree.path = NULL }

  #Initial checks
  if (is.null(alignment.directory) == T){ stop("A folder of alignments is needed.") }
  if (is.null(output.directory) == T){ stop("An output directory is needed.") }

  #So I don't accidentally delete everything while testing resume
  if (resume == TRUE & overwrite == TRUE){
    overwrite = FALSE
    stop("Error: resume = T and overwrite = T, cannot resume if you are going to delete everything!")
  }

  if (dir.exists(output.directory) == TRUE) {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.directory))
      dir.create(output.directory)
    }
  } else { dir.create(output.directory) }

  #gathers files and resumes
  locus.names = list.files(alignment.directory)

  #Checks for alignments already done and removes from the to-do list
  if (resume == TRUE){
    done = list.files(output.directory)
    locus.names = locus.names[!locus.names %in% gsub(".phy$", "", done)]
  }
  if (length(locus.names) == 0){ quit() }

  #Figures out start and end of subset
  sub.start = floor(subset.start * length(locus.names))
  if (sub.start == 0){ sub.start = 1}
  sub.end = floor(subset.end * length(locus.names))

  #Loops through each locus and writes each species to end of file
  for (i in sub.start:sub.end) {

    #Checks if the tree file exists
    if (file.exists(paste0(output.directory, "/", locus.names[i], ".treefile")) == T){ next }

    #create output place
    system(paste0("cp ", alignment.directory, "/", locus.names[i], " ",
                  output.directory, "/", locus.names[i]))

    #Run IQ tree on tree
    system(paste0(iqtree.path, "iqtree2 -s ", output.directory,"/", locus.names[i],
                  " -bb 1000 -nt ", threads, " -m MFP+MERGE -rcluster 20 -msub nuclear"), ignore.stdout = quiet)

    #Delete extra files
    if (file.exists(paste0(output.directory, "/", locus.names[i], ".treefile")) == F){
      print(paste0(locus.names[i], " tree failed, will retry with AUTO threads."))
      #Run IQ tree on tree
      system(paste0(iqtree.path, "iqtree2 -s ", output.directory,"/", locus.names[i],
                    " -bb 1000 -nt AUTO -m MFP+MERGE -rcluster 20 -msub nuclear"), ignore.stdout = quiet)
    }#end if

    del.files = dir(path = output.directory)
    temp.delete = del.files[grep(pattern = locus.names[i], x = del.files, invert=F)]
    to.delete = temp.delete[grep(pattern = "treefile$", x = temp.delete, invert=T)]
    if (length(to.delete) != 0){ unlink(paste0(output.directory, "/", to.delete)) }

  } #i loopo end

}#end function

### END SCRIPT
