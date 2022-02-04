#' @title hyphy.SLAC
#'
#' @description Function for batch trimming a folder of alignments, with the various trimming functions available to select from
#'
#' @param tree.directory path to a folder of sequence alignments in phylip format.
#'
#' @param alignment.directory available input alignment formats: fasta or phylip
#'
#' @param dataset.name contigs are added into existing alignment if algorithm is "add"
#'
#' @param threads path to a folder of sequence alignments in phylip format.
#'
#' @param memory give a save name if you wnat to save the summary to file.
#'
#' @param overwrite TRUE to supress mafft screen output
#'
#' @param resume TRUE to supress mafft screen output
#'
#' @param hyphy.path TRUE to supress mafft screen output
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

hyphy.SLAC = function(tree.directory = NULL,
                      alignment.directory = NULL,
                      dataset.name = "SLAC",
                      threads = 1,
                      memory = 1,
                      overwrite = FALSE,
                      resume = TRUE,
                      quiet = FALSE,
                      hyphy.path = NULL) {


  #Read in basic genome info
  # library(PhyloCap)
  # setwd("/Users/chutter/Dropbox/Research/1_Main-Projects/0_Working-Projects/Rodent_Mitochondrial")
  # tree.directory= "/Users/chutter/Dropbox/Research/1_Main-Projects/0_Working-Projects/Rodent_Mitochondrial/Align-Trees/nuclear/trees_oxphos-coding"
  # alignment.directory = "/Users/chutter/Dropbox/Research/1_Main-Projects/0_Working-Projects/Rodent_Mitochondrial/Align-Trees/nuclear/alignments_oxphos-coding"
  # #metadata.file = "/Volumes/Rodents/Australian_Rodents/Data_Processing/Mus-selected-sequences_metadata_final.csv"
  # dataset.name = "Busted"
  # threads = 4
  # memory = 4
  # resume = T
  # overwrite = F
  # quiet = TRUE
  # hyphy.path = "/usr/local/bin"

  #Same adds to bbmap path
  if (is.null(hyphy.path) == FALSE){
    b.string = unlist(strsplit(hyphy.path, ""))
    if (b.string[length(b.string)] != "/") {
      hyphy.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { hyphy.path = "" }

  if (is.null(tree.directory) == T){ stop("A directory of trees is needed.") }
  if (is.null(alignment.directory) == T){ stop("A directory of alignments is needed.") }

  #So I don't accidentally delete everything while testing resume
  if (resume == TRUE & overwrite == TRUE){
    overwrite = FALSE
    stop("Error: resume = T and overwrite = T, cannot resume if you are going to delete everything!")
  }

  #Creates main analysis directory if its not already there
  if (dir.exists("hyphy") == FALSE) { dir.create("hyphy") }

  output.directory = paste0("hyphy/", dataset.name)
  if (dir.exists(output.directory) == TRUE){
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.directory))
      dir.create(output.directory)
    }
  } else { dir.create(output.directory) }

  align.files = list.files(alignment.directory)
  tree.files = list.files(tree.directory)
  #meta.data = read.csv(metadata.file)

  #Resumes file download
  if (resume == TRUE){
    done.files = list.files(output.directory)
    tree.files = tree.files[!gsub("\\..*", "", tree.files) %in% done.files]
  }

  #Run analyses through all trees
  for (i in 1:length(tree.files)){

    #Find correct gene tree
    locus.name = gsub("\\..*", "", tree.files[i])
  #  temp.meta = meta.data[meta.data$gene_id %in% locus.name,]
   # gene.name = unique(temp.meta$gene_id)
  #  prot.name = unique(temp.meta$protein_id)
   # prot.name = prot.name[is.na(prot.name) != T]

    #if (length(prot.name) == 0){
    #  print(paste0(locus.name, " protein not found in the alignments. Skipped."))
    #  next }

    align.name = align.files[grep(locus.name, align.files)]

    if (length(align.name) == 0){
      print(paste0(locus.name, " not found in the alignments. Skipped."))
      next }

    #Read in alignment
    temp.align = Biostrings::readDNAStringSet(paste0(alignment.directory, "/", align.name))
    #Prune gene tree to match alignment
    temp.tree = ape::read.tree(paste0(tree.directory, "/", tree.files[i]))

    #First remove sample from alignment not in tree
    temp.align = temp.align[names(temp.align) %in% temp.tree$tip.label]
    #Remove tips from tree not in alignment
    temp.tree = ape::drop.tip(temp.tree, temp.tree$tip.label[!temp.tree$tip.label %in% names(temp.align)])

    if (length(temp.align) <= 3){
      print(paste0(locus.name, " alignment has too few taxa. Skipped."))
      next }

    if (length(temp.tree$tip.label) <= 3){
      print(paste0(locus.name, " alignment has too few taxa. Skipped."))
      next }

    if (is.null(temp.tree) == TRUE){
      print(paste0(locus.name, " alignment has too few taxa. Skipped."))
      next }

    dir.create(paste0(output.directory, "/", locus.name))
    #Save both as temp files to be used as input in hyphy
    write.temp = strsplit(as.character(temp.align), "")
    aligned.set = as.matrix(ape::as.DNAbin(write.temp) )
    writePhylip(alignment = aligned.set,
                file=paste0(output.directory, "/", locus.name, "/", locus.name, ".phy"),
                interleave = F,
                strict = F)

    ape::write.tree(temp.tree, file = paste0(output.directory, "/", locus.name, "/", locus.name, ".tre"))

    #Runs the busted
    system(paste0(hyphy.path, "hyphy slac --alignment ", output.directory, "/", locus.name, "/", locus.name, ".phy",
                                     " --tree ", output.directory, "/", locus.name, "/", locus.name, ".tre",
                                     " --output ", output.directory, "/", locus.name, "/SLAC-results.json"),
           ignore.stderr = quiet, ignore.stdout = quiet)
  }#end i loop


}#end function


