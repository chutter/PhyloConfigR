#' @title formatAlignmentFolder
#'
#' @description Function for easily converting alignments between different R alignment classes
#'
#' @param alignment.path a folder of individual alignments in nexus or phylip format
#'
#' @param out.dir the name of the new directory with reformatted alignments
#'
#' @param overwrite if TRUE overwrites file if it exists; FALSE the dataset is skipped
#'
#' @param from.format the format to convert the alignment from (phylip, nexus, fasta)
#'
#' @param to.format the format to convert the alignment to (phylip, nexus, fasta)
#'
#' @return saves to file a folder of alignments for the desired format
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

#function
formatAlignmentFolder = function(alignment.path = NULL,
                                 out.dir = NULL,
                                 overwrite = FALSE,
                                 from.format = "nexus",
                                 to.format = "phylip") {

  #Gets list of alignments from path
  if (overwrite == TRUE){
    if (dir.exists(out.dir) == TRUE){ system(paste0("rm -r ", out.dir)) }
    dir.create(out.dir)
  } else {
    if (dir.exists(out.dir) == TRUE){ return("Directory exists and overwrite == FALSE") }
    dir.create(out.dir)
  }#end overwrite if

  if (alignment.path == out.dir){ stop("Error: output directory cannot be the same as input directory.")}

  align.names = list.files(alignment.path)

  #Loops through each alignment to gather statistics
  for (x in 1:length(align.names)){
    #Reads in alignment

    if (from.format == "phylip"){
      save.name = gsub(".phy$", "", align.names[x])
      align = ape::read.dna(paste0(alignment.path, "/", align.names[x]),
                            format = "sequential")
    }

    if (from.format == "nexus"){
      save.name = gsub(".nexus$|.nex$", "", align.names[x])
      l.align = ape::read.nexus.data(paste0(alignment.path, "/", align.names[x]))
      align = ape::as.DNAbin(matrix(unlist(l.align), ncol = length(l.align[[1]]), byrow = TRUE))
      rownames(align) = names(l.align)
    }

    if (from.format == "fasta"){
      save.name = gsub(".fasta$|.fa$", "", align.names[x])
      align = ape::read.FASTA(paste0(alignment.path, "/", align.names[x]), type = "DNA")
    }

    if (to.format == "phylip"){
      writePhylip(alignment = align, file = paste0(out.dir, "/",  save.name, ".phy"))
    }

    if (to.format == "nexus"){
      ape::write.nexus.data(align,  file = paste0(out.dir, "/",  save.name, ".nexus"), interleaved = FALSE)
    }

    #Saves as a fasta file
    if (to.format == "fasta"){
      temp.align = as.character(as.list(align))
      temp.align2 = lapply(temp.align, FUN = function(x) paste(x, collapse = ""))
      align.out = Biostrings::DNAStringSet(unlist(temp.align2))
      write.loci = as.list(as.character(align.out))
      seqinr::write.fasta(sequences = write.loci, names = names(write.loci),
                          paste0(out.dir, "/",  save.name, ".fa"),
                          nbchar = 1000000, as.string = T)
    }#end fasta if

  }#end x loop

}#end function

