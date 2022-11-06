#' @title alignmentConversion
#'
#' @description Function for easily converting alignments between different R alignment classes
#'
#' @param input.alignment alignment file read into R in almost any format
#'
#' @param end.format the format to convert the alignment to. Options include: DNABin, DNAStringSet, matrix, list
#'
#' @return returns a new alignment from the desired format
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

#Converts alignments to various different formats to be used with a variety of functions
alignmentConversion = function(input.alignment = NULL,
                                end.format = NULL) {

  #input.alignment = alignment
  #end.format = "matrix"

  #Parameter checks
  if(is.null(input.alignment) == TRUE){ stop("Error: an input alignment is needed.") }
  if(is.null(end.format) == TRUE){ stop("Error: an ending format is needed.") }

  align.class = class(input.alignment)

  #Converts from DNAStringSets to DNABin Objects
  if (align.class == "DNAStringSet" && end.format == "DNAbin"){
    new.align = strsplit(as.character(input.alignment), "")
    mat.align = lapply(new.align, tolower)
    align.out = as.matrix(as.DNAbin(mat.align))
    return(align.out)
  }# end if

  #Converts from DNABin to DNAStringSet Objects
  if (align.class == "DNAbin" && end.format == "DNAStringSet"){

    # temp.loci = alignmentConversion(input.alignment = align, end.format = "DNAStringSet")
    #  This might be the correct way to do it if fails again
    # temp.align = as.character(as.list(align))
    # temp.align2 = lapply(temp.align, FUN = function(x) paste(x, collapse = ""))
    # align.out = Biostrings::DNAStringSet(unlist(temp.align2))
    # write.loci = as.list(as.character(temp.loci))

    temp.align = as.list(data.frame(t(input.alignment)))
    temp.align2 = lapply(temp.align, FUN = function(x) paste(x, collapse = ""))
    align.out = DNAStringSet(unlist(temp.align2))
    return(align.out)
  }# end if

  #Converts from DNABin to DNAStringSet Objects
  if (align.class == "DNAStringSet" && end.format == "matrix"){
    temp.align = strsplit(as.character(input.alignment), "")
    align.out = as.matrix(as.DNAbin(temp.align) )
    return(align.out)
  }# end if

  #Converts from DNABin to DNAStringSet Objects
  if (align.class == "DNAbin" && end.format == "matrix"){
    align.out = as.matrix(input.alignment)
    return(align.out)
  }# end if

  ##### FINISH THIS PART
  stop("Alignment class or conversion format not found.")

}#end function
