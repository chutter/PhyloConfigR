#' @title writePhylip
#'
#' @description Function for writing alignments in R to phylip format
#'
#' @param alignment summary data file from filterSummary
#'
#' @param file summary data file from alignmentSummary
#'
#' @param interleave folder of alignments to be filtered
#'
#' @param strict save format for alignments
#'
#' @return saves the alignment as a phylip file
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

##### saves alignment as a phylip file
writePhylip = function(alignment = NULL,
                       file = NULL,
                       interleave = FALSE,
                       strict = FALSE){

  x = as.matrix(alignment)
  #file = "test"

  #Parameter checks
  if(is.null(alignment) == TRUE){ stop("Error: an input alignment is needed.") }
  if(is.null(file) == TRUE){ stop("Error: a file name is needed.") }
  if(is.null(row.names(alignment)) == TRUE){ stop("Error: no row names found on matrix.") }

  #converts alignment
  # new.align = alignmentConversion(input.alignment = alignment, end.format = "matrix")
  #
  # ntax = nrow(new.align)
  # nchar = ncol(new.align)
  # phy.header = paste(ntax, nchar)

  # #Now have to reformat somehow
  # #Prep sample names to all have the same length padded with spaces
  # name.length = max(nchar(sample.names)) + 5
  # for (x in 1:length(sample.names)){
  #   temp.name = concat.data$Sample[x]
  #   space.pad = name.length - nchar(temp.name)
  #   space.add = paste(rep(" ", space.pad), collapse = "")
  #   new.name = paste0(temp.name, space.add, collapse = "")
  #   concat.data$Sample[x] = new.name
  # }#end x

  str2cha <- function(x) {
    unlist(strsplit(x, ""))
  }
  datatype <- ifelse(is.numeric(x[1, 1]), "continuous", "nc")
  ntax <- nrow(x)
  nchar <- ncol(x)
  taxnames <- rownames(x)
  if (strict) {
    taxnames <- substring(taxnames, 1, truncate)
    missing <- 10 - unlist(lapply(strsplit(taxnames, ""),
                                  length))
    for (i in seq(along = taxnames)) taxnames[i] <- paste(taxnames[i],
                                                          paste(rep("*", missing[i]), collapse = ""), sep = "")
    if (any(duplicated(taxnames)))
      cat("WARNING: Truncation of taxon names created",
          "identical strings.")
  }
  else {
    xx <- nchar(taxnames)
    diff <- max(xx) - xx + 3
    for (i in 1:ntax) taxnames[i] <- paste(taxnames[i], paste(rep(" ",
                                                                  diff[i]), collapse = ""), sep = "")
  }
  if (!interleave)
    interleave <- nchar
  nbpart <- ceiling(nchar/interleave)
  pt <- matrix(nrow = nbpart, ncol = 2)
  pt[1, ] <- c(1, interleave)
  if (nbpart > 1)
    for (i in 2:(dim(pt)[1])) {
      pt[i, ] <- c(pt[i - 1, 2] + 1, pt[i - 1, 2] + interleave)
      pt[nbpart, 2] <- nchar
    }
  phy <- paste(ntax, nchar)
  for (i in seq(along = pt[, 1])) {
    sm <- as.character(x[, pt[i, 1]:pt[i, 2]])
    if (is.null(dim(sm)))
      sm <- as.matrix(sm, ncol = 1)
    sm <- apply(sm, 1, paste, collapse = "")
    if (i == 1)
      sm <- paste(taxnames, sm)
    if (i < max(seq(along = pt[, 1])))
      sm <- c(sm, "")
    phy <- c(phy, sm)
  }
  if (file == "") {
    cat(phy, sep = "\n")
  }
  else {
    write(phy, file = file)
  }
}



