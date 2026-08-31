#' @title informativeSites
#'
#' @description Calculates the number or proportion of parsimony informative sites in an alignment
#'
#' @param alignment alignment in ape DNABin or a matrix format
#'
#' @param count Whethe to return the count of parsimoney informative sites (TRUE) or the proportion (FALSE)
#'
#' @param ambiguities Whether to consider ambiguities (TRUE) or not (FALSE)
#'
#' @return if count = TRUE, returns the integer count of parsimony informative
#'   sites; if count = FALSE, returns the proportion (0-1)
#'
#' @examples
#'
#' align = ape::read.dna("path/to/alignment.phy", format = "sequential")
#' pis_count = informativeSites(alignment = align, count = TRUE, ambiguities = TRUE)
#'
#' @export

#Calculates informative sites
informativeSites = function(alignment = NULL,
                            count = TRUE,
                            ambiguities = TRUE) {

  #Helper function to use with apply
  column.pars = function(x) {
    x = table(x)
    x = x[x > 1]
    if (length(x[!names(x) %in% n]) > 1){ return(TRUE) } else { return(FALSE) }
  }#end function

  #characters to exclude
  n = c("-", "?", "n")
  if (ambiguities == FALSE){
    n = append(n, c( "r", "y", "k", "m", "s", "w", "b", "d", "h", "v")) }

  #alignment length
  x.len = dim(alignment)[2]
  #goes through each column and sees if they are different
  alignment = as.character(alignment)
  col.pis = apply(alignment, 2, column.pars)
  out = length(col.pis[col.pis == TRUE])

  if (count != TRUE){ out = round(out/x.len, digits = 3) }
  return(out)

}#end informative sites function
