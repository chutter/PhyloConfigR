#' @title heterozygousSites
#'
#' @description Counts columns in an alignment that contain heterozygous (IUPAC
#'   ambiguity) characters, treating each such column as a heterozygous site.
#'   Columns with only missing data markers are not counted.
#'
#' @param alignment alignment in ape DNAbin or matrix format
#'
#' @param count if TRUE (default) returns the integer count of heterozygous
#'   sites; if FALSE returns the proportion relative to alignment length
#'
#' @param ambiguities if TRUE (default) IUPAC ambiguity codes (R, Y, K, M, S,
#'   W, B, D, H, V) are treated as heterozygous; if FALSE they are excluded
#'
#' @return integer count or numeric proportion of heterozygous sites
#'
#' @examples
#'
#' align = ape::read.dna("path/to/alignment.phy", format = "sequential")
#' het_count = heterozygousSites(alignment = align, count = TRUE, ambiguities = TRUE)
#'
#' @export

#Calculates informative sites
heterozygousSites = function(alignment = NULL,
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
