#' @title hyphy.plotSLAC
#'
#' @description Plots SLAC omega (dN/dS) values on a species tree using a
#'   continuous color mapping. Tip-level mean omega values are mapped onto
#'   the phylogeny via phytools::contMap.
#'
#' @param slac.spreadsheet path to the CSV file produced by hyphy.Parser with
#'   hyphy.analysis = "SLAC" (the "_SLAC_by-branch.csv" file)
#'
#' @param species.tree path to a species tree file readable by ape::read.tree
#'
#' @param outgroups character vector of outgroup taxon names used to root the tree
#'
#' @param log.transform if TRUE (default), omega values are log10-transformed
#'   before plotting
#'
#' @param lwd line width for tree branches passed to phytools::contMap
#'
#' @param res resolution of the color gradient passed to phytools::contMap
#'
#' @return a contMap object (returned invisibly); the plot is drawn as a side effect
#'
#' @examples
#'
#' hyphy.plotSLAC(slac.spreadsheet = "slac_results_SLAC_by-branch.csv",
#'                species.tree = "my_species_tree.tre",
#'                outgroups = c("Outgroup_species"),
#'                log.transform = TRUE)
#'
#' @export

hyphy.plotSLAC = function(slac.spreadsheet = NULL,
                          species.tree = NULL,
                          outgroups = NULL,
                          log.transform = TRUE,
                          lwd = 4,
                          res = 100) {

  if (is.null(slac.spreadsheet) == T){ stop("A SLAC spreadsheet is needed.") }
  if (is.null(species.tree) == T){ stop("A species tree file path is needed.") }
  if (is.null(outgroups) == T){ stop("Outgroups are needed to root the tree.") }

  if (!file.exists(slac.spreadsheet)){ stop("SLAC spreadsheet file could not be found.") }
  if (!file.exists(species.tree)){ stop("Species tree file could not be found.") }

  slac.results = data.table::fread(slac.spreadsheet)
  slac.results = slac.results[!is.na(slac.results$omega),]
  slac.results = slac.results[!is.infinite(slac.results$omega),]

  mean.slac = aggregate(x = slac.results, by = list(slac.results$sample), FUN = mean)
  sample.omega = mean.slac$omega
  names(sample.omega) = mean.slac$Group.1

  slac.tree = ape::read.tree(species.tree)
  slac.tree = ape::root(slac.tree, outgroup = outgroups, resolve.root = TRUE)
  slac.tree$tip.label = gsub("-", "_", slac.tree$tip.label)
  sample.omega = sample.omega[names(sample.omega) %in% slac.tree$tip.label]

  if (log.transform == TRUE){
    plot.vals = log10(sample.omega)
  } else {
    plot.vals = sample.omega
  }

  cm = phytools::contMap(slac.tree, plot.vals, res = res, fsize = NULL,
                         ftype = NULL, lwd = lwd, legend = NULL,
                         lims = NULL, sig = 3, type = "phylogram",
                         direction = "rightwards")

  return(invisible(cm))

}#end function
