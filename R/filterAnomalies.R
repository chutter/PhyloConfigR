#' @title filterAnomalies
#'
#' @description Applies the anomaly zone calculation across all filtered datasets and combines with filter summary data for downstream analyses
#'
#' @param astral.directory directory of filtered astral results
#'
#' @param outgroups outgroups to root your tree
#'
#' @param filter.data your master filtered dataset summary stats
#'
#' @return a data.table of anomaly zone data calculated for all nodes in your tree
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

filterAnomalies = function(astral.directory = NULL,
                           outgroups = NULL,
                           filter.data = NULL){


  #astral.directory = "filtered-Astral"
  #outgroups = outgroup.taxa
  #filter.data = filt.summary

    #Input checks
  if (is.null(astral.directory) == TRUE){ stop("Error: No directory to filtered astral datasets provided.") }
  if (is.null(outgroups) == TRUE){ stop("Error: Outgroups were not provided.") }
  if (is.null(filter.data) == TRUE){ stop("Error: No filter summary data provided.") }

  #Load in astral data
  astral.files = list.files(astral.directory)

  #Go through each filtered dataset, calculate anomaly zone, and merge AZ data with filtered data
  all.data = data.table()
  for (x in 1:length(astral.files)){
    #Read in tree

    filt.tree = AstralPlane::readAstral(astral.tree = paste0(astral.directory, "/", astral.files[x]),
                                        outgroups = outgroups,
                                        tip.length = 1)

    #anomaly zone calculation for this filtered dataset
    anom.data = anomalyZone(tree = filt.tree,
                            outgroups = outgroups,
                            print.node = FALSE)

    #Find the matching filtered data and combine together
    filt.data = data.table(filter.data[filter.data$filter_file == gsub("_astral.tre", "", astral.files[x]),])
    combined.data = cbind(filt.data, anom.data)

    #Combine together
    all.data = rbind(all.data, combined.data)

  }#end x loop

  return(all.data)

}#end function

