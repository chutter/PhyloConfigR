#' @title filterConcordance
#'
#' @description Function for obtaining filtered concordance factor data
#'
#' @param input.dir directory of concordance factor data generated from the filtered datasets
#'
#' @param clade.list a named list of clades of interest to test for concordance factors
#'
#' @param outgroups outgroups to root the tree
#'
#' @return data.frame with concordance factor data for each filtered replicate
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


# Filter concordance function
filterConcordance = function(input.dir = NULL,
                             clade.list = NULL,
                             outgroups = NULL) {

  #Debug
 # input.dir = "concordance-factors"
#  clade.list = NULL
 # outgroups  = outgroup.taxa

  if (is.null(input.dir) == TRUE){ stop("Error: No input directory provided.") }
  if (is.null(outgroups) == TRUE){ stop("Error: No outgroups provided.") }

  cf.files = list.files(input.dir)
  cf.files = unique(gsub("\\..*", "", cf.files))

  save.data = c()
  for (x in 1:length(cf.files)){

    #Reads in the plane data
    plane.data = AstralPlane::createAstralPlaneCF(cf.file.name = paste0(input.dir, "/", cf.files[x]),
                                                  outgroups = outgroups,
                                                  tip.length = 1)

    #Merges the CF and node data
    merge.data = merge(plane.data@concordanceFactorData,
                       plane.data@nodeData,
                       by = "node")

    #Sets up new columns
    merge.data = cbind(dataset = gsub(".*/", "", cf.files[x]), merge.data)

    #only marks the clades if given a clade list
    if (is.null(clade.list) == FALSE){
      #Adds the clade columns
      merge.data = cbind(monophyletic = NA, merge.data)
      merge.data = cbind(clade = NA, merge.data)
      for (y in 1:length(clade.list)){

        #Overall target tree
        mrca.node = ape::getMRCA(plane.data@phylo, clade.list[[y]])

        if (ape::is.monophyletic(plane.data@phylo,  clade.list[[y]]) == T){
          mono.clade = TRUE
        } else{
          mono.clade = FALSE
        }#end else

        #What is the goal here?
        merge.data[merge.data$node %in% mrca.node,]$clade = names(clade.list)[y]
        merge.data[merge.data$node %in% mrca.node,]$monophyletic = mono.clade

        #Obtain other close nodes
        temp.a = plane.data@edgeData
        node.data1 = temp.a[temp.a$node1 %in% mrca.node,]
        node.data2 = temp.a[temp.a$node2 %in% mrca.node,]
        node.data = rbind(node.data1, node.data2)
        more.nodes = unique(c(node.data$node1, node.data$node2))
        merge.data[merge.data$node %in% more.nodes,]$clade = names(clade.list)[y]
      }#y loop
    }#End if statement

    save.data = rbind(save.data, merge.data)

  }#end x loop

  #Returns all the data
  save.data$node = as.numeric(save.data$node)
  return(save.data)

}#end function

