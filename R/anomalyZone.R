#' @title anomalyZone
#'
#' @description Applies the anomaly zone calculation across the entire tree (Degnan & Rosenberg 2006)
#'
#' @param tree phylogenetic tree read into R in phylo format
#'
#' @param outgroups outgroups to root your tree
#'
#' @param print.node prints the anomaly zone nodes if found
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

anomalyZone = function(tree = NULL,
                       outgroups = NULL,
                       print.node = TRUE) {

  #Parameter checks
  if (is.null(tree) == TRUE){ stop("Error: No phylogenetic tree provided.") }
  if (is.null(outgroups) == TRUE){ stop("Error: No outgroups provided.") }

  #If imported as a file path
  if (class(tree) == "character"){
    spp.tree = ape::read.tree(tree)
    spp.tree = ape::unroot(spp.tree)
  }

  #If imported as a tree
  if (class(tree) == "phylo"){
    spp.tree = ape::unroot(tree)
  }

  if (ape::is.monophyletic(spp.tree, outgroups) == T){
    spp.tree = ape::root(spp.tree, outgroups, resolve.root = T)
  } else{ spp.tree = ape::root(spp.tree, outgroups[1], resolve.root = T) }

  edge.node = AstralPlane::edgeLengthTable(spp.tree, tips = FALSE)

  #Gets the new tree fiels it made
  header.data = c("parent_node", "child_node", "parent_branch", "child_branch",
                  "x_branch_length", "y_branch_length", "a_s", "anomaly_zone")

  collect.data = data.table(matrix(as.numeric(0),
                                   nrow = nrow(edge.node),
                                   ncol = length(header.data)))
  setnames(collect.data, header.data)

  row.index = as.integer(1)
  for (x in 1:nrow(edge.node)){

    #Gets parent and child node
    parent.node = edge.node$node1[x]
    child.node = edge.node$node2[x]

    #Gets parent and child edge number
    par.branch = edge.node[edge.node$node2 == parent.node,]$edge
    child.branch = edge.node[edge.node$node2 == child.node,]$edge

    #child node branch length
    x.node = spp.tree$edge.length[par.branch]
    y.node = spp.tree$edge.length[child.branch]

    #function a(x) from Degnan and Rosenberg (2006)
    ax = log(2.0/3+((3*exp(2*x.node)-2)/(18*(exp(3*x.node)-exp(2*x.node)))))

    #if there no is value, if its NA, or if its Infinite, it skips it
    if (length(ax) == 0){ next }
    if (is.na(ax) == T){ next }
    if (ax == Inf){ next }

    #No anomaly zone, unless lambda2 is less than a(x)
    anom = 0
    if (y.node < ax){ anom = 1 } else {anom = 0 }

    if (anom != 0){
      if (print.node == TRUE){
        print(paste0("anomaly zone found at node ", parent.node))
      }#end if 2
    }#end if 1

    set(collect.data, i = as.integer(x), j = match("parent_node", header.data), value = parent.node)
    set(collect.data, i = as.integer(x), j = match("child_node", header.data), value = child.node)
    set(collect.data, i = as.integer(x), j = match("parent_branch", header.data), value = par.branch )
    set(collect.data, i = as.integer(x), j = match("child_branch", header.data), value = child.branch )

    set(collect.data, i = as.integer(x), j = match("x_branch_length", header.data), value = x.node )
    set(collect.data, i = as.integer(x), j = match("y_branch_length", header.data), value = y.node )
    set(collect.data, i = as.integer(x), j = match("a_s", header.data), value = ax )
    set(collect.data, i = as.integer(x), j = match("anomaly_zone", header.data), value = anom )
  }#end x loop

  save.data = collect.data[collect.data$parent_node != 0,]
  return(save.data)

} #end function

