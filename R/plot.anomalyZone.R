#' @title plot.anomalyZone
#'
#' @description Function for plotting data collected for the anomaly zone
#'
#' @param tree Tree that you wish to plot
#'
#' @param data data.table from the anomalyZone function
#'
#' @param outgroups provide outgroups for rooting plotted tree
#'
#' @param save.file if you wish to save to file, put file name. Saves as PDF
#'
#' @param tip.label.size size of the tip labels, passed to cex in plotting function
#'
#' @param node.label.size size of the node label circles, passed to cex in plotting function
#'
#' @param edge.width size of the branch edges, passed to edge.width in plot.phylo
#'
#' @return plots the phylogenetic tree and selected data associated with an AstralPlane object. Can optionally be saved to file as a PDF by giving save.file a file name.
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

### Plots the anomaly zone on the tree
plot.anomalyZone = function(tree = NULL,
                            data = NULL,
                            outgroups = NULL,
                            save.file = NULL,
                            tip.label.size = 0.5,
                            node.label.size = 1,
                            edge.width = 3) {

  #Checks for tree
  if (is.null(tree) == TRUE){ stop("Error: No phylogenetic tree provided.") }
  if (is.null(data) == TRUE){ stop("Error: No anomaly zone data provided.") }
  if (is.null(outgroups) == TRUE){ stop("Error: No outgroups provided.") }

  #If imported as a file path
  if (class(tree) == "character"){
    tree = ape::read.tree(tree)
    tree = ape::unroot(tree)
  }

  #If imported as a tree
  if (class(tree) == "phylo"){
    tree = ape::unroot(tree)
  }

  if (ape::is.monophyletic(tree, outgroups) == T){
    tree = ape::root(tree, outgroups, resolve.root = T)
  } else{ tree = ape::root(tree, outgroups[1], resolve.root = T) }

  #Grabs the node numbers
  node.no = as.numeric(length(tree$tip.label)+1:tree$Nnode)
  tree$edge.length[is.na(tree$edge.length) == T] = 1

  #Creates dataset for nodes to be plotted
  #########################################
  yes.states = data[data$anomaly_zone == 1,]
  yes.states = c(yes.states$parent_node, yes.states$child_node)
  yes.states = yes.states[duplicated(yes.states) != T]
  plot.states = rep(1, length(yes.states))
  names(plot.states) = as.character(yes.states)

  missing = as.character(node.no)[!as.character(node.no) %in% names(plot.states)]
  miss.state = rep(0, length(missing))
  names(miss.state) = missing
  plot.states = append(plot.states, miss.state)
  plot.states = plot.states[order(names(plot.states))]

  node.colors = plot.states
  node.colors[node.colors == 0] = "#7BC143"
  node.colors[node.colors == 1] = "#DE3293"

  #Creates dataset for branches to be plotted
  #########################################
  branch.no = AstralPlane::edgeLengthTable(tree, tips = TRUE)$edge
  yes.states = data[data$anomaly_zone == 1,]
  yes.states = c(yes.states$parent_branch, yes.states$child_branch)
  yes.states = yes.states[duplicated(yes.states) != T]
  b.states = rep(1, length(yes.states))
  names(b.states) = as.character(yes.states)

  missing = as.character(branch.no)[!as.character(branch.no) %in% names(b.states)]
  miss.state = rep(0, length(missing))
  names(miss.state) = missing
  b.states = append(b.states, miss.state)
  b.states = b.states[order(match(names(b.states), branch.no))]

  branch.colors = b.states
  branch.colors[branch.colors == 0] = "#7BC143"
  branch.colors[branch.colors == 1] = "#DE3293"

  if(is.null(save.file) != T){ pdf(file = save.file, width = 10, height = 8) }

  #Plots the stuff
  ape::plot.phylo(tree, cex = tip.label.size, edge.color = branch.colors, edge.width = edge.width)
  labels = rep("", length(node.colors))
  ape::nodelabels(labels, frame = "circle", cex = node.label.size, bg = node.colors)

  if(is.null(save.file) != T){  dev.off() }

}#end function

