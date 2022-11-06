#' @title filterGeneTrees
#'
#' @description Filters gene trees using provided filtration and alignment data
#'
#' @param filter.summary summary data file from filterSummary
#'
#' @param alignment.data summary data file from alignmentSummary
#'
#' @param genetree.folder your target folder of gene trees that correspond to the alignments being filtered
#'
#' @param format save format for genetrees
#'
#' @param overwrite if TRUE overwrites file if it exists; FALSE the dataset is skipped
#'
#' @param taxa.remove species that you would like removed from each gene tree
#'
#' @param min.trees mimimum number of trees to keep filtered set
#'
#' @param min.n.samples the minimum number of samples to keep a gene tree
#'
#' @param min.sample.prop the minimum proportion of samples to keep a gene tree
#'
#' @param make.polytomy collapses polytomies in the gene trees
#'
#' @param polytomy.limit the value at which to collapse a node into a polytomy
#'
#' @param remove.node.labels strips trees of node labels if downstream analyses give you trouble (not recommended)
#'
#' @return filters gene trees either in a folder or concatenated set of trees
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

filterGeneTrees = function(filter.summary = NULL,
                           alignment.data = NULL,
                           genetree.folder = NULL,
                           format = c("folder", "concatenated"),
                           overwrite = FALSE,
                           taxa.remove = NULL,
                           min.trees = 5,
                           min.n.samples = 4,
                           min.sample.prop = NULL,
                           make.polytomy = TRUE,
                           polytomy.limit = 0,
                           remove.node.labels = FALSE) {
#
  # filter.summary = filt.summary
  # alignment.data = align.summary
  # genetree.folder = dataset.trees
  # format = "concatenated"
  # overwrite = FALSE
  # taxa.remove = NULL
  # min.trees = 5
  # min.n.samples = 4
  # min.sample.prop = NULL
  # make.polytomy = TRUE
  # polytomy.limit = 10
  # remove.node.labels = FALSE

  if(is.null(alignment.data) == TRUE){ stop("Error: a table of alignment data is needed.") }
  if(is.null(filter.summary) == TRUE){ stop("Error: a filter.summary file is needed.") }
  if(is.null(genetree.folder) == TRUE){ stop("Error: a folder of gene trees in genetree.folder is needed.") }
  if(is.null(format) == TRUE){ stop("Error: an output format (folder or concatenated) is needed.") }
  if(length(format) != 1){ stop("Error: only one output format can be provided.") }
  if(min.n.samples <= 3){ stop("Error: too few samples selected. Must be 4 or greater")}

  #Check if files exist or not
  if (dir.exists(genetree.folder) == F){
    return(paste0("Directory of geen trees could not be found. Exiting."))
  }#end file check

  #Sets up directory for output
  if (format == "folder"){
    if (dir.exists("filtered-genetrees-folders") == F){ dir.create("filtered-genetrees-folders") }
    #Checks for output directory and creates it if not found
    if (overwrite == TRUE){
      if (dir.exists("filtered-genetrees-folders") == T){
        unlink("filtered-genetrees-folders", recursive = T)
        dir.create("filtered-genetrees-folders")
      }#end dir exist
    }#end overwrite
  }#end folder format

  if (format == "concatenated"){
    #Sets up directory for output
    if (dir.exists("filtered-genetrees-concatenated") == F){ dir.create("filtered-genetrees-concatenated") }
    #Checks for output directory and creates it if not found
    if (overwrite == TRUE){
      if (dir.exists("filtered-genetrees-concatenated") == T){
        unlink("filtered-genetrees-concatenated", recursive = T)
        dir.create("filtered-genetrees-concatenated")
      }#end dir exist
    }#end overwrite
  } #end format

  #Read in alignment data and set up
  if (length(alignment.data) == 1) {
    alignment.stats = data.table::fread(alignment.data, header = T)
  } else {
    alignment.stats = alignment.data
  }

  #Gets list of gene trees
  gene.trees = list.files(genetree.folder)

  if (length(gene.trees) == 0){ stop("Error: no gene trees found.") }

  #Goes through each filter and applies the filter from the raw data
  for (y in 1:nrow(filter.summary)){

    #Gets filter subset
    temp.filter = filter.summary[y,]

    if (temp.filter$no_trees < min.trees){
      print(paste0(temp.filter$filter_file, " does not have enough trees. skipping..."))
      next
    }#end if

    #overwriting
    if (overwrite == FALSE){
        if (file.exists(paste0("filtered-genetrees-concatenated/", temp.filter$filter_file, "_genetrees.tre")) == TRUE){
          print(paste0(temp.filter$filter_file, "_genetrees.tre already exists and overwrite = FALSE. skipping..."))
          next
        }#end file exists
    }#end overwrite if

    #Applies filters
    filt.data = alignment.stats[alignment.stats$alignment_length >= temp.filter$filter_length,]
    filt.data = filt.data[filt.data$proportion_samples >= temp.filter$filter_sample,]
    filt.data = filt.data[filt.data$proportion_pis >= temp.filter$filter_prop_pis,]
    filt.data = filt.data[filt.data$count_pis >= temp.filter$filter_count_pis,]
    filt.data = filt.data[filt.data$number_samples >= min.n.samples,]

    #skips minimum number of trees for dataset
    if (nrow(filt.data) < min.trees){ next }

    #Obtains list of markers and associated gene trees
    marker.list = filt.data$file
    marker.list = gsub("\\..*", "", marker.list)
    tree.files = gene.trees[gsub("\\..*", "", gene.trees) %in% marker.list]

    #Save foldername as dataset name
    if (length(format[format == "folder"]) == 1){
      dir.create(paste0("filtered-genetrees-folders/", temp.filter$filter_file))
    } #end folder if

    #Obtains total number of samples
    max.taxa = c()
    if (is.null(min.sample.prop) != TRUE & length(max.taxa) != 0){
      total.taxa = c()
      for (x in 1:length(tree.files)){
        temp.tree = ape::read.tree(paste0(genetree.folder, "/", tree.files[x]))
        total.taxa = append(total.taxa, length(temp.tree$tip.label))
      }#end x loop
      max.taxa = max(total.taxa)
    }#end if for sample prop

    #Loops through and places trees
    tree.counter = 0
    for (x in 1:length(tree.files)){
      #read in tree file
      temp.tree = ape::read.tree(paste0(genetree.folder, "/", tree.files[x]))

      #Removes taxa if it needs to
      if (is.null(taxa.remove) != TRUE){
        temp.tree = ape::drop.tip(temp.tree, taxa.remove)
      }

      #Skips if less than 4 taxa
      if (length(temp.tree$tip.label) < min.n.samples){
        #print(paste0(tree.files[x], " skipped, less than ", min.n.samples, " samples."))
        next
      }#end min sample check

      if (is.null(min.sample.prop) != TRUE){
        #Skips if less than the desire sampling proportion
        if (length(temp.tree$tip.label)/max.taxa < min.sample.prop){
          print(paste0(gene.trees[x], " skipped, less than ",
                       min.sample.prop, " proportion samples."))
          next
        }#end min sample check
      }

      #Does the polytomy check
      if (make.polytomy == TRUE){
        #Checks for node labels
        if (length(temp.tree$node.label) == 0) {
          new.tree = temp.tree
        } else {
          #Find and collapse nodes with bl close to 0 from above
          temp.tree$node.label[temp.tree$node.label == ""] = "100"
          new.tree = AstralPlane::makePolytomy(tree = temp.tree, polytomy.limit = polytomy.limit)
        }#end else
      }#make polytomy end if

      if (class(new.tree) == "character"){ next }

      if (length(format[format == "concatenated"]) == 1){
        #writes tree to a single file
        ape::write.tree(new.tree, file = paste0("filtered-genetrees-concatenated/",
                                                temp.filter$filter_file, "_genetrees.tre"), append = T)
      }#end if

      #Save foldername as dataset name
      if (length(format[format == "folder"]) == 1){
        #Writes to a separate file
        ape::write.tree(new.tree, file = paste0("filtered-genetree-folders/",
                                      temp.filter$filter_file, "/",
                                      tree.files[x]))
      } #end folder if
      tree.counter = tree.counter + 1
    }#end x loop

    print(paste0(temp.filter$filter_file, " complete! Wrote ", tree.counter,
                 " gene trees to file: ", temp.filter$filter_file, "_genetrees.tre"))
  }#end y loop

}#end function
