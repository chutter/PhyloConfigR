#' @title filterStats
#'
#' @description Function for gather filtration statistics, used internally by filterSummary
#'
#' @param data alignment summmary data from alignmentSummary function
#'
#' @param filter.name Name of the filter being passed. Only one can be used at a time.
#'
#' @param filter.values the filter values for the given scheme
#'
#' @param align.dataset select three colors to plot your pie.plot
#'
#' @return Alignment statistics for the specific filteration scheme provided
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

filterStats = function(data = NULL,
                       filter.name = c("alignment_length", "count_pis", "proportion_pis", "proportion_samples"),
                       filter.values = NULL,
                       align.dataset = NULL) {


  # data = alignment.stats
  # filter.name = "proportion_samples"
  # filter.values = sample.filters
  # align.dataset = "guibemantis"

  #Parameter checks
  if(is.null(data) == TRUE){ stop("Error: no filter data provided.") }
  if (length(filter.name) != 1){ stop("Error: please provide only 1 filter") }

  #Gets filtered dataset stats together
  header.data = c("align_dataset", "filter_file", "filter_name",  "no_trees",
                  "mean_length", "mean_sample", "mean_prop_pis", "mean_count_pis",
                  "filter_length", "filter_sample", "filter_prop_pis", "filter_count_pis")

  collect.data = data.table(matrix(as.numeric(0),
                                   nrow = length(filter.values),
                                   ncol = length(header.data)))
  setnames(collect.data, header.data)
  collect.data[, filter_file:=as.character(filter_file)]
  collect.data[, filter_name:=as.character(filter_name)]
  collect.data[, align_dataset:=as.character(align_dataset)]

  for (x in 1:length(filter.values)){

    #Apply size filter
    if(filter.name == "alignment_length") {
      filt.data = data[data$alignment_length >= filter.values[x],]
      sample.value = 0
      count.value = 0
      length.value = filter.values[x]
      prop.value = 0
      #Makes file name
      file.name = paste0(align.dataset, "_", filter.name,
                         "_filters-S0", "-L", filter.values[x],
                         "-CPIS0", "-PPIS0")
    }

    #Apply size filter
    if(filter.name == "count_pis") {
      filt.data = data[data$count_pis >= filter.values[x],]
      sample.value = 0
      count.value = filter.values[x]
      length.value = 0
      prop.value = 0
      file.name = paste0(align.dataset, "_", filter.name,
                         "_filters-S0", "-L0",
                         "-CPIS", filter.values[x], "-PPIS0")
    }

    #Apply size filter
    if(filter.name == "proportion_pis") {
      filt.data = data[data$proportion_pis >= filter.values[x],]
      sample.value = 0
      count.value = 0
      length.value = 0
      prop.value = filter.values[x]
      file.name = paste0(align.dataset, "_", filter.name,
                         "_filters-S0", "-L0",
                         "-CPIS0", "-PPIS", filter.values[x]*100)
    }

    #Apply size filter
    if(filter.name == "proportion_samples") {
      filt.data = data[data$proportion_samples >= filter.values[x],]
      sample.value = filter.values[x]
      count.value = 0
      length.value = 0
      prop.value = 0
      file.name = paste0(align.dataset, "_", filter.name,
                         "_filters-S", filter.values[x] * 100, "-L0",
                         "-CPIS0", "-PPIS0")
    }

    #Collect data
    set(collect.data, i = as.integer(x), j = match("align_dataset", header.data), value = align.dataset )
    set(collect.data, i = as.integer(x), j = match("filter_file", header.data), value = file.name )
    set(collect.data, i = as.integer(x), j = match("filter_name", header.data), value = filter.name )
    set(collect.data, i = as.integer(x), j = match("no_trees", header.data), value = nrow(filt.data) )

    if (nrow(filt.data) == 0){ next }

    set(collect.data, i = as.integer(x), j = match("mean_length", header.data), value = mean(filt.data$alignment_length) )
    set(collect.data, i = as.integer(x), j = match("mean_sample", header.data), value = mean(filt.data$proportion_samples) )
    set(collect.data, i = as.integer(x), j = match("mean_prop_pis", header.data), value = mean(filt.data$proportion_pis) )
    set(collect.data, i = as.integer(x), j = match("mean_count_pis", header.data), value = mean(filt.data$count_pis) )

    set(collect.data, i = as.integer(x), j = match("filter_length", header.data), value = length.value )
    set(collect.data, i = as.integer(x), j = match("filter_sample", header.data), value = sample.value )
    set(collect.data, i = as.integer(x), j = match("filter_prop_pis", header.data), value = prop.value )
    set(collect.data, i = as.integer(x), j = match("filter_count_pis", header.data), value = count.value )

  }#end x loop

  return(collect.data)

} #end filterStats

