#' @title filterAlignments
#'
#' @description Function for filtering alignments based on filter summary data
#'
#' @param filter.summary summary data file from filterSummary
#'
#' @param alignment.data summary data file from alignmentSummary
#'
#' @param alignment.folder folder of alignments to be filtered
#'
#' @param format save format for alignments
#'
#' @param min.alignments minimum number of alignments for a filtered set of data
#'
#' @param overwrite if TRUE overwrites file if it exists; FALSE the dataset is skipped
#'
#' @return filters the alignments
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

##### Create filtered alignment datasets
filterAlignments = function(filter.summary = NULL,
                            alignment.data = NULL,
                            alignment.folder = NULL,
                            format = c("folder", "concatenated"),
                            min.alignments = 5,
                            min.n.samples = 4,
                            overwrite = FALSE ) {

  # filter.summary = filt.summary
  # alignment.data = align.summary
  # alignment.folder = align.dir
  # format = "concatenated"
  # min.alignments = 5
  # min.n.samples = 5
  # overwrite = FALSE

  #Parameter checks
  if(is.null(filter.summary) == TRUE){ stop("Error: a filter.summary file is needed.") }
  if(is.null(alignment.data) == TRUE){ stop("Error: an alignment.data file is needed.") }
  if(is.null(alignment.folder) == TRUE){ stop("Error: a folder of alignments is needed.") }
  if(is.null(format) == TRUE){ stop("Error: an output format (folder or concatenated) is needed.") }
  if(min.n.samples <= 3){ stop("Error: too few samples selected. Must be 4 or greater")}

  #Check if files exist or not
  if (dir.exists(alignment.folder) == F){
    return(paste0("Directory of alignments could not be found. Exiting."))
  }#end file check

  #Sets up directory for output
  #Deletes folder if not wanting to save
  if (length(format[format == "folder"]) == 1){
    #Checks if it exists
    if (dir.exists("filtered-alignments") == F){ dir.create("filtered-alignments") }
    #Checks for output directory and creates it if not found
    if (overwrite == TRUE){
      if (dir.exists("filtered-alignments") == T){
        unlink("filtered-alignments", recursive = T)
        dir.create("filtered-alignments")
      }#end dir exist
    }#end overwrite
  }#end folder if

  #Sets up directory for output
  if (dir.exists("filtered-alignments-concatenated") == F){ dir.create("filtered-alignments-concatenated") }
  #Checks for output directory and creates it if not found
  if (overwrite == TRUE){
    if (dir.exists("filtered-alignments-concatenated") == T){
      unlink("filtered-alignments-concatenated", recursive = T)
      dir.create("filtered-alignments-concatenated")
    }#end dir exist
  }#end overwrite

  #Read in alignment data and set up
  if (length(alignment.data) == 1) {
    alignment.stats = data.table::fread(alignment.data, header = T)
  } else {
    alignment.stats = alignment.data
  }
  #Gets list of alignment files
  alignment.files = list.files(alignment.folder)

  #Goes through each filter and applies the filter from the raw data
  for (x in 1:nrow(filter.summary)){

    #Gets filter subset
    temp.filter = filter.summary[x,]
    #Checks for too few trees
    if (temp.filter$no_trees < min.alignments){
      print(paste0(temp.filter$filter_file, " does not have enough trees. skipping..."))
      next
    }#end if

    #overwriting
    if (overwrite == FALSE){
      if (file.exists(paste0("filtered-alignments-concatenated/", temp.filter$filter_file, ".phy")) == TRUE){
        print(paste0(temp.filter$filter_file, ".phy already exists and overwrite = FALSE. skipping"))
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
    if (nrow(filt.data) < min.alignments){ next }

    #Obtains list of markers and associated gene trees
    marker.list = filt.data$file
    marker.list = gsub("\\..*", "", marker.list)
    align.files = alignment.files[gsub("\\..*", "", alignment.files) %in% marker.list]

    #Deletes folder if not wanting to save
    if (length(format[format == "concatenated"]) == 1){
      concatenateAlignments(alignment.path = alignment.folder,
                            alignment.names = align.files,
                            file.name = temp.filter$filter_file,
                            output.dir = "filtered-alignments-concatenated",
                            partition.format = c("raxml"))
    } #end if

    #Deletes folder if not wanting to save
    if (length(format[format == "folder"]) == 1){
      #Save foldername as dataset name
      dir.create(paste0("filtered-alignments/", temp.filter$filter_file))

      #Loops through and places trees
      for (y in 1:length(align.files)){
        #Save foldername as dataset name
          system(paste0("cp ", alignment.folder, "/", align.files[y],
                        " filtered-alignments/", temp.filter$filter_file))
      }#end y loop for each alignment within filtered dataset
    }

    print(paste0(temp.filter$filter_file, " complete!"))


  }#end x loop for each filtered dataset

}#end function



