#' @title summarizeAlignments
#'
#' @description Calculates per-alignment summary statistics (number of samples,
#'   proportion of maximum sampling, alignment length, count and proportion of
#'   parsimony informative sites, and missing data) across a folder of alignments.
#'   Results are returned as a data.table and optionally saved as a CSV.
#'
#' @param alignment.path path to a folder of alignment files
#'
#' @param file.export base file name (without extension) for saving the CSV
#'   summary; if NULL the summary is only returned in R
#'
#' @param overwrite if TRUE overwrites an existing CSV; if FALSE and the file
#'   exists, the existing file is loaded and returned
#'
#' @param dataset.name a label for the dataset (e.g. "exons", "UCEs") stored
#'   in the "dataset" column of the output table
#'
#' @param alignment.format format of the input alignments: "phylip" or "nexus"
#'
#' @return a data.table with one row per alignment containing: dataset, file,
#'   number_samples, proportion_samples, alignment_length, count_pis,
#'   proportion_pis, count_missing_bp, proportion_missing_bp
#'
#' @examples
#'
#' align.summary = summarizeAlignments(alignment.path = "path/to/alignments",
#'                                     dataset.name = "exons",
#'                                     file.export = "alignment_summary",
#'                                     alignment.format = "phylip")
#'
#' @export


summarizeAlignments = function(alignment.path = NULL,
                               file.export = NULL,
                               overwrite = FALSE,
                               dataset.name = NULL,
                               alignment.format = c("phylip", "nexus")) {

  #alignment.path = "/Volumes/LaCie/Anax/data-analysis/alignments/untrimmed_all-markers"
  #dataset.name = "test"
  #file.export = "test.csv"
  #alignment.format = "phylip"
  #overwrite = FALSE

  if(is.null(alignment.path) == TRUE){ stop("Error: no alignment path provided.") }
  if(is.null(dataset.name) == TRUE){ stop("Error: a dataset name is needed.") }

  #Check if files exist or not
  if (dir.exists(alignment.path) == F){
    return(paste0("Directory of alignments could not be found. Exiting."))
  }#end file check

  #Overwrite checker
  if (overwrite == TRUE){
    if (file.exists(paste0(file.export, ".csv")) == T){
      #Checks for output directory and creates it if not found
      system(paste0("rm ", file.export, ".csv"))
    }#end file exists
  } else {
    if (file.exists(paste0(file.export, ".csv")) == T){
      print(paste0("File exists for ", file.export, " and overwrite = FALSE. Exiting."))
      save.data = read.csv(paste0(file.export, ".csv"))
      return(save.data)
    }#end file check
  }#end else

  #Gets list of alignments from path
  align.names = list.files(alignment.path)

  #Collects the super cool data
  header.data = c("dataset", "file", "number_samples", "proportion_samples", "alignment_length",
                  "count_pis", "proportion_pis", "count_missing_bp", "proportion_missing_bp")
  #Sets up data collection data.frame
  collect.data = data.table::data.table(matrix(as.numeric(0),
                                   nrow = length(align.names),
                                   ncol = length(header.data)))
  data.table::setnames(collect.data, header.data)
  collect.data[, file:=as.character(file)]
  collect.data[, dataset:=as.character(dataset)]

  #Loops through each alignment to gather statistics
  for (x in 1:length(align.names)){
    #Reads in alignment

    if (alignment.format == "phylip"){
      align = ape::read.dna(paste0(alignment.path, "/", align.names[x]),
                       format = "sequential")
    }
    if (alignment.format == "nexus"){
      align = ape::read.nexus.data(paste0(alignment.path, "/", align.names[x]))
      align = ape::as.DNAbin(matrix(unlist(align), ncol = length(align[[1]]), byrow = TRUE))
    }

    #Collect data
    data.table::set(collect.data, i = as.integer(x), j = match("dataset", header.data), value = dataset.name )
    data.table::set(collect.data, i = as.integer(x), j = match("file", header.data), value = align.names[x] )
    #Sample data
    data.table::set(collect.data, i = as.integer(x), j = match("number_samples", header.data), value = nrow(align) )
    #Length data
    data.table::set(collect.data, i = as.integer(x), j = match("alignment_length", header.data), value = ncol(align) )

    count.pis = PhyloConfigR::informativeSites(align, count = T, ambiguities = T)
    prop.pis = round(count.pis/ncol(align),3)
    data.table::set(collect.data, i = as.integer(x), j = match("count_pis", header.data), value = count.pis)
    data.table::set(collect.data, i = as.integer(x), j = match("proportion_pis", header.data), value = prop.pis)

    #Removes samples that too short individually
    len.temp = as.character(as.list(align))
    len.loci = lapply(len.temp, function (x) x[x != "-"])
    len.loci = lapply(len.loci, function (x) x[x != "n"])
    len.loci = lapply(len.loci, function (x) x[x != "?"])
    spp.len = unlist(lapply(len.loci, function (x) length(x)))
    miss.total = (max(spp.len) - spp.len)
    miss.prop = round(sum(miss.total)/(max(spp.len)*nrow(align)), 3)

    #Get missing bp data
    data.table::set(collect.data, i = as.integer(x), j = match("count_missing_bp", header.data), value =  sum(miss.total) )
    data.table::set(collect.data, i = as.integer(x), j = match("proportion_missing_bp", header.data), value = miss.prop )

  } #x loop

  save.data = collect.data[collect.data$file != 0,]
  save.data[, proportion_samples:=round(number_samples/max(number_samples), 3)]

  if (is.null(file.export) != TRUE){
    write.csv(save.data, file = paste0(file.export, ".csv"), row.names = F)
  }#end if

  return(save.data)

}#End function summarizeAlignments

