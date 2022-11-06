#' @title concatenateAlignments
#'
#' @description Function for concatenating a large number of alignments
#'
#' @param alignment.folder folder that contains aligmnents to be concatenated
#'
#' @param file.name output file name
#'
#' @param output.dir directory to save file to
#'
#' @param partition.format partition file format. Can select both.
#'
#' @return saves to file concatenated alignments and partition files delimiting the coordinates of each indidividual marker
#'
#' @examples
#'
#' your.tree = ape::read.tree(file = "file-path-to-tree.tre")
#' astral.data = astralPlane(astral.tree = your.tree,
#'                           outgroups = c("species_one", "species_two"),
#'                           tip.length = 1)
#'
#' @export


concatenateAlignments = function(alignment.path = NULL,
                                 alignment.names = NULL,
                                 file.name = NULL,
                                 output.dir = NULL,
                                 partition.format = c("raxml", "table", "none")) {

  # alignment.path = alignment.folder
  # alignment.names = align.files
  # file.name = temp.filter$filter_file
  # output.dir = "filtered-alignments-concatenated"
  # partition.format = "raxml"

  #Gets list of alignments
  alignment.files = list.files(alignment.path, full.names = T)
  alignment.files = alignment.files[gsub(".*\\/", "", alignment.files) %in% alignment.names]

  if (length(alignment.files) != length(alignment.names)){ stop("Not all alignment names provided are in the alignment folder.")}

  align.list = lapply(alignment.files, function (x) fread(x, header = T))
  sample.names = unique(unlist(lapply(align.list, function (x) x[,1])))
  sample.names = sample.names[order(sample.names)]

  #sets up data collection for proportion
  feat.headers = c("Marker", "Length", "Start", "End")
  feature.data = data.table(matrix(as.numeric(0),
                                   nrow = length(alignment.files),
                                   ncol = length(feat.headers)))
  setnames(feature.data, feat.headers)
  feature.data[, Marker:=as.character(Marker)]

  #Gets the names and number of columns
  names(align.list) = gsub(".*/", "", alignment.files)
  n.cols = length(align.list)

  #For the concatenated alignment
  concat.headers = c("Sample", names(align.list))
  concat.data = data.table(matrix(as.character("n"), nrow = length(sample.names),
                                  ncol = length(concat.headers)))
  setnames(concat.data, concat.headers)
  concat.data[, Sample:=as.character(sample.names)]

  # Add columns with generated name with standardized length
  start = 1
  for (x in 2:ncol(concat.data)) {

    align = data.table(align.list[[x-1]])
    align.len = as.numeric(colnames(align)[2])
    colnames(align) = c("Sample", "Sequence")

    #This fills in missing samples
    current.samples = unlist(c(align[,1]))
    names(current.samples) = NULL
    miss.samples = sample.names[!sample.names %in% current.samples]

    #Fill in missing samples
    if (length(miss.samples) != 0){
      #Adds Ns in for each missing sample
      for (j in 1:length(miss.samples)){
        temp.seq = paste0(rep("n", align.len), collapse = "")
        save.seq = data.table(Sample = miss.samples[j],
                              Sequence = temp.seq)
        align = rbind(align, save.seq)
      }#end j loop
    } #end if

    align = align[order(align$Sample),]
    align = align[duplicated(align$Sample) != T,]

    concat.data = concat.data[order(concat.data$Sample),]

    set(concat.data, i = match(concat.data$Sample, align$Sample),
        j = as.integer(x), value = align$Sequence)

    align.name = colnames(concat.data)[x]

    if (start == 1){ end = start + align.len - 1 } else { end = start + align.len - 1 }

    #Saves location data
    set(feature.data, i = as.integer(x-1), j = match("Marker", feat.headers), value =  align.name)
    set(feature.data, i = as.integer(x-1), j = match("Length", feat.headers), value =  align.len)
    set(feature.data, i = as.integer(x-1), j = match("Start", feat.headers), value =  start)
    set(feature.data, i = as.integer(x-1), j = match("End", feat.headers), value =  end)
    start = end + 1
  }#end x loop

  #Now have to reformat somehow
  ntax = nrow(concat.data)
  nchar = sum(feature.data$Length)
  phy.header = paste(ntax, nchar)

  #Prep sample names to all have the same length padded with spaces
  name.length = max(nchar(sample.names)) + 5
  for (x in 1:length(sample.names)){
    temp.name = concat.data$Sample[x]
    space.pad = name.length - nchar(temp.name)
    space.add = paste(rep(" ", space.pad), collapse = "")
    new.name = paste0(temp.name, space.add, collapse = "")
    concat.data$Sample[x] = new.name
  }#end x

  #Concatenated each sample into 1 sequence and write it to file as a text file
  # cat(phy, sep = "\n")
  #write(phy, file = file)

  if (file.exists(paste0(output.dir, "/", file.name, ".phy")) == TRUE){
    unlink(paste0(output.dir, "/", file.name, ".phy"))}

  fileConn = file(paste0(output.dir, "/", file.name, ".phy"), open = "a")
  writeLines(phy.header, con = fileConn, sep = "\n")

  for (x in 1:nrow(concat.data)){

    sequence.data = paste0(concat.data[x,], collapse = "")
    writeLines(sequence.data, con = fileConn, sep = "\n")

  }#end x loop

  close(fileConn)

  print(paste0(file.name, ".phy concatenated alignment written in phylip format."))

  ##########################
  ##### Saves partition file
  ##########################
  #Saves in a simple table format
  if (length(partition.format[partition.format == "table"]) == 1){
    #Checks to overwrite
    if (file.exists(paste0(output.dir, "/", file.name, "_partitions_table.txt")) == TRUE){
      unlink(paste0(output.dir, "/", file.name, "_partitions_table.txt"))}

    write.table(feature.data, file = paste0(output.dir, "/", file.name, "_partitions_table.txt"),
                sep = "\t", row.names = FALSE, quote = FALSE)

    print(paste0(file.name, "_partitions_table.txt written in table format."))
  }#ends table if

  #Save in raxml format
  if (length(partition.format[partition.format == "raxml"]) == 1){
    #Checks to overwrite
    if (file.exists(paste0(output.dir, "/", file.name, "_partitions_raxml.txt")) == TRUE){
      unlink(paste0(output.dir, "/", file.name, "_partitions_raxml.txt"))}

    fileConn = file(paste0(output.dir, "/", file.name, "_partitions_raxml.txt"), open = "a")

    for (x in 1:nrow(feature.data)){

      writeLines(paste0("DNA, ", feature.data$Marker[x], " = ",
                        feature.data$Start[x], "-", feature.data$End[x]),
                 con = fileConn, sep = "\n")

    }#end x loop

    close(fileConn)

    print(paste0(file.name, "_partitions_raxml.txt written in raxml format."))

  }#end raxml if

}#end function

