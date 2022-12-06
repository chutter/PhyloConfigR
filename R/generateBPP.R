#' @title generateBPP
#'
#' @description Function for generating BPP input alignment and IMAP files.
#'
#' @param alignment.directory path to a folder of sequence alignments in phylip format.
#'
#' @param tree.file a tree file that can be read by the ape R package function read.tree
#'
#' @param outgroups provide outgroups for rooting the tree, can use multiple like so: c("taxa_1", "taxa_2", "taxa_3)
#'
#' @param output.name name for output files
#'
#' @param population.file "provided" if you already have a population file and do not need one created. "letter" if you want to create a population file for each sample is a population and given a 3 letter code. "sample" if you want to have the population named after each sample.
#'
#' @param overwrite if TRUE overwrites previous runs of the function
#'
#' @return a single BPP formatted concatenated alignment and optionally an IMAP file
#'
#' @examples
#'
#'Situation 1. Population file provided so an IMAP file is not created.
#'
#'generateBPP(alignment.directory = directory/to/alignments,
#'            tree.file = path/to/tree/file,
#'            outgroups = c("taxa_1"),
#'            output.name = "bpp-alignment",
#'            population.file = "provided",
#'            overwrite = FALSE)
#'
#'Situation 2. No population file, name each population with the sample name. The sample and population name will be the same in the IMAP.
#'
#'generateBPP(alignment.directory = directory/to/alignments,
#'            tree.file = path/to/tree/file,
#'            outgroups = c("taxa_1", "taxa_2"),
#'            output.name = "bpp-alignment",
#'            population.file = "sample",
#'            overwrite = FALSE)
#'
#'Situation 3. No population file, each sample is a single population. Each sample is named with a 3 letter long random letter code to match the population IMAP file.
#'
#'generateBPP(alignment.directory = directory/to/alignments,
#'            tree.file = path/to/tree/file,
#'            outgroups = c("taxa_1", "taxa_2", "taxa_3"),
#'            output.name = "bpp-alignment",
#'            population.file = "letter",
#'            overwrite = FALSE)
#'
#' @export

generateBPP = function(alignment.directory = NULL,
                       tree.file = NULL,
                       outgroups = NULL,
                       output.name = "bpp-alignment",
                       population.file = c("letter", "sample", "provided"),
                       overwrite = FALSE) {

  # #Debug setup
  # library(PhyloCap)
  # setwd("/Volumes/LaCie/Hylidae/BPP/Hylidae")
  # alignment.directory = "/Volumes/LaCie/Hylidae/Alignments/locus-combined"
  # output.name = "gene"
  # population.file = "sample"
  # overwrite = TRUE
  #
  # tree.file = "/Volumes/LaCie/Hylidae/Alignments/Filtered/Best/all-markers_trimmed_filters_S70-L3000-PIN10-PIX50.phy.treefile"
  # outgroups = c("Phyllomedusa_tomopterna_WED_55506", "Nyctimystes_infrafrenatus_SLT_771")

  #Gets list of alignments
  align.files = list.files(alignment.directory, full.names = T)
  align.list = lapply(align.files, function (x) data.table::fread(x, header = T))
  sample.names = unique(unlist(lapply(align.list, function (x) x[,1])))
  sample.names = sample.names[order(sample.names)]

  if (population.file == "letter"){
    sample.table = data.frame(original = sample.names,
                              tag = rep(1:length(sample.names)),
                              bpp = paste0(sample.names, "^", rep(1:length(sample.names))),
                              pop = stringi::stri_rand_strings(length(sample.names), 3, "[A-Z]"))
  }#end letter

  if (population.file == "sample"){
    sample.table = data.frame(original = sample.names,
                              tag = sample.names,
                              bpp = paste0( "^", sample.names),
                              pop = sample.names)
  }#end letter

  #Gets the names and number of columns
  names(align.list) = gsub(".*/", "", align.files)
  n.loci = length(align.list)

  if (overwrite == TRUE){
    if (file.exists(paste0(output.name, "_alignment.txt")) == TRUE){ unlink(paste0(output.name, "_alignment.txt")) }
    if (file.exists(paste0(output.name, "_species-tree.txt")) == TRUE){ unlink(paste0(output.name, "_species-tree.txt")) }
    if (file.exists(paste0(output.name, "_count-table.txt")) == TRUE){ unlink(paste0(output.name, "_count-table.txt")) }
  }#end if

  if (overwrite == FALSE){
    if (file.exists(paste0(output.name, "_alignment.txt")) == TRUE){ stop("overwrite == FALSE and files present.") }
    if (file.exists(paste0(output.name, "_species-tree.txt")) == TRUE){ stop("overwrite == FALSE and files present.") }
    if (file.exists(paste0(output.name, "_count-table.txt")) == TRUE){ stop("overwrite == FALSE and files present.") }
  }#end if

  fileConn = file(paste0(paste0(output.name, "_alignment.txt")), open = "a")

  #Loops through each locus and does operations on them
  for (i in 1:length(align.list)){

    align.data = align.list[[i]]
    ###################
    #Now have to reformat
    ntax = as.character(nrow(align.data))
    #ntax = length(sample.names)
    nchar = colnames(align.data)[2]
    phy.header = paste0(ntax, " ", nchar)

    #Prep sample names to all have the same length padded with spaces
    name.length = max(nchar(sample.names)) + 4
    for (x in 1:nrow(align.data)){
      temp.name = align.data[x,1]
      #adds space padding
      int.name = sample.table[sample.table$original == as.character(temp.name),]$bpp
      save.name = sample.table[sample.table$original == as.character(temp.name),]$pop
      space.pad = name.length - nchar(int.name)
      space.add = paste(rep(" ", space.pad), collapse = "")
      new.name = paste0(int.name, space.add, collapse = "")
      align.data[x,1] = new.name

      #Uppercase
      temp.align = align.data[x,2]
      temp.align = gsub("-", "N", temp.align)
      temp.align = gsub("\\?", "N", temp.align)
      align.data[x,2] = toupper(temp.align)
    }#end x

    writeLines(phy.header, con = fileConn, sep = "\n")
    #writeLines("\n", con = fileConn, sep = "")

    #Saves each line
    for (x in 1:nrow(align.data)){
      sequence.data = paste0(align.data[x,], collapse = "")
      writeLines(sequence.data, con = fileConn, sep = "\n")
    }#end x loop

    # #Writes missin gtaxa
    # missing.taxa = sample.names[!sample.names %in% save.names]
    # #Saves each line
    # for (x in 1:length(missing.taxa)){
    #   space.pad = name.length - nchar(temp.name)
    #   space.add = paste(rep(" ", space.pad), collapse = "")
    #   new.name = paste0(missing.taxa[x], space.add, collapse = "")
    #   sequence.data = paste0(new.name, paste0(rep("N", as.numeric(nchar)), collapse = ""))
    #   writeLines(sequence.data, con = fileConn, sep = "\n")
    # }#end x loop

    writeLines("\n", con = fileConn, sep = "")
  }#end loop

  #Close and end
  close(fileConn)

  #Creates Imap file
  imap.save = data.frame(sample.table$tag, sample.table$pop)
  colnames(imap.save) = NULL
  write.table(imap.save, file = paste0(output.name, "_Imap.txt"),
              quote = FALSE, sep = " ", row.names = F)

  #Configures tree
  input.tree = ape::read.tree(file = tree.file)
  input.tree = ape::root(input.tree, outgroups, resolve.root = T)
  input.tree$node.label = NULL
  input.tree$edge.length = NULL

  #Change names
  for (x in 1:nrow(sample.table)){
    input.tree$tip.label[input.tree$tip.label %in% sample.table$original[x]] = sample.table$pop[x]
  }

  ape::write.tree(input.tree, file = paste0(output.name, "_species-tree.txt"))

  #Make count table
  count.table = data.frame(table(sample.table$pop))
  colnames(count.table) = NULL
  count.table = transpose(count.table)
  write.table(count.table, file = paste0(output.name, "_count-table.txt"),
              quote = FALSE, sep = " ", row.names = F, col.names = F)

  print(paste0(output.name, " concatenated alignment written BPP format. Generic IMAP written."))

} #end function

#########################
###### END SCRIPT
#########################


