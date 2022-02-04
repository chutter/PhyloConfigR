#' @title hyphyOmega
#'
#' @description Function for batch trimming a folder of alignments, with the various trimming functions available to select from
#'
#' @param genome.directory path to a folder of sequence alignments in phylip format.
#'
#' @param output.directory available input alignment formats: fasta or phylip
#'
#' @param threads contigs are added into existing alignment if algorithm is "add"
#'
#' @param threads path to a folder of sequence alignments in phylip format.
#'
#' @param memory give a save name if you wnat to save the summary to file.
#'
#' @param overwrite TRUE to supress mafft screen output
#'
#' @param resume TRUE to supress mafft screen output
#'
#' @return an alignment of provided sequences in DNAStringSet format. Also can save alignment as a file with save.name
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

hyphyOmega = function(tree.directory = NULL,
                      alignment.directory = NULL,
                      dataset.name = "omega",
                      threads = 1,
                      memory = 1,
                      overwrite = FALSE,
                      resume = TRUE,
                      quiet = TRUE,
                      hyphy.path = NULL,
                      mg94.path = NULL) {


  #Read in basic genome info
  # library(PhyloCap)
  # setwd("/Volumes/Rodents/Australian_Rodents/Data_Processing")
  # tree.directory= "/Volumes/Rodents/Australian_Rodents/Data_Processing/Trees/Ausfull/genes_trimmed_trees"
  # alignment.directory = "/Volumes/Rodents/Australian_Rodents/Data_Processing/Alignments/Ausfull/coding_trimmed/nt"
  # metadata.file = "/Volumes/Rodents/Australian_Rodents/Data_Processing/Mus-selected-sequences_metadata_final.csv"
  # dataset.name = "Omega"
  # threads = 4
  # memory = 4
  # resume = T
  # overwrite = F
  # hyphy.path = "/usr/local/bin"
  # mg94.path = "/Users/chutter/hyphy-analyses/FitMG94"

  #Directoires
  tree.directory = "/Users/chutter/Dropbox/Research/1_Main-Projects/0_Working-Projects/Rodent_Mitochondrial/Align-Trees/nuclear/trees_oxphos-coding"
  alignment.directory = "/Users/chutter/Dropbox/Research/1_Main-Projects/0_Working-Projects/Rodent_Mitochondrial/Align-Trees/nuclear/alignments_oxphos-coding"
  dataset.name = "Omega"
  threads = 4
  memory = 8
  resume = T
  overwrite = F
  quiet = T
  hyphy.path = "/usr/local/bin/"


  #Same adds to bbmap path
  if (is.null(hyphy.path) == FALSE){
    b.string = unlist(strsplit(hyphy.path, ""))
    if (b.string[length(b.string)] != "/") {
      hyphy.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { hyphy.path = "" }

  #Same adds to bbmap path
  if (is.null(mg94.path) == FALSE){
    b.string = unlist(strsplit(mg94.path, ""))
    if (b.string[length(b.string)] != "/") {
      mg94.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { mg94.path = "" }

  if (is.null(tree.directory) == T){ stop("A directory of trees is needed.") }
  if (is.null(alignment.directory) == T){ stop("A directory of alignments is needed.") }

  #So I don't accidentally delete everything while testing resume
  if (resume == TRUE & overwrite == TRUE){
    overwrite = FALSE
    stop("Error: resume = T and overwrite = T, cannot resume if you are going to delete everything!")
  }

  #Creates main analysis directory if its not already there
  if (dir.exists("hyphy") == FALSE) { dir.create("hyphy") }

  output.directory = paste0("hyphy/", dataset.name)
  if (dir.exists(output.directory) == TRUE){
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.directory))
      dir.create(output.directory)
    }
  } else { dir.create(output.directory) }

  align.files = list.files(alignment.directory)
  tree.files = list.files(tree.directory)

  #Resumes file download
  if (resume == TRUE){
    done.files = list.files(output.directory)
    done.files = done.files[duplicated(done.files) == TRUE]
    tree.files = tree.files[!gsub("\\..*", "", tree.files) %in% gsub("\\..*", "", done.files)]
  }

  #Stats table prepare
  header.data = c("sample", "dn", "ds", "omega")
  collect.data = data.table::data.table(matrix(as.numeric(0), nrow = length(tree.files), ncol = length(header.data)))
  data.table::setnames(collect.data, header.data)
  collect.data[, sample:=as.character(sample)]

  #Run analyses through all trees
  for (i in 1:length(tree.files)){

    #Find correct gene tree
    locus.name = gsub("\\..*", "", tree.files[i])
    # temp.meta = meta.data[meta.data$gene_id %in% locus.name,]
    # gene.name = unique(temp.meta$gene_id)
    # prot.name = unique(temp.meta$protein_id)
    # prot.name = prot.name[is.na(prot.name) != T]
    #
    # if (length(prot.name) == 0){
    #   print(paste0(locus.name, " protein not found in the alignments. Skipped."))
    #   next }
    #

    align.name = align.files[grep(locus.name, align.files)]

    if (length(align.name) == 0){
      print(paste0(locus.name, " not found in the alignments. Skipped."))
      next }

    #Read in alignment
    temp.align = Biostrings::readDNAStringSet(paste0(alignment.directory, "/", align.name))
    #Prune gene tree to match alignment
    temp.tree = ape::read.tree(paste0(tree.directory, "/", tree.files[i]))

    #First remove sample from alignment not in tree
    temp.align = temp.align[names(temp.align) %in% temp.tree$tip.label]
    #Remove tips from tree not in alignment
    temp.tree = ape::drop.tip(temp.tree, temp.tree$tip.label[!temp.tree$tip.label %in% names(temp.align)])

    if (length(temp.align) <= 3){
      print(paste0(locus.name, " alignment has too few taxa. Skipped."))
      next }

    if (length(temp.tree$tip.label) <= 3){
      print(paste0(locus.name, " alignment has too few taxa. Skipped."))
      next }

    if (is.null(temp.tree) == TRUE){
      print(paste0(locus.name, " alignment has too few taxa. Skipped."))
      next }

    dir.create(paste0(output.directory, "/", locus.name))
    #Save both as temp files to be used as input in hyphy
    write.temp = strsplit(as.character(temp.align), "")
    aligned.set = as.matrix(ape::as.DNAbin(write.temp) )
    writePhylip(alignment = aligned.set,
                file=paste0(output.directory, "/", locus.name, "/", locus.name, ".phy"),
                interleave = F,
                strict = F)

    ape::write.tree(temp.tree, file = paste0(output.directory, "/", locus.name, "/", locus.name, ".tre"))

    #Gets Dn/Ds from experiment hyphy function will probably have to change later
    system(paste0(hyphy.path, "hyphy FitMG94.bf",
                  " --alignment ", output.directory, "/", locus.name, "/", locus.name, ".phy",
                  " --tree ", output.directory, "/", locus.name, "/", locus.name, ".tre",
                  " --output ", output.directory, "/", locus.name, "/", dataset.name, "-fitmg94-results.json",
                  " --type local --lrt Yes CPU=", threads))

    json.data = jsonlite::fromJSON(paste0(output.directory, "/", locus.name, "/", dataset.name, "-fitmg94-results.json"))
    branch.att = json.data$`branch attributes`
    branch.att = unlist(branch.att[[1]])

    #Gather stats
    dn.stats = branch.att[grep("dN", names(branch.att) )]
    ds.stats = branch.att[grep("dS", names(branch.att) )]

    #Make table from stats
    ds.data = data.frame(Sample = gsub(".dS", "", names(ds.stats)),
               dS = as.numeric(ds.stats))

    dn.data = data.frame(Sample = gsub(".dN", "", names(dn.stats)),
               dN = as.numeric(dn.stats))

    w.data = merge(ds.data, dn.data, by = "Sample")
    w.data$omega = w.data$dN/w.data$dS

    write.csv(w.data, file = paste0(output.directory, "/", locus.name, "/dn-ds_stats.csv"), row.names = F, quote = F)

    print(paste0(locus.name, " finished omega calculation!"))

    #Gets Dn/Ds from experiment hyphy function will probably have to change later
    # system(paste0(hyphy.path, "hyphy ", mg94.path, "FitMG94.bf",
    #               " --alignment ", output.directory, "/", gene.name, "/", gene.name, ".phy",
    #               " --tree ", output.directory, "/", gene.name, "/", gene.name, ".tre",
    #               " --output ", output.directory, "/", gene.name, "/", dataset.name, "-fitmg94-results.json",
    #               " --type local --lrt Yes CPU=", threads))
    #


    #Runs hyphy
# #    system(paste0(hyphy.path, "hyphy busted --alignment ", output.directory, "/", gene.name, "/", gene.name, ".phy",
#                   " --tree ", output.directory, "/", gene.name, "/", gene.name, ".tre",
#                   " --output ", output.directory, "/", gene.name, "/", dataset.name, "-results.json"))
#
#     json.data = readLines(paste0(output.directory, "/", gene.name, "/", dataset.name, "-results.json"), warn = F)
#     x = grep("dN & dS for each branch", mlc)
#     y = grep("Time used", mlc)
#
#
#
#   hyphy <method_name> --alignment <path_to_alignment_file> <additional_method_specific_arguments>
#
#
#     echo `(echo "10"; echo "4"; echo "1"; echo "/path/to/downloaded/example/files/HIV.nex"; echo "Y"; echo "1"; echo "d") | HYPHYMP`


  }#end i loop


}#end function


