#' @title hyphy.Parser
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

hyphy.Parser = function(results.directory = NULL,
                        slac.spreadsheet = NULL,
                        output.name = NULL,
                        tips.only = TRUE,
                        threads = 1,
                        memory = 1,
                        overwrite = FALSE,
                        quiet = TRUE,
                        hyphy.path = NULL) {


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
  setwd("/Users/chutter/Dropbox/Research/1_Main-Projects/0_Working-Projects/Rodent_Mitochondrial")
  results.directory = "/Volumes/Rodents/Australian_Rodents/Data_Processing/hyphy/slac_all"
  output.name = "slac_all"
  slac.spreadsheet = paste0(output.name, "_SLAC_by-branch.csv")
  species.tree = "/Users/chutter/Dropbox/Research/1_Main-Projects/0_Working-Projects/Rodent_Mitochondrial/concat_tree.tre"
  threads = 4
  memory = 8
  overwrite = F
  quiet = T
  tips.only = T
  hyphy.path = "/usr/local/bin/"
  outgroups = c("Rattus_morotaiensis_ASAM29")

  #Same adds to bbmap path
  if (is.null(hyphy.path) == FALSE){
    b.string = unlist(strsplit(hyphy.path, ""))
    if (b.string[length(b.string)] != "/") {
      hyphy.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { hyphy.path = "" }

  if (is.null(results.directory) == T){ stop("A directory of results is needed.") }

  slac.results = data.table::fread(slac.spreadsheet)

  slac.results = slac.results[is.na(slac.results$omega) != T,]
  slac.results = slac.results[is.infinite(slac.results$omega) != T,]
  gene.results = aggregate(x = slac.results, by = list(slac.results$file), FUN = mean)

  mean.slac = aggregate(x = slac.results, by = list(slac.results$sample), FUN = mean)
  sample.omega = mean.slac$omega
  names(sample.omega) = mean.slac$Group.1

  slac.tree = ape::read.tree(species.tree)
  slac.tree = ape::root(slac.tree, outgroup = outgroups, resolve.root = TRUE)
  slac.tree$tip.label = gsub("-", "_", slac.tree$tip.label)
  sample.omega = sample.omega[names(sample.omega) %in% slac.tree$tip.label]

  phytools::contMap(slac.tree, log10(sample.omega), res=100, fsize=NULL, ftype=NULL, lwd=4, legend=NULL,
          lims=NULL, sig=3, type="phylogram", direction="rightwards")


  #OLDDDD compare dn/ds branch among groups
  hill.dir = list.files("/Volumes/Rodents/Australian_Rodents/Data_Processing/hyphy/busted_hill")
  gen.dir = list.files("/Volumes/Rodents/Australian_Rodents/Data_Processing/hyphy/busted_genbank")

  hill.results = gene.results[gene.results$Group.1 %in% hill.dir,]
  bg.results = gene.results[!gene.results$Group.1 %in% hill.dir,]
  gen.results = gene.results[gene.results$Group.1 %in% gen.dir,]

  mean(bg.results$omega)
  mean(hill.results$omega)
  mean(gen.results$omega)

  gg.data.a = data.frame(dataset = "Mito-interactors", omega = hill.results$omega)
  #gg.data.b = data.frame(dataset = "GenBank", omega = gen.results$omega)
  gg.data.c = data.frame(dataset = "Background", omega = bg.results$omega)
  plot.data = rbind(gg.data.a, gg.data.c)

  library(ggplot2)
  ggplot(plot.data, aes(x=omega, fill=dataset)) +
    geom_density(alpha = 0.4)

  ggplot(plot.data, aes(x=dataset, y=omega)) +
    geom_boxplot()

  #Trees of those
  phytools::contMap(slac.tree, log10(sample.omega), res=100, fsize=NULL, ftype=NULL, lwd=4, legend=NULL,
                    lims=NULL, sig=3, type="phylogram", direction="rightwards")


  ####
  samples.sheet = read.csv("/Users/chutter/Dropbox/Research/1_Main-Projects/0_Working-Projects/Rodent_Mitochondrial/Sample_Data.csv")
  monsoon.samples = samples.sheet[samples.sheet$Habitat == "Monsoon",]
  arid.samples = samples.sheet[samples.sheet$Habitat == "Arid",]
  mesic.samples = samples.sheet[samples.sheet$Habitat == "Mesic",]

  mean.slac = aggregate(x = slac.results, by = list(slac.results$sample), FUN = mean)

  monsoon.data = mean.slac[mean.slac$Group.1 %in% monsoon.samples$Sample_Name,]
  gg.data.a = data.frame(dataset = "Monsoon", omega = monsoon.data$omega)
  arid.data = mean.slac[mean.slac$Group.1 %in% arid.samples$Sample_Name,]
  gg.data.b = data.frame(dataset = "Arid", omega = arid.data$omega)
  mesic.data = mean.slac[mean.slac$Group.1 %in% mesic.samples$Sample_Name,]
  gg.data.c = data.frame(dataset = "Mesic", omega = mesic.data$omega)
  plot.data = rbind(gg.data.a, gg.data.b, gg.data.c)


  ggplot(plot.data, aes(x=dataset, y=omega)) +
    geom_boxplot() + geom_dotplot(binaxis='y', stackdir='center', dotsize=1)

}#end function


