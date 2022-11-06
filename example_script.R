
#Load in packages
devtools::install_github("chutter/AstralPlane")
library(AstralPlane)
library(data.table)

#To do:
### Make name simpler
### Use hyphen
#### Use astral results to calcualte anomaly zone across everything
#### Collect node data from anomaly zone

######################################################################################
##### Testing the anomaly zone
work.dir = "/Volumes/Armored/Hylidae/Dataset-single"
astral.path = "/usr/local/bin/Astral-5-14/astral.5.14.2.jar"
dir.create(work.dir)
setwd(work.dir)

tree.file = "/Users/chutter/Dropbox/Research/1_Main-Projects/0_Working-Projects/Hylidae/Trees/Tree_Grid/Unfiltered/Astral/uce.tre"
outgroup.taxa = c("Phyllomedusa_tomopterna_WED_55506", "Nyctimystes_infrafrenatus_SLT_771")

uce.tree = ape::read.tree(tree.file)
anom.data = anomalyZone(tree = uce.tree,
                        outgroups = outgroup.taxa)

plot.anomalyZone(tree = uce.tree,
                 data = anom.data,
                 outgroups = outgroup.taxa,
                 save.file = NULL,
                 tip.label.size = 0.5,
                 node.label.size = 1,
                 edge.width = 3)

######################################################################################
##### Filtering demonstration

#Set up your directories
align.dir = "/Volumes/Armored/Hylidae/Alignments/all-markers_trimmed"
tree.dir = "/Users/chutter/Dropbox/Research/1_Main-Projects/0_Working-Projects/Hylidae/Trees/Gene_Trees/all-markers_trimmed"
astral.path = "/usr/local/bin/Astral-5-14/astral.5.14.2.jar"
dir.create(work.dir)
setwd(work.dir)

#Set up your filters
filter.length = c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
                  1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000,
                  2100, 2200, 2300, 2400, 2500) #number of base pairs
filter.sample = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1) #proportion samples
filter.prop.pis = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1) #proportion pis
filter.count.pis = c(10, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500) #count of pis

#Alignment summary function [20 minutes]
align.summary = summarizeAlignments(alignment.path = align.dir,
                                    file.export = "Hylidae_Alignment_stats",
                                    alignment.format = "phylip",
                                    dataset.name = "Hylidae",
                                    overwrite = TRUE)

#Apply filters and create summary table of filters [< 1 min]
filt.summary = filterSummary(alignment.data = align.summary,
                             alignment.folder = align.dir,
                             dataset.name  = "Hylidae",
                             file.out = "filter_summary",
                             length.filters = filter.length,
                             sample.filters = filter.sample,
                             prop.pis.filters = filter.prop.pis,
                             count.pis.filters = filter.count.pis,
                             overwrite = T)

#Make filtered alignments datasets [10 minutes]
filterAlignments(filter.summary = filt.summary,
                 alignment.data = align.summary,
                 alignment.folder = align.dir,
                 format = "concatenated",
                 min.alignments = 5,
                 overwrite = TRUE)

#Make filtered gene trees datasets [5 minutes]
filterGeneTrees(filter.summary = filt.summary,
                alignment.data = align.summary,
                genetree.folder = tree.dir,
                format = "concatenated",
                overwrite = TRUE,
                taxa.remove = NULL,
                min.trees = 5,
                min.n.samples = 4,
                min.sample.prop = NULL,
                make.polytomy = TRUE,
                polytomy.limit = 10,
                remove.node.labels = FALSE)

#Runs astral across all filtered gene tree sets
AstralPlane::astralRunner(concat.genetree.folder = "filtered-genetrees-concatenated",
                          output.dir = "filtered-Astral",
                          overwrite = FALSE,
                          astral.path = astral.path,
                          astral.t = 2,
                          quiet = FALSE,
                          multi.thread = TRUE,
                          memory = "8g")

#Concordance factor analysis across all filtered datasets
AstralPlane::concordanceRunner(alignment.dir = "filtered-alignments-concatenated",
                               species.tree.dir = "filtered-Astral",
                               genetree.dir = "filtered-genetrees-concatenated",
                               output.dir = "concordance-factors",
                               overwrite = FALSE,
                               quiet = TRUE,
                               threads  = 6)

######################################################################################
##### Testing the effectiveness of the anomaly zone

#Set up your directories
work.dir = "/Volumes/Armored/Hylidae/Dataset-single"
align.dir = "/Volumes/Armored/Hylidae/Alignments/all-markers_trimmed"
tree.dir = "/Users/chutter/Dropbox/Research/1_Main-Projects/0_Working-Projects/Hylidae/Trees/Gene_Trees/all-markers_trimmed"
setwd(work.dir)

align.summary = read.csv("Hylidae_Alignment_stats.csv")
filt.summary = read.csv("filter_summary.csv")

#Set up a list of your clades of interest
outgroup.taxa = c("Phyllomedusa_tomopterna_WED_55506", "Nyctimystes_infrafrenatus_SLT_771")
taxa.set = list()
taxa.set[[1]] = c("Colomascirtus_staufferorum_LAC_2153", "Hyloscirtus_phyllognathus_WED_58378", "Hypsiboas_boans_WED_57877")
taxa.set[[2]] = c("Sphaenorhynchus_lacteus_WED_54090", "Scinax_garbei_LAC_1787", "Scinax_ruber_WED_56140" )
taxa.set[[3]] = c("Trachycephalus_jordani_WED_53658", "Phyrnohyas_venulosa_WED_55450", "Osteocephalus_taurinus_WED_55452")
taxa.set[[4]] = c("Scarthyla_goinorum_WED_58246" , "Lysapsus_laevis_CAS_257655", "Pseudis_paradoxus_CAS_245053")
taxa.set[[5]] = c("Dendropsophus_leucophyllatus_WED_59288", "Dendropsophus_parviceps_MZUTI_1357", "Dendropsophus_koechlini_WED_57879")
taxa.set[[6]] = c("Plectrohyla_quecchi_MVZ_251534" , "Ptychohyla_salvadorensis_EBG_518","Smilisca_phaeota_LAC_2299",
                  "Hyla_sarda_WED_54544", "Hyla_walkeri_MVZ_263408", "Dryophtes_cinerea_WED_56355")
#taxa.set[[7]] = c("Acris_blanchardi_DSM_2012", "Pseudacris_triseriata_BLO_006","Hyliola_cadaverina_WED_54461")
names(taxa.set) = c("node1", "node2", "node3", "node4", "node5", "node6")


#### Function start
anomaly.data = filterAnomalies(astral.directory = "filtered-Astral",
                                outgroups = outgroup.taxa,
                                filter.data = filt.summary)

#Obtains the concordance factors data for all filtered replicates
concord.data = filterConcordance(input.dir = "concordance-factors",
                                 clade.list = taxa.set,
                                 outgroups  = outgroup.taxa)

#### Pull out best tree
bestFilterTrees(anomaly.zone.data = anomaly.data,
                concordance.factors.data = concord.data,
                output.dir = "best-trees",
                min.trees = 20,
                fewest.anomaly.zones = TRUE,
                highest.gene.cf = TRUE,
                highest.post.prob = TRUE,
                all.datasets = TRUE,
                top.best = 1)

# General plotting function to get a sense of the results without specific node
plot.filterZone(anomaly.zone.data = anomaly.data,
                concordance.factors.data = concord.data,
                save.plots = TRUE,
                output.dir = "Filter-Plots",
                dataset.name = "all",
                plot.gcf = TRUE,
                plot.scf = TRUE,
                az.colors = c("#7BC143", "#DE3293"),
                m.shape = c(22, 21),
                min.trees = 10)

#Plot alignment length, node 2
plot.filterNode(anomaly.zone.data = anomaly.data,
                concordance.factors.data = concord.data,
                output.dir = "Filter-Plots",
                focal.node = "node2",
                filter.name = "alignment_length",
                dataset.name = "all",
                plot.gcf = TRUE,
                plot.scf = TRUE,
                az.colors = c("#7BC143", "#DE3293"),
                m.shape = c(22, 21),
                min.trees = 10 )

#Plot alignment length, node 6
plot.filterNode(anomaly.zone.data = anomaly.data,
                concordance.factors.data = concord.data,
                output.dir = "Filter-Plots",
                focal.node = "node6",
                filter.name = "alignment_length",
                dataset.name = "all",
                plot.gcf = TRUE,
                plot.scf = TRUE,
                az.colors = c("#7BC143", "#DE3293"),
                m.shape = c(22, 21),
                min.trees = 10 )

######################################################################################
######################################################################################
######################################################################################
######################################################################################
######################################################################################
##### Filtering demonstration with multi-dataset loop

#Load in packages
devtools::install_github("chutter/AstralPlane")
library(AstralPlane)
library(data.table)

#Set up your directories
align.dir = "/Volumes/Armored/Hylidae/Alignments"
tree.dir = "/Users/chutter/Dropbox/Research/1_Main-Projects/0_Working-Projects/Hylidae/Trees/Gene_Trees"
work.dir = "/Volumes/Armored/Hylidae/Dataset-filtering"
astral.path = "/usr/local/bin/Astral-5-14/astral.5.14.2.jar"
dir.create(work.dir)
setwd(work.dir)

#If you have many different datasets to set up and run, you can give it the folder
#of various datasets.
batchAstral(genetree.datasets = tree.dir,
            astral.t = 2,
            output.dir = "test-dataset",
            min.n.samples = 4,
            min.sample.prop = 0.1,
            taxa.remove = NULL,
            overwrite = TRUE,
            quiet = F,
            astral.path = astral.path,
            make.polytomy = TRUE,
            polytomy.limit = 10,
            multi.thread = TRUE,
            memory = "8g")

#Obtains dataset names
datasets = list.dirs(tree.dir, full.names = F, recursive = F)

for (i in 1:length(datasets)){

  #Read in the astral data and tree and organize it into different slots
  dataset.name = paste0("test-dataset/", datasets[i], "_astral.tre")
  astral.data = createAstralPlane(astral.tree = dataset.name,
                                  outgroups = outgroups,
                                  tip.length = 1)

  #Plots the astral data
  astralProjection(astral.plane = astral.data,
                   local.posterior = TRUE,
                   pie.plot = TRUE,
                   pie.data = "qscore",
                   save.file = paste0("test-dataset/", datasets[i], ".pdf"),
                   pie.colors = c("purple", "blue", "green"),
                   node.color.text = c("white"),
                   node.color.bg = c("black"),
                   tip.label.size = 0.75,
                   pie.chart.size = 1)

  species.tree = ape::read.tree(paste0(astral.dir, "/", astral.trees[i]))
  anom.data = anomalyZone(tree = species.tree,
                          outgroups = outgroup.taxa)

  plot.anomalyZone(tree = species.tree,
                   data = anom.data,
                   outgroups = outgroup.taxa,
                   save.file = paste0(out.dir, "/", astral.trees[i],"_anomaly_zone.pdf"),
                   tip.label.size = 0.4,
                   node.label.size = 0.3,
                   edge.width = 2)

}#end i loop


##### Filtering demonstration with multi-dataset loop
#######

#Set up your filters
filter.length = c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
                  1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000,
                  2100, 2200, 2300, 2400, 2500) #number of base pairs
filter.sample = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1) #proportion samples
filter.prop.pis = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1) #proportion pis
filter.count.pis = c(10, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500) #count of pis

#Gather datasets
datasets = c("all-markers_trimmed", "exon-only_trimmed",
             "intron-only_trimmed", "legacy-markers_trimmed",
             "locus-combined", "uce_trimmed")

master.align.summary = data.table()
master.filt.summary = data.table()
for (i in 1:length(datasets)){

  #makes the alignment directory name from the dataset
  dataset.align = paste0(align.dir, "/", datasets[i])
  dataset.trees = paste0(tree.dir, "/", datasets[i])

  #Alignment summary function [20 minutes]
  #Do not want to save a file. Do this at the end.
  #Use datasetalign here
  align.summary = summarizeAlignments(alignment.path = dataset.align,
                                      dataset.name = datasets[i],
                                      file.export = datasets[i],
                                      alignment.format = "phylip",
                                      overwrite = FALSE)

  #Apply filters and create summary table of filters [< 1 min]
  #Do not want to save a file. Do this at the end.
  filt.summary = filterSummary(alignment.data = align.summary,
                               alignment.folder = align.dir,
                               dataset.name  = datasets[i],
                               file.out = NULL,
                               length.filters = filter.length,
                               sample.filters = filter.sample,
                               prop.pis.filters = filter.prop.pis,
                               count.pis.filters = filter.count.pis)

  #Make filtered alignments datasets [10 minutes]
  #Use datasetalign here; overwrite = FALSE
  filterAlignments(filter.summary = filt.summary,
                   alignment.data = align.summary,
                   alignment.folder = dataset.align,
                   format = "concatenated",
                   min.alignments = 5,
                   min.n.samples = 4,
                   overwrite = FALSE)

  #Make filtered gene trees datasets [5 minutes]
  #Use datasetalign here; overwrite = FALSE
  filterGeneTrees(filter.summary = filt.summary,
                  alignment.data = align.summary,
                  genetree.folder = dataset.trees,
                  format = "concatenated",
                  overwrite = FALSE,
                  taxa.remove = NULL,
                  min.trees = 5,
                  min.n.samples = 5,
                  min.sample.prop = NULL,
                  make.polytomy = TRUE,
                  polytomy.limit = 10,
                  remove.node.labels = FALSE)

  master.align.summary = rbind(master.align.summary, align.summary)
  master.filt.summary = rbind(master.filt.summary, filt.summary)
  print(paste0(datasets[i], " complete!"))
}#end i loop

##### Save the files here since it too awhile to generate
write.csv(master.align.summary, file = "master_alignment_summary.csv", row.names = F)
write.csv(master.filt.summary, file = "master_filter_summary.csv", row.names = F)

#Runs astral across all filtered gene tree sets
AstralPlane::astralRunner(concat.genetree.folder = "filtered-genetrees-concatenated",
                          output.dir = "filtered-astral",
                          overwrite = TRUE,
                          astral.path = astral.path,
                          astral.t = 2,
                          quiet = FALSE,
                          multi.thread = TRUE,
                          memory = "8g")


#Concordance factor analysis across all filtered datasets
AstralPlane::concordanceRunner(alignment.dir = paste0(work.dir, "/filtered-alignments-concatenated"),
                               species.tree.dir = paste0(work.dir, "/filtered-astral"),
                               genetree.dir = paste0(work.dir, "/filtered-genetrees-concatenated"),
                               output.dir = "concordance-factors",
                               overwrite = TRUE,
                               quiet = TRUE,
                               threads  = 6)


######################################################################################
######################################################################################
######################################################################################
######################################################################################
######################################################################################
##### Analyze data all together

#Load in packages
devtools::install_github("chutter/AstralPlane")
library(AstralPlane)
library(data.table)

#Set up your directories
align.dir = "/Volumes/Armored/Hylidae/Alignments"
work.dir = "/Volumes/Armored/Hylidae/Dataset-filtering"
astral.dir = "filtered-astral"
setwd(work.dir)

#outgroups
outgroup.taxa = c("Phyllomedusa_tomopterna_WED_55506", "Nyctimystes_infrafrenatus_SLT_771")

#Set up a list of your clades of interest
taxa.set = list()
taxa.set[[1]] = c("Colomascirtus_staufferorum_LAC_2153", "Hyloscirtus_phyllognathus_WED_58378", "Hypsiboas_boans_WED_57877")
taxa.set[[2]] = c("Sphaenorhynchus_lacteus_WED_54090", "Scinax_garbei_LAC_1787", "Scinax_ruber_WED_56140" )
taxa.set[[3]] = c("Trachycephalus_jordani_WED_53658", "Phyrnohyas_venulosa_WED_55450", "Osteocephalus_taurinus_WED_55452")
taxa.set[[4]] = c("Scarthyla_goinorum_WED_58246" , "Lysapsus_laevis_CAS_257655", "Pseudis_paradoxus_CAS_245053")
taxa.set[[5]] = c("Dendropsophus_leucophyllatus_WED_59288", "Dendropsophus_parviceps_MZUTI_1357", "Dendropsophus_koechlini_WED_57879")
taxa.set[[6]] = c("Plectrohyla_quecchi_MVZ_251534" , "Ptychohyla_salvadorensis_EBG_518","Smilisca_phaeota_LAC_2299",
                  "Hyla_sarda_WED_54544", "Hyla_walkeri_MVZ_263408", "Dryophtes_cinerea_WED_56355")
taxa.set[[7]] = c("Pseudacris_triseriata_BLO_006","Hyliola_cadaverina_WED_54461")
names(taxa.set) = c("node1", "node2", "node3", "node4", "node5", "node6", "node7")

align.summary = read.csv("master_alignment_summary.csv")
filt.summary = read.csv("master_filter_summary.csv")

#### Collect anamaly zone data from all filtered replicates
anomaly.data = filterAnomalies(astral.directory = astral.dir,
                               outgroups = outgroup.taxa,
                               filter.data = filt.summary)

#Obtains the concordance factors data for all filtered replicates
concord.data = filterConcordance(input.dir = "concordance-factors",
                                 clade.list = taxa.set,
                                 outgroups  = outgroup.taxa)

# General plotting function to get a sense of the results without specific node
plot.filterZone(anomaly.zone.data = anomaly.data,
                concordance.factors.data = concord.data,
                save.plots = TRUE,
                output.dir = "Filter-Plots",
                dataset.name = "all",
                plot.gcf = TRUE,
                plot.scf = TRUE,
                az.colors = c("#7BC143", "#DE3293"),
                m.shape = c(22, 21),
                min.trees = 10)

#### Pull out best tree
bestFilterTrees(anomaly.zone.data = anomaly.data,
                concordance.factors.data = concord.data,
                output.dir = "best-trees",
                min.trees = 20,
                fewest.anomaly.zones = TRUE,
                highest.gene.cf = TRUE,
                highest.post.prob = TRUE,
                all.datasets = TRUE,
                top.best = 1)

#Plot alignment length, node 2
plot.filterNode(anomaly.zone.data = anomaly.data,
                concordance.factors.data = concord.data,
                output.dir = "Filter-Plots",
                focal.node = "node2",
                filter.name = "alignment_length",
                dataset.name = "all",
                plot.gcf = TRUE,
                plot.scf = TRUE,
                az.colors = c("#7BC143", "#DE3293"),
                m.shape = c(22, 21),
                min.trees = 10 )

#Plot alignment length, node 6
plot.filterNode(anomaly.zone.data = anomaly.data,
                concordance.factors.data = concord.data,
                output.dir = "Filter-Plots",
                focal.node = "node6",
                filter.name = "alignment_length",
                dataset.name = "all",
                plot.gcf = TRUE,
                plot.scf = TRUE,
                az.colors = c("#7BC143", "#DE3293"),
                m.shape = c(22, 21),
                min.trees = 10 )

#Plot count PIS node 2
plot.filterNode(anomaly.zone.data = anomaly.data,
                concordance.factors.data = concord.data,
                output.dir = "Filter-Plots",
                focal.node = "node2",
                filter.name = "count_pis",
                dataset.name = "all",
                plot.gcf = TRUE,
                plot.scf = TRUE,
                az.colors = c("#7BC143", "#DE3293"),
                m.shape = c(22, 21),
                min.trees = 10 )

#Plot count PIS  node 6
plot.filterNode(anomaly.zone.data = anomaly.data,
                concordance.factors.data = concord.data,
                output.dir = "Filter-Plots",
                focal.node = "node6",
                filter.name = "count_pis",
                dataset.name = "all",
                plot.gcf = TRUE,
                plot.scf = TRUE,
                az.colors = c("#7BC143", "#DE3293"),
                m.shape = c(22, 21),
                min.trees = 10 )

###############################################################################
###############################################################################
###############  PLOT THE NODES OF INTEREST!           ########################
###############################################################################
###############################################################################

#Set up your directories
align.dir = "/Volumes/Armored/Hylidae/Alignments"
species.tree.path = "/Users/chutter/Dropbox/Research/1_Main-Projects/0_Working-Projects/Hylidae/Trees/Astral/Astral_everything.tre"
species.tree.path = "/Users/chutter/Dropbox/Research/1_Main-Projects/0_Working-Projects/Hylidae/Trees/Concat/all-markers_trimmed_50.phy.tre"
tree.dir = "/Users/chutter/Dropbox/Research/1_Main-Projects/0_Working-Projects/Hylidae/Trees/Gene_Trees/all-markers_trimmed"
work.dir = "/Volumes/Armored/Hylidae/Dataset-filtering"
setwd(work.dir)

align.summary = read.csv("master_alignment_summary.csv")
filt.summary = read.csv("master_filter_summary.csv")

outgroups = c("Phyllomedusa_tomopterna_WED_55506", "Nyctimystes_infrafrenatus_SLT_771")


### discordance versus alignment length

erroneousZone = function(species.tree.file = NULL,
                         outgroup.taxa = NULL,
                         genetree.directory = NULL,
                         alignment.summary = align.summary,
                         low.support.value = 94,
                         remove.taxa = NULL,
                         output.dir = "Erroneous-Zone",
                         dataset.name = "all",
                         filter.name = c("alignment_length", "count_pis", "proportion_pis", "proportion_sample"),
                         zone.colors = c("#7BC143", "#DE3293") ) {

  #Debug
  remove.taxa = c("Phyllomedusa_tomopterna_WED_55506", "Nyctimystes_infrafrenatus_SLT_771")
  low.support.value = 94
  outgroup.taxa = outgroups
  species.tree.file = species.tree.path
  genetree.directory = tree.dir
  alignment.summary = align.summary
  alignment.summary = alignment.summary[alignment.summary$dataset %in% "all-markers_trimmed",]
  output.dir = "Erroneous-Zone"
  dataset.name = "all"
  filter.name = c("alignment_length", "count_pis", "proportion_pis", "proportion_sample")
  zone.colors = c("#7BC143", "#DE3293")

  #nitial checks
  if (is.null(dataset.name) == TRUE){ stop("Error: a dataset name is needed.")}
  if (is.null(filter.name) == TRUE){ stop("Error: filter choice in filter.name is needed.")}
  #if (length(filter.name) != 1){ stop("Error: Only 1 filter can be plotted at a time.")}
  if (is.null(dataset.name) == TRUE){ stop("Error: a dataset name is needed.")}

  #Create working directory
  if(file.exists(output.dir) == F) { dir.create(output.dir) }

  #Get  trees
  gene.trees = list.files(genetree.directory)

  #Read in species tree
  species.tree = ape::read.tree(species.tree.file)
  species.tree = ape::drop.tip(species.tree, remove.taxa)
  species.tree = ape::unroot(species.tree)

  #CReates empty datasets
  header.data = c("file", "marker_name", "rf_dist_abs", "rf_dist_norm", "tree_discordant",
                  "tree_n_samples", "mean_support", "min_support", "max_support",
                  "n_low_support", "prop_low_support",
                  "align_n_samples", "align_proportion_samples", "alignment_length",
                  "alignment_count_pis", "alignment_proportion_pis", "alignment_count_missing_bp",
                  "alignment_proportion_missing_bp")

  collect.data = data.table(matrix(as.numeric(0),
                                   nrow = length(gene.trees),
                                   ncol = length(header.data)))
  setnames(collect.data, header.data)
  collect.data[, file:=as.character(file)]
  collect.data[, marker_name:=as.character(marker_name)]

  index.val = 1
  for (i in 1:length(gene.trees)){

    #Read in gene tree
    tree.gene = ape::read.tree(paste0(genetree.directory, "/", gene.trees[i]) )
    tree.gene = ape::drop.tip(tree.gene, remove.taxa)
    tree.gene = ape::unroot(tree.gene)
    gene.name = gsub(".phy.*", "", gene.trees[i])

    #Drop taxa from each to make the same tree
    drop.taxa = species.tree$tip.label[!species.tree$tip.label %in% tree.gene$tip.label]
    red.tree = ape::drop.tip(species.tree, drop.taxa)
    drop.taxa = tree.gene$tip.label[!tree.gene$tip.label %in% species.tree$tip.label]
    tree.gene = ape::drop.tip(tree.gene, drop.taxa)
    #Skips if removed too much
    if (length(tree.gene$tip.label) <= 3){ next }

    #Obtains alignment summary data
    gene.summary = alignment.summary[gsub("\\..*", "", alignment.summary$file) %in% gene.name,]

    #Gets gene tree stats
    gene.node.mean = mean(as.numeric(tree.gene$node.label), na.rm = T)
    gene.node.low = as.numeric(tree.gene$node.label)[as.numeric(tree.gene$node.label) <= low.support.value]
    gene.node.low = length(gene.node.low[is.na(gene.node.low) != T])
    gene.node.min = min(as.numeric(tree.gene$node.label), na.rm = T)
    gene.node.max = max(as.numeric(tree.gene$node.label), na.rm = T)

    #Gets rf distances
    gene.rf.norm = phangorn::RF.dist(red.tree, tree.gene, normalize = T)
    gene.rf.abs = phangorn::RF.dist(red.tree, tree.gene, normalize = F)
    discord.value = 0
    if (gene.rf.norm != 0){ discord.value = 1 }

    #Saves data
    set(collect.data, i = as.integer(index.val), j = match("file", header.data), value = gene.trees[i] )
    set(collect.data, i = as.integer(index.val), j = match("marker_name", header.data), value = gene.name )
    set(collect.data, i = as.integer(index.val), j = match("rf_dist_abs", header.data), value = gene.rf.abs )
    set(collect.data, i = as.integer(index.val), j = match("rf_dist_norm", header.data), value = gene.rf.norm )
    set(collect.data, i = as.integer(index.val), j = match("tree_discordant", header.data), value = discord.value )
    set(collect.data, i = as.integer(index.val), j = match("tree_n_samples", header.data), value = length(tree.gene$tip.label) - length(remove.taxa) )
    set(collect.data, i = as.integer(index.val), j = match("mean_support", header.data), value = gene.node.mean )
    set(collect.data, i = as.integer(index.val), j = match("min_support", header.data), value = gene.node.min )
    set(collect.data, i = as.integer(index.val), j = match("max_support", header.data), value = gene.node.max )
    set(collect.data, i = as.integer(index.val), j = match("n_low_support", header.data), value = gene.node.low )
    set(collect.data, i = as.integer(index.val), j = match("prop_low_support", header.data), value = gene.node.low/(length(tree.gene$node.label)-1) )
    set(collect.data, i = as.integer(index.val), j = match("align_n_samples", header.data), value = gene.summary$number_samples )
    set(collect.data, i = as.integer(index.val), j = match("align_proportion_samples", header.data), value = gene.summary$proportion_samples )
    set(collect.data, i = as.integer(index.val), j = match("alignment_length", header.data), value = gene.summary$alignment_length )
    set(collect.data, i = as.integer(index.val), j = match("alignment_count_pis", header.data), value = gene.summary$count_pis )
    set(collect.data, i = as.integer(index.val), j = match("alignment_proportion_pis", header.data), value = gene.summary$proportion_pis )
    set(collect.data, i = as.integer(index.val), j = match("alignment_count_missing_bp", header.data), value = gene.summary$count_missing_bp )
    set(collect.data, i = as.integer(index.val), j = match("alignment_proportion_missing_bp", header.data), value = gene.summary$proportion_missing_bp )

    #Adds oup counter
    index.val = index.val + 1
  }#end i loop

  #All regression
  summary(lm((collect.data$rf_dist_norm) ~ collect.data$alignment_length ))


  #concordant = 0
  c.data = collect.data[collect.data$rf_dist_norm == 0,]
  c.data$tree_discordant = factor(c.data$tree_discordant)
  c.data = c.data[is.na(c.data$tree_discordant) != T,]

  #Discordant = 1
  d.data = collect.data[collect.data$rf_dist_norm >= 0.5,]
  d.data$tree_discordant = factor(d.data$tree_discordant)
  d.data = d.data[is.na(d.data$tree_discordant) != T,]

  #Box plot and t-test
  boxplot(c.data$alignment_length, d.data$alignment_length)
  og.model = t.test(c.data$alignment_length, d.data$alignment_length, var.equal = F)
  og.model
  #cohensD(red.orig$rfDist, red.exon$rfDist)

  plot.data = collect.data
  plot.data = plot.data[plot.data$tree_n_samples >= 21,]
  plot.data$tree_discordant = factor(plot.data$tree_discordant)
  ggplot(plot.data, aes(x=alignment_length, color = tree_discordant)) +
    geom_density()




  boxplot(plot.data$alignment_length~plot.data$tree_discordant)
  model = glm(Discord ~ alignment_length, data = plot.data, family = "binomial")
  summary(model)
  exp(cbind(Odds=coef(model), confint(model)))


}#END FUNCTION




  #Loop through each gene tree and compare exon trees
  index.val = 1
  for (i in 1:length(gene.trees)){

    gene.name = gsub(".phy.*", "", gene.trees[i])
    gene.exons = exon.trees[grep(gene.name, exon.trees)]
    tree.gene = read.tree(paste0(gene.tree.dir, "/", gene.trees[i]) )
    tree.gene = drop.tip(tree.gene, outgroups)
    if (length(tree.gene$tip.label) <= 3){ next }

    gene.node.mean = mean(as.numeric(tree.gene$node.label), na.rm = T)
    gene.low.node = as.numeric(tree.gene$node.label)[as.numeric(tree.gene$node.label) <= 94]
    gene.low.node = gene.low.node[is.na(gene.low.node) != T]

    drop.taxa = species.tree$tip.label[!species.tree$tip.label %in% tree.gene$tip.label]
    red.tree = drop.tip(species.tree, drop.taxa)

    red.tree = unroot(red.tree)
    gene.dist = RF.dist(red.tree, tree.gene, normalize = T)

    if (length(gene.exons) == 0){ next }

    #Loops through each exon
    for (j in 1:length(gene.exons)){

      tree.exon = read.tree(paste0(exon.tree.dir, "/", gene.exons[j]) )
      tree.exon = drop.tip(tree.exon, outgroups)
      if (length(tree.exon$tip.label) <= 3){ next }

      drop.taxa = tree.gene$tip.label[!tree.gene$tip.label %in% tree.exon$tip.label]
      red.tree = drop.tip(tree.gene, drop.taxa)

      #Drops taxa not on the exon tree but on the gene tree
      drop.taxa = tree.exon$tip.label[!tree.exon$tip.label %in% red.tree$tip.label]
      red.exon = drop.tip(tree.exon, drop.taxa)
      exon.dist = RF.dist(red.tree, red.exon, normalize = T)

      #Gets bootstrap
      exon.node.mean = mean(as.numeric(tree.exon$node.label), na.rm = T)
      exon.low.node = as.numeric(tree.exon$node.label)[as.numeric(tree.exon$node.label) <= 94]
      exon.low.node = exon.low.node[is.na(exon.low.node) != T]

      #Gets species tree data
      drop.taxa = species.tree$tip.label[!species.tree$tip.label %in% tree.exon$tip.label]
      red.species = drop.tip(species.tree, drop.taxa)
      red.species = unroot(red.species)
      spp.dist = RF.dist(red.species, tree.exon, normalize = T)

      #Saves data
      exon.name = gsub(".phy.*", "", gene.exons[j])
      set(collect.data, i = as.integer(index.val), j = match("Exon", header.data), value = exon.name )
      set(collect.data, i = as.integer(index.val), j = match("Gene", header.data), value = gene.name )
      set(collect.data, i = as.integer(index.val), j = match("nSamples", header.data), value = length(red.species$tip.label) )
      set(collect.data, i = as.integer(index.val), j = match("rfDist_Exon_Gene", header.data), value = exon.dist )
      set(collect.data, i = as.integer(index.val), j = match("rfDist_Exon_Spp", header.data), value = spp.dist )
      set(collect.data, i = as.integer(index.val), j = match("rfDist_Gene_Spp", header.data), value = gene.dist )
      set(collect.data, i = as.integer(index.val), j = match("Gene_lowBS_prop", header.data), value = length(gene.low.node)/(length(tree.gene$node.label)-1) )
      set(collect.data, i = as.integer(index.val), j = match("Gene_meanBS", header.data), value = gene.node.mean )
      set(collect.data, i = as.integer(index.val), j = match("Exon_lowBS_prop", header.data), value = length(exon.low.node)/(length(tree.exon$node.label)-1) )
      set(collect.data, i = as.integer(index.val), j = match("Exon_meanBS", header.data), value = exon.node.mean )

      index.val = index.val + 1

    } #end j loop

  } #end i loop






  ###############



  #Start function
  concord.red = concordance.factors.data[concordance.factors.data$clade %in% focal.node,]
  filter.datasets = unique(concord.red$dataset)

  collect.data = c()
  for (x in 1:length(filter.datasets)){

    temp.a = anomaly.zone.data[anomaly.zone.data$filter_file %in% filter.datasets[x],]
    temp.c = concord.red[concord.red$dataset %in% filter.datasets[x],]
    mono.p = unique(temp.c$monophyletic)
    mono.p = mono.p[is.na(mono.p) != T]
    temp.c = temp.c[temp.c$gCF %in% min(temp.c$gCF),]
    temp.c$monophyletic = mono.p

    if (length(unique(temp.a$anomaly_zone)) == 1){
      temp.collect = cbind(temp.c, temp.a[1,])
      collect.data = rbind(collect.data, temp.collect)
      next
    }

    #Match up anomaly zone or not
    node.data1 = temp.a[temp.a$parent_node %in% temp.c$node,]
    node.data2 = temp.a[temp.a$child_node %in% temp.c$node,]
    node.data3 = temp.a[temp.a$parent_node %in% node.data1$child_node,]
    node.data = rbind(node.data1, node.data2, node.data3)

    az.data = node.data[node.data$anomaly_zone == 1,]

    #Picks the best
    if (nrow(az.data) >= 2){ az.data = az.data[1,] }
    if (nrow(az.data) == 0){ az.data = node.data[1,] }

    temp.collect = cbind(temp.c, az.data)
    collect.data = rbind(collect.data, temp.collect)

  }#x loop end

  ###################################################
  #Filters down to the target dataset and filter set
  ###################################################
  sub.data = collect.data[collect.data$filter_name %in% filter.name,]

  if (filter.name == "alignment_length"){ sub.data$final_filter = sub.data$filter_length }
  if (filter.name == "count_pis"){ sub.data$final_filter = sub.data$filter_count_pis }
  if (filter.name == "proportion_pis"){ sub.data$final_filter = sub.data$filter_prop_pis }
  if (filter.name == "proportion_sample"){ sub.data$final_filter = sub.data$filter_sample }

  if (dataset.name != "all"){
    sub.data = sub.data[sub.data$align_dataset %in% dataset.name,]
    dataset.names = dataset.name
  } else{
    dataset.names = unique(sub.data$align_dataset)
  }

  plot.list = vector("list", length(dataset.names))
  for (x in 1:length(dataset.names)){
    red.data = sub.data[sub.data$align_dataset %in% dataset.names[x],]

    if (plot.gcf == TRUE){
      plot.data = data.frame(Type = "gCF",
                             CF = red.data$gCF,
                             Trees = red.data$no_trees,
                             Filter = red.data$final_filter,
                             AZ = red.data$anomaly_zone,
                             M = red.data$monophyletic)
    }
    if (plot.scf == TRUE){
      plot.data1 = data.frame(Type = "sCF",
                              CF = red.data$sCF,
                              Trees = red.data$no_trees,
                              Filter = red.data$final_filter,
                              AZ = red.data$anomaly_zone,
                              M = red.data$monophyletic)
      plot.data = rbind(plot.data, plot.data1)
    }

    #plot.data = rbind(plot1.data, plot2.data)
    plot.data = plot.data[plot.data$Trees > min.trees,]
    plot.data = plot.data[plot.data$Filter != 0,]
    plot.data = plot.data[is.na(plot.data$CF) != T,]
    plot.data$Filter = factor(plot.data$Filter)
    plot.data$M = factor(plot.data$M)
    plot.data$AZ = factor(plot.data$AZ)

    #Sets up colors in case there is just one
    if (length(unique(plot.data$AZ)) == 1){

      if(unique(plot.data$AZ) == 1){
        temp.az.colors = az.colors[2]
      }

      if(unique(plot.data$AZ) == 0){
        temp.az.colors = az.colors[1]
      }
    } else {temp.az.colors = az.colors }

    #Sets up shapes in case there is just one
    if (length(unique(plot.data$M)) == 1){

      if(unique(plot.data$M) == TRUE){
        temp.m.shape = m.shape[2]
      }

      if(unique(plot.data$M) == FALSE){
        temp.m.shape = m.shape[1]
      }
    } else { temp.m.shape = m.shape }


    #Have to display 1 state correctly for each of the four
    #Pick a better shape?

    #Plot gene concordance factor
    if (plot.gcf == TRUE){
      gcf.plot = plot.data[plot.data$Type == "gCF",]

      #Getst the best y limit
      best.y = signif(max(gcf.plot$CF), digits = 0)
      if (best.y >= 80){
        best.y = best.y + 10
        #If greater than 100
        if(best.y >= 100){ best.y = 100 }
      } else { best.y = 80 }
      #End the y limit

      #Final one for paper
      main.plot = ggplot2::ggplot(gcf.plot, ggplot2::aes(x=Filter, y=CF, fill = AZ, shape = M)) +
        ggplot2::geom_point(ggplot2::aes(color = AZ, shape = M), size = 10, color = "black") +
        ggplot2::ylim(limits=c(0, best.y)) +
        ggplot2::scale_shape_manual(values = temp.m.shape) +
        ggplot2::scale_fill_manual(values=temp.az.colors)

      main.plot = main.plot +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.title = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank())

      main.plot

      #Saves the plot
      if(is.null(save.plots) != TRUE) {
        ggplot2::ggsave(paste0(output.dir, "/GCF_", focal.node, "_", filter.name, "_", dataset.names[x], ".pdf"),
                        width = 8, height = 4, device = "pdf")
        print(paste0(output.dir, "/GCF_", focal.node, "_", filter.name, "_", dataset.names[x], ".pdf",
                     " saved to file!"))
      }#end if

    }#end GCF plot

    #Plot site concordance factor
    if (plot.scf == TRUE){
      scf.plot = plot.data[plot.data$Type == "sCF",]

      #Getst the best y limit
      best.y = signif(max(scf.plot$CF), digits = 0)
      if (best.y >= 80){
        best.y = best.y + 10
        #If greater than 100
        if(best.y >= 100){ best.y = 100 }
      } else { best.y = 80 }
      #End the y limit

      #Final one for paper
      main.plot = ggplot2::ggplot(scf.plot, ggplot2::aes(x=Filter, y=CF, fill = AZ, shape = M)) +
        ggplot2::geom_point(ggplot2::aes(color = AZ, shape = M), size = 10, color = "black") +
        ggplot2::ylim(limits=c(0, best.y)) +
        ggplot2::scale_shape_manual(values = temp.m.shape) +
        ggplot2::scale_fill_manual(values= temp.az.colors)

      main.plot = main.plot +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.title = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank())

      main.plot
      plot.list[[x]] = main.plot

      #Saves the plot
      if(is.null(save.plots) != TRUE) {
        ggplot2::ggsave(paste0(output.dir, "/SCF_", focal.node, "_", filter.name, "_", dataset.names[x], ".pdf"),
                        width = 8, height = 4, device = "pdf")    }#end if
      print(paste0(output.dir, "/SCF_", focal.node, "_", filter.name, "_", dataset.names[x], ".pdf",
                   " saved to file!"))
    }#end GCF plot

  }#end dataset loop

  #### COMBINE INTO 1 FIGURE

  #ggpubr::ggarrange(plot.list[[1]], plot.list[[2]], plot.list[[3]], plot.list[[4]], plot.list[[5]], nrow = 2)
  #ggsave(paste(output.dir, "/GCF_", focal.node, "_", filter.name, "_ALL.pdf"), width = 28, height = 4, device = "pdf")

}#end function





### discordance versus alignment length after filtration




###############################################################################
###############################################################################
###############  PLOT THE NODES OF INTEREST!           ########################
###############################################################################
###############################################################################

#### Ideas: Randomization test
#Using concordance factors, we develop an alternative test of whether the species tree
#is in the erroneous zone, predicting that gCF will be less than sCF because EGTs will be
#unable estimate the true species tree because of insufficient informative sites.
#Conversely, to support the anomaly zone, we would expect gCF and sCF to be about equal
#because ILS would occur on the underlying sites and be reflected on gene trees given accurate
#estimation.



