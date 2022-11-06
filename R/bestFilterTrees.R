#' @title bestFilterTrees
#'
#' @description Function for saving the best trees that resulted from dataset filtering
#'
#' @param anomaly.zone.data table output from the filterAnomalies function
#'
#' @param concordance.factors.data table output from the filterConcordance function
#'
#' @param output.dir the name of the output directory to save the plots if TRUE above
#'
#' @param single.best saves a pdf of the single best tree
#'
#' @param top.best saves the top five best trees
#'
#' @param all.datasets TRUE for top five for each datasets FALSE for top five across all datasets
#'
#' @param fewest.anomaly.zones TRUE to return the trees with the fewest anomaly zones
#'
#' @param highest.post.prob TRUE to return the trees with the highest mean posterior probability
#'
#' @param highest.gene.cf TRUE to return the trees with the highest mean gene concordance factors
#'
#' @param min.trees minimum number of trees to keep a filtration replicate. Default: 10
#'
#' @return a folder saved as output.dir that contains the options for your best trees
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

bestFilterTrees = function(anomaly.zone.data = anomaly.data,
                           concordance.factors.data = concord.data,
                           output.dir = NULL,
                           single.best = TRUE,
                           top.best = 5,
                           all.datasets = TRUE,
                           fewest.anomaly.zones = TRUE,
                           highest.post.prob = TRUE,
                           highest.gene.cf = TRUE,
                           min.trees = 10) {

  #Debug
  # anomaly.zone.data = anomaly.data
  # concordance.factors.data = concord.data
  # output.dir = "best"
  # min.trees = 20
  # fewest.anomaly.zones = TRUE
  # highest.gene.cf = TRUE
  # highest.post.prob = TRUE
  # all.datasets = TRUE
  # top.best = 1

  #parameter checks
  if(is.null(anomaly.zone.data) == TRUE){ stop("Error: no anomaly.zone.data provided.") }
  if(is.null(concordance.factors.data) == TRUE){ stop("Error: no concordance.factors.data provided.") }
  if(is.null(output.dir) == TRUE){ stop("Error: no output save directory name provided.") }

  #nitial checks
  if(file.exists(output.dir) == F) { dir.create(output.dir) }

  #Filters data further
  anomaly.red = anomaly.zone.data[anomaly.zone.data$no_trees > min.trees,]
  if (nrow(anomaly.red) == 0){ stop("Error: no datasets greater than set minimum trees value.") }
  concord.red = concordance.factors.data[concordance.factors.data$dataset %in% anomaly.red$filter_file,]
  if (nrow(concord.red) == 0){ stop("Error: no datasets greater than set minimum trees value.") }

  filter.datasets = unique(concord.red$dataset)

  collect.data = c()
  for (x in 1:length(filter.datasets)){

    temp.a = anomaly.zone.data[anomaly.zone.data$filter_file %in% filter.datasets[x],]
    temp.c = concord.red[concord.red$dataset %in% filter.datasets[x],]

    #Data to collect:
    #Have: anomaly zone data and concord data for every tree
    #Summary of trees
    # - Count AZ
    # - get average gCF / sCF?
    # - min, max gCF / sCF
    # - gather filter information, do this to everything at once?
    # - gather alignment mean data too

    filter.name = as.character(unique(temp.a$filter_name))
    #Gets filter value
    if (filter.name == "alignment_length"){ filter.value = mean(temp.a$filter_length, na.rm = T) }
    if (filter.name == "count_pis"){ filter.value = mean(temp.a$filter_count_pis, na.rm = T) }
    if (filter.name == "proportion_pis"){ filter.value = mean(temp.a$filter_prop_pis, na.rm = T) }
    if (filter.name == "proportion_samples"){ filter.value = mean(temp.a$filter_sample, na.rm = T) }
    #Gets alignment value
    if (filter.name == "alignment_length"){ align.value = mean(temp.a$mean_length, na.rm = T) }
    if (filter.name == "count_pis"){ align.value = mean(temp.a$mean_count_pis, na.rm = T) }
    if (filter.name == "proportion_pis"){ align.value = mean(temp.a$mean_prop_pis, na.rm = T) }
    if (filter.name == "proportion_samples"){ align.value = mean(temp.a$mean_sample, na.rm = T) }

    #Summarizes anomaly zone and alignment data
    az.data = data.frame(dataset = unique(temp.a$align_dataset),
                         file_name = filter.datasets[x],
                         filter_name = filter.name,
                         filter_value = filter.value,
                         alignment_value = align.value,
                         no_trees = unique(temp.a$no_trees),
                         count_AZ = sum(temp.a$anomaly_zone, na.rm = T),
                         prop_AZ = sum(temp.a$anomaly_zone, na.rm = T)/nrow(temp.a))

    #Summarizes concordance factors data
    cf.data = data.frame(mean_gCF = mean(temp.c$gCF, na.rm = T),
                         min_gCF = min(temp.c$gCF, na.rm = T),
                         mean_sCF = mean(temp.c$sCF, na.rm = T),
                         min_sCF = min(temp.c$sCF, na.rm = T),
                         mean_pp = mean(temp.c$pp1, na.rm = T),
                         min_pp = min(temp.c$pp1, na.rm = T))

    temp.collect = cbind(az.data, cf.data)
    collect.data = rbind(collect.data, temp.collect)

  }#x loop end

  ### Collect the data
  ###########################

  save.data.table = c()
  if (fewest.anomaly.zones == TRUE){
    #to include one data or none
    if (all.datasets == FALSE){
      fewest.az = collect.data[order(collect.data$count_AZ),]
      fewest.az = fewest.az[1:top.best,]
      fewest.az = cbind(best_tree_category = "fewest-anomaly-zone", fewest.az)
      #Save the data
      save.data.table = rbind(save.data.table, fewest.az)
      #Save the tree
      for (y in 1:nrow(fewest.az)){
        system(paste0("cp ", "filtered-astral/", fewest.az$file_name[y], "_astral.tre ",
                      output.dir, "/", fewest.az$file_name[y], "_fewest-az.tre"))
      }#end y loop
    } else {
      #Gets all the datasets instead of just the best
      dataset.name = unique(collect.data$dataset)
      for (i in 1:length(dataset.name)){
        fewest.az = collect.data[collect.data$dataset %in% dataset.name[i],]
        fewest.az = fewest.az[order(fewest.az$count_AZ),]
        fewest.az = fewest.az[1:top.best,]
        fewest.az = cbind(best_tree_category = "fewest-anomaly-zone", fewest.az)
        #Save the tree
        for (y in 1:nrow(fewest.az)){
          system(paste0("cp ", "filtered-astral/", fewest.az$file_name[y], "_astral.tre ",
                        output.dir, "/", fewest.az$file_name[y], "_fewest-az.tre"))
        }#end y loop

        #Save the data
        save.data.table = rbind(save.data.table, fewest.az)
      }#end i loop
    }#end else
  }#fewest az if

  #Highest posterior probability
  if (highest.post.prob == TRUE){

    if (all.datasets == FALSE){
      #Orders and then gets top best
      highest.pp = collect.data[order(collect.data$mean_pp, decreasing = F),]
      highest.pp = highest.pp[1:top.best,]
      highest.pp = cbind(best_tree_category = "highest-posterior-probability", highest.pp)
      #Save the tree
      for (y in 1:nrow(highest.pp)){
        system(paste0("cp ", "filtered-astral/", highest.pp$file_name[y], "_astral.tre ",
                      output.dir, "/", highest.pp$file_name[y], "_highest-pp.tre"))
      }#end y loop
      #Save the data
      save.data.table = rbind(save.data.table, highest.pp)
    } else {
      #Gets the top best count instead of just one
      dataset.name = unique(collect.data$dataset)
      for (i in 1:length(dataset.name)){
        highest.pp = collect.data[collect.data$dataset %in% dataset.name[i],]
        #Orders and then gets top best
        highest.pp = highest.pp[order(highest.pp$mean_pp, decreasing = F),]
        highest.pp = highest.pp[1:top.best,]
        highest.pp = cbind(best_tree_category = "highest-posterior-probability", highest.pp)
        #Save the tree
        for (y in 1:nrow(highest.pp)){
          system(paste0("cp ", "filtered-astral/", highest.pp$file_name[y], "_astral.tre ",
                        output.dir, "/", highest.pp$file_name[y], "_highest-pp.tre"))
        }#end y loop
        #Save the data
        save.data.table = rbind(save.data.table, highest.pp)
      }#end i loop
    }#end else
  }#end if

  #Highest gene concordance
  if (highest.gene.cf == TRUE){

    if (all.datasets == FALSE){
      #Orders and then gets top best
      highest.gcf = collect.data[order(collect.data$mean_gCF, decreasing = T),]
      highest.gcf = highest.gcf[1:top.best,]
      highest.gcf = cbind(best_tree_category = "highest-gene-concordance", highest.gcf)
      #Save the tree
      for (y in 1:nrow(highest.gcf)){
        system(paste0("cp ", "filtered-astral/", highest.gcf$file_name[y], "_astral.tre ",
                      output.dir, "/", highest.gcf$file_name[y], "_highest-gcf.tre"))
      }#end y loop
      #Save the data
      save.data.table = rbind(save.data.table, highest.gcf)
    } else {
      #Gets the top best count instead of just one
      dataset.name = unique(collect.data$dataset)
      for (i in 1:length(dataset.name)){
        highest.gcf = collect.data[collect.data$dataset %in% dataset.name[i],]
        #Orders and then gets top best
        highest.gcf = highest.gcf[order(highest.gcf$mean_gCF, decreasing = T),]
        highest.gcf = highest.gcf[1:top.best,]
        highest.gcf = cbind(best_tree_category = "highest-gene-concordance", highest.gcf)
        #Save the tree
        for (y in 1:nrow(highest.gcf)){
          system(paste0("cp ", "filtered-astral/", highest.gcf$file_name[y], "_astral.tre ",
                        output.dir, "/", highest.gcf$file_name[y], "_highest-gcf.tre"))
        }#end y loop
        #Save the data
        save.data.table = rbind(save.data.table, highest.gcf)
      }#end i loop
    }#end else

  }#end if

  write.table(save.data.table,
              file = paste0(output.dir, "/summary-best-datasets.txt"),
              sep = "\t",
              row.names = F)

}#end function
