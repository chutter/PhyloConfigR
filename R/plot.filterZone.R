#' @title plot.filterZone
#'
#' @description Function for plotting gene tree filtration results for a focal node
#'
#' @param anomaly.zone.data table output from the filterAnomalies function
#'
#' @param concordance.factors.data table output from the filterConcordance function
#'
#' @param save.plots if you wish to save to file select TRUE
#'
#' @param output.dir the name of the output directory to save the plots if TRUE above
#'
#' @param focal.node the node annotated from the filterConcordance function
#'
#' @param filter.name the filter to plot. Options include: alignment_length, count_pis, proportion_pis, proportion_sample
#'
#' @param dataset.name from your main sets of data (i.e. exons, introns). All = plots of all them.
#'
#' @param plot.gcf should the gene concordance factor be plotted?
#'
#' @param plot.scf should the site concordance factor be plotted?
#'
#' @param az.colors colors to indicate the anomaly zone. Default: Green: presence; Purple: absence
#'
#' @param m.shape monophyly shape on the graph; circle = monophyletic; square paraphyletic
#'
#' @param min.trees minimum number of trees to keep a filtration replicate. Default: 10
#'
#' @return plots a dot plot of each filter and the concordance factors for the focal node. Monophyly of that node is also shown.
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


plot.filterZone = function(anomaly.zone.data = NULL,
                           concordance.factors.data = NULL,
                           save.plots = TRUE,
                           output.dir = "Filter-Plots",
                           dataset.name = NULL,
                           plot.gcf = TRUE,
                           plot.scf = TRUE,
                           plot.pp = TRUE,
                           az.colors = c("#7BC143", "#DE3293"),
                           m.shape = c(22, 21),
                           min.trees = 10){
  #Debug
  # anomaly.zone.data = anomaly.data
  # concordance.factors.data = concord.data
  # output.dir = "Filter-Plots"
  # dataset.name = "100kb_regions"
  # plot.gcf = TRUE
  # plot.scf = TRUE
  # az.colors = c("#7BC143", "#DE3293")
  # m.shape = c(22, 21)
  # min.trees = 10
  # save.plots = TRUE

  if(is.null(anomaly.zone.data) == TRUE){ stop("Error: no anomaly.zone.data provided.") }
  if(is.null(concordance.factors.data) == TRUE){ stop("Error: no concordance.factors.data provided.") }
  if(is.null(output.dir) == TRUE){ stop("Error: no output save directory name provided.") }

  #nitial checks
  if(file.exists(output.dir) == F) { dir.create(output.dir) }

  #Start function
  if (dataset.name == "all" || is.null(dataset.name) == TRUE){
    concord.red = concordance.factors.data
  } else{
    concord.red = concordance.factors.data[grep(dataset.name, concordance.factors.data$dataset),]
  }
  if (nrow(concord.red) == 0){ stop("Error: dataset.name does not exist in data.") }

  #Filters data further
  anomaly.red = anomaly.zone.data[anomaly.zone.data$no_trees > min.trees,]
  if (nrow(anomaly.red) == 0){ stop("Error: no datasets greater than set minimum trees value.") }
  concord.red = concord.red[concord.red$dataset %in% anomaly.red$filter_file,]
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

  ###################################################
  #Plots separately
  ###################################################

  data.names = as.character(unique(collect.data$filter_name))

  plot.list = vector("list", length(data.names))
  for (x in 1:length(data.names)){
    red.data = collect.data[collect.data$filter_name %in% data.names[x],]

    if (plot.gcf == TRUE){
      plot.data = data.frame(Type = "gCF",
                             dataset = red.data$dataset,
                             CF = red.data$mean_gCF,
                             PP = red.data$mean_pp,
                             Trees = red.data$no_trees,
                             Filter = red.data$filter_value,
                             Align = red.data$alignment_value,
                             AZ = red.data$count_AZ)
    }
    if (plot.scf == TRUE){
      plot.data1 = data.frame(Type = "sCF",
                              dataset = red.data$dataset,
                              CF = red.data$mean_sCF,
                              PP = red.data$mean_pp,
                              Trees = red.data$no_trees,
                              Filter = red.data$filter_value,
                              Align = red.data$alignment_value,
                              AZ = red.data$count_AZ)
      plot.data = rbind(plot.data, plot.data1)
    }

    if (plot.pp == TRUE){
      plot.data1 = data.frame(Type = "pp",
                              dataset = red.data$dataset,
                              CF = red.data$mean_sCF,
                              PP = red.data$mean_pp,
                              Trees = red.data$no_trees,
                              Filter = red.data$filter_value,
                              Align = red.data$alignment_value,
                              AZ = red.data$count_AZ)
      plot.data = rbind(plot.data, plot.data1)
    }

    #plot.data = rbind(plot1.data, plot2.data)
    plot.data = plot.data[plot.data$Trees > min.trees,]
    plot.data = plot.data[plot.data$Filter != 0,]
    plot.data = plot.data[is.na(plot.data$CF) != T,]
    plot.data = plot.data[order(plot.data$Filter),]
    plot.data$Filter = factor(plot.data$Filter)

    #Sets up colors in case there is just one
    if (length(unique(plot.data$AZ)) == 1){

      if(unique(plot.data$AZ) == 1){
        temp.az.colors = az.colors[2]
      }

      if(unique(plot.data$AZ) == 0){
        temp.az.colors = az.colors[1]
      }
    } else {temp.az.colors = az.colors }


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
      main.plot = ggplot2::ggplot(gcf.plot, ggplot2::aes(x=Filter, y=CF, Fill = dataset)) +
        ggplot2::geom_point(size = 3, ggplot2::aes(color = dataset)) +
        ggplot2::ylim(limits=c(0, best.y))

      main.plot = main.plot +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.title = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank())

      main.plot

      #Saves the plot
      if(is.null(save.plots) != TRUE) {
        ggplot2::ggsave(paste0(output.dir, "/GCF_Filtered_", data.names[x], ".pdf"),
                        width = 8, height = 4, device = "pdf")
        print(paste0(output.dir, "/GCF_Filtered_", data.names[x], ".pdf",
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
      main.plot = ggplot2::ggplot(scf.plot, ggplot2::aes(x=Filter, y=CF, Fill = dataset)) +
        ggplot2::geom_point(size = 3, ggplot2::aes(color = dataset)) +
        ggplot2::ylim(limits=c(0, best.y))

      main.plot = main.plot +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.title = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank())

      main.plot

      plot.list[[x]] = main.plot

      #Saves the plot
      if(is.null(save.plots) != TRUE) {
        ggplot2::ggsave(paste0(output.dir, "/SCF_Filtered_", data.names[x], ".pdf"),
                        width = 8, height = 4, device = "pdf")    }#end if
      print(paste0(output.dir, "/SCF_Filtered_", data.names[x], ".pdf",
                   " saved to file!"))
    }#end GCF plot

    #Plot site concordance factor
    if (plot.pp == TRUE){
      pp.plot = plot.data[plot.data$Type == "pp",]

      #Getst the best y limit
      best.y = signif(max(pp.plot$PP), digits = 0)
      if (best.y >= .8){
        best.y = best.y + 1
        #If greater than 100
        if(best.y >= 1){ best.y = 1 }
      } else { best.y = .8 }
      #End the y limit

      #Final one for paper
      main.plot = ggplot2::ggplot(pp.plot, ggplot2::aes(x=Filter, y=PP, Fill = dataset)) +
        ggplot2::geom_point(size = 3, ggplot2::aes(color = dataset)) +
        ggplot2::ylim(limits=c(0, best.y))

      main.plot = main.plot +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.title = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank())

      main.plot

      plot.list[[x]] = main.plot

      #Saves the plot
      if(is.null(save.plots) != TRUE) {
        ggplot2::ggsave(paste0(output.dir, "/PP_Filtered_", data.names[x], ".pdf"),
                        width = 8, height = 4, device = "pdf")    }#end if
      print(paste0(output.dir, "/PP_Filtered_", data.names[x], ".pdf",
                   " saved to file!"))
    }#end GCF plot

  }#end dataset loop

  #### COMBINE INTO 1 FIGURE

  #ggpubr::ggarrange(plot.list[[1]], plot.list[[2]], plot.list[[3]], plot.list[[4]], plot.list[[5]], nrow = 2)
  #ggsave(paste(output.dir, "/GCF_", focal.node, "_", filter.name, "_ALL.pdf"), width = 28, height = 4, device = "pdf")

}#end function


