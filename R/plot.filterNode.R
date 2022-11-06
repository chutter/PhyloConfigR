#' @title plot.filterNode
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


plot.filterNode = function(anomaly.zone.data = NULL,
                           concordance.factors.data = NULL,
                           save.plots = TRUE,
                           output.dir = "Filter-Plots",
                           focal.node = NULL,
                           filter.name = c("alignment_length", "count_pis", "proportion_pis", "proportion_sample"),
                           dataset.name = NULL,
                           plot.gcf = TRUE,
                           plot.scf = TRUE,
                           az.colors = c("#7BC143", "#DE3293"),
                           m.shape = c(22, 21),
                           min.trees = 10){

  #nitial checks
  if (is.null(dataset.name) == TRUE){ stop("Error: a dataset name is needed.")}
  if (is.null(filter.name) == TRUE){ stop("Error: filter choice in filter.name is needed.")}
  if (length(filter.name) != 1){ stop("Error: Only 1 filter can be plotted at a time.")}
  if (is.null(dataset.name) == TRUE){ stop("Error: a dataset name is needed.")}
  if (is.null(focal.node) == TRUE){ stop("Error: a focal node is needed.")}

  if(file.exists(output.dir) == F) { dir.create(output.dir) }

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


