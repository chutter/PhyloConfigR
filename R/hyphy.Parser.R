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
                      hyphy.analysis = c("absrel", "BUSTED"),
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
  hyphy.analysis = "SLAC"
  threads = 4
  memory = 8
  resume = T
  overwrite = F
  quiet = T
  tips.only = T
  hyphy.path = "/usr/local/bin/"


  #Same adds to bbmap path
  if (is.null(hyphy.path) == FALSE){
    b.string = unlist(strsplit(hyphy.path, ""))
    if (b.string[length(b.string)] != "/") {
      hyphy.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { hyphy.path = "" }

  if (is.null(results.directory) == T){ stop("A directory of results is needed.") }

  results.files = list.dirs(results.directory, full.names = F)
  results.files = results.files[results.files != ""]

  #################################################################################
  ################################# SLAC ############################
  ######################################################
  if (hyphy.analysis == "SLAC"){

    #Stats table prepare
    #header.data = c("sample", "constrained", "unconstrained", "gtr_model", "mg94xrev")
    #collect.data = data.table::data.table(matrix(as.numeric(0), nrow = length(results.files), ncol = length(header.data)))
    #data.table::setnames(collect.data, header.data)
    #collect.data[, sample:=as.character(sample)]

    #Run analyses through all results output folders
    all.data = data.frame()
    for (i in 1:length(results.files)){

      json.data = jsonlite::fromJSON(paste0(results.directory, "/", results.files[i], "/SLAC-results.json"))
      branch.att = json.data$`branch attributes`
      branch.att = unlist(branch.att[[2]])

      #Gets header info
      mle.data = json.data$MLE
      data.headers = mle.data$headers[,1]

      #Gets data
      bysite.resolved = mle.data$content$`0`$`by-site`$RESOLVED
      colnames(bysite.resolved) = data.headers
      bysite.average = mle.data$content$`0`$`by-site`$AVERAGED
      colnames(bysite.average) = data.headers

      bybr.resolved = mle.data$content$`0`$`by-branch`$RESOLVED
      colnames(bybr.resolved) = data.headers
      bybr.averaged = mle.data$content$`0`$`by-branch`$AVERAGED
      colnames(bybr.averaged) = data.headers

      taxa.names = mle.data$content$`0`$`by-branch`$NAMES

      ### * Want by-site average for alignment-based stats
     # data.headers = c("file", "sample", "ES", "EN", "S", "N", "PS",
    #                   "dS", "dN", "dN-dS", "P>1", "P<1", "branchLength")
    #  branch.data = data.frame(file = results.files[i], sample = taxa.names, bybr.averaged)
     # colnames(branch.data) = data.headers
      #branch.data$omega = branch.data$dN/branch.data$dS

      # Write the output headers
      # ES = Expected synonymous sites
      # EN = Expected non-synonymous sites
      # S = Inferred synonymous substitutions
      # N = Inferred non-synonymous substitutions

      ### * Want by-branch average for the branches
      data.headers = c("file", "sample", "ES", "EN", "S", "N", "PS",
                       "dS", "dN", "dN-dS", "P>1", "P<1", "branchLength")
      branch.data = data.frame(file = results.files[i], sample = taxa.names, bybr.averaged)
      colnames(branch.data) = data.headers
      branch.data$omega = branch.data$dN/branch.data$dS

      all.data = rbind(all.data, branch.data)

    } #end i loop

    #removes nodes
    if (tips.only == TRUE){
      final.data = all.data[grep("Node.*", all.data$sample, invert = T),]
    } else { final.data = all.data }

    write.csv(final.data, file = paste0(output.name, "_SLAC_by-branch.csv"), row.names = F, quote = F)
    print(paste0("Finished ", hyphy.analysis, " data summary!"))

  }#end busted if
  #################################################################################



  #################################################################################
  ################################# BUSTED ############################
  ######################################################
  if (hyphy.analysis == "BUSTED"){

    #Stats table prepare
    #header.data = c("sample", "constrained", "unconstrained", "gtr_model", "mg94xrev")
    #collect.data = data.table::data.table(matrix(as.numeric(0), nrow = length(results.files), ncol = length(header.data)))
    #data.table::setnames(collect.data, header.data)
    #collect.data[, sample:=as.character(sample)]

    #Run analyses through all results output folders
    all.data = data.frame()
    for (i in 1:length(results.files)){

      json.data = jsonlite::fromJSON(paste0(results.directory, "/", results.files[i], "/BUSTED-results.json"))
      branch.att = json.data$`branch attributes`
      branch.att = unlist(branch.att[[1]])

      #Gather stats
      con.stats = branch.att[grep("\\.constrained", names(branch.att) )]
      uncon.stats = branch.att[grep("\\.unconstrained", names(branch.att) )]
      nucl.stats = branch.att[grep("\\.Nucleotide", names(branch.att) )]
      omega.stats = branch.att[grep("\\.MG94xREV", names(branch.att) )]

      #Make table from stats
      con.data = data.frame(Sample = gsub("\\.constrained", "", names(con.stats)),
                 null_constrained = as.numeric(con.stats))

      uncon.data = data.frame(Sample = gsub("\\.unconstrained", "", names(uncon.stats)),
                 unconstrained = as.numeric(uncon.stats))

      nucl.data = data.frame(Sample = gsub("\\.Nucleotide GTR", "", names(nucl.stats)),
                              gtr_model = as.numeric(nucl.stats))

      omega.data = data.frame(Sample = gsub("\\.MG94xREV.*", "", names(omega.stats)),
                              mg94xrev = as.numeric(omega.stats))

      if (nrow(con.data) != 0) { a.data = merge(con.data, uncon.data, by = "Sample") } else {
        a.data = uncon.data
        a.data$null_constrained = NA
      }
      b.data = merge(nucl.data, omega.data, by = "Sample")
      w.data = merge(a.data, b.data, by = "Sample")

      save.data = cbind(Locus = results.files[i], w.data)
      all.data = rbind(all.data, save.data)
    } #end i loop

  if (tips.only == TRUE){
    all.data = all.data[grep("Node.*", all.data$Sample, invert = T),]
  }

  write.csv(all.data, file = paste0(output.name, "_BUSTED-stats.csv"), row.names = F, quote = F)
  print(paste0("Finished ", hyphy.analysis, " data summary!"))

  }#end busted if
  #################################################################################


  #Directoires
  setwd("/Users/chutter/Dropbox/Research/1_Main-Projects/0_Working-Projects/Rodent_Mitochondrial")
  results.directory = "/Volumes/Rodents/Australian_Rodents/Data_Processing/hyphy/absrel_all"
  output.name = "absrel_all"
  hyphy.analysis = "absrel"
  threads = 4
  memory = 8
  resume = T
  overwrite = F
  quiet = T
  tips.only = T
  hyphy.path = "/usr/local/bin/"
  results.files = list.dirs(results.directory, full.names = F)
  results.files = results.files[results.files != ""]

  #################################################################################
  ################################# ABSRESL ############################
  ######################################################
  if (hyphy.analysis == "absrel"){

    #Stats table prepare
    #header.data = c("sample", "constrained", "unconstrained", "gtr_model", "mg94xrev")
    #collect.data = data.table::data.table(matrix(as.numeric(0), nrow = length(results.files), ncol = length(header.data)))
    #data.table::setnames(collect.data, header.data)
    #collect.data[, sample:=as.character(sample)]

    #Run analyses through all results output folders
    all.data = data.frame()
    for (i in 1:length(results.files)){

      json.data = jsonlite::fromJSON(paste0(results.directory, "/", results.files[i], "/absrel-results.json"))
      branch.att = json.data$`branch attributes`
      branch.att = unlist(branch.att[[1]])

      #Gather stats
      full.stats = branch.att[grep("\\.Full adaptive model$", names(branch.att) )]
      fs.stats = branch.att[grep("\\.Full adaptive model \\(synonymous subs/site\\)", names(branch.att) )]
      fns.stats = branch.att[grep("\\.Full adaptive model \\(non-synonymous subs/site\\)", names(branch.att) )]
      mg94.stats = branch.att[grep("\\.Baseline MG94xREV$", names(branch.att) )]
      omega.stats = branch.att[grep("\\.Baseline MG94xREV omega ratio", names(branch.att) )]
      upval.stats = branch.att[grep("\\.Uncorrected P-value", names(branch.att) )]
      cpval.stats = branch.att[grep("\\.Corrected P-value", names(branch.att) )]
      rate1.stats = branch.att[grep("\\.Rate Distributions1", names(branch.att) )]
      rate2.stats = branch.att[grep("\\.Rate Distributions2", names(branch.att) )]
      ngtr.stats = branch.att[grep("\\.Nucleotide GTR", names(branch.att) )]
      rcl.stats = branch.att[grep("\\.Rate classes", names(branch.att) )]
      lrt.stats = branch.att[grep("\\.LRT", names(branch.att) )]

      #Make table from stats
      full.data = data.frame(Sample = gsub("\\.Full adaptive model$", "", names(full.stats)),
                            full_adaptive_model = as.numeric(full.stats))

      fs.data = data.frame(Sample = gsub("\\.Full adaptive model \\(synonymous subs/site\\)", "", names(fs.stats)),
                              full_adaptive_model_syn = as.numeric(fs.stats))

      fns.data = data.frame(Sample = gsub("\\.Full adaptive model \\(non-synonymous subs/site\\)", "", names(fns.stats)),
                            full_adaptive_model_nonsyn = as.numeric(fns.stats))

      mg94.data = data.frame(Sample = gsub("\\.Baseline MG94xREV$", "", names(mg94.stats)),
                              baseline_mg94xrev = as.numeric(mg94.stats))

      omega.data = data.frame(Sample = gsub("\\.Baseline MG94xREV omega ratio", "", names(omega.stats)),
                            baseline_omega_ratio = as.numeric(omega.stats))

      upval.data = data.frame(Sample = gsub("\\.Uncorrected P-value", "", names(upval.stats)),
                              uncorrected_pval = as.numeric(upval.stats))

      cpval.data = data.frame(Sample = gsub("\\.Corrected P-value", "", names(cpval.stats)),
                             corrected_pval = as.numeric(cpval.stats))

      rate1.data = data.frame(Sample = gsub("\\.Rate Distributions1", "", names(rate1.stats)),
                              rate_distributions1 = as.numeric(rate1.stats))

      rate2.data = data.frame(Sample = gsub("\\.Rate Distributions2", "", names(rate2.stats)),
                              rate_distributions2 = as.numeric(rate2.stats))

      ngtr.data = data.frame(Sample = gsub("\\.Nucleotide GTR", "", names(ngtr.stats)),
                             nuclotide_gtr = as.numeric(ngtr.stats))

      rcl.data = data.frame(Sample = gsub("\\.Rate classes", "", names(rcl.stats)),
                              rate_classes = as.numeric(rcl.stats))

      lrt.data = data.frame(Sample = gsub("\\.LRT", "", names(lrt.stats)),
                              lrt = as.numeric(lrt.stats))


      # if (nrow(con.data) != 0) { a.data = merge(con.data, uncon.data, by = "Sample") } else {
      #   a.data = uncon.data
      #   a.data$null_constrained = NA
      # }
      a.data = merge(full.data, fs.data, by = "Sample")
      b.data = merge(a.data, fns.data, by = "Sample")
      c.data = merge(b.data, mg94.data, by = "Sample")
      d.data = merge(c.data, omega.data, by = "Sample")
      e.data = merge(d.data, upval.data, by = "Sample")
      f.data = merge(e.data, cpval.data, by = "Sample")
      g.data = merge(f.data, rate1.data, by = "Sample")
      h.data = merge(g.data, rate2.data, by = "Sample")
      i.data = merge(h.data, ngtr.data, by = "Sample")
      j.data = merge(i.data, rcl.data, by = "Sample")
      k.data = merge(j.data, lrt.data, by = "Sample")
      save.data = cbind(Locus = results.files[i], k.data)
      all.data = rbind(all.data, save.data)
    } #end i loop

    if (tips.only == TRUE){
      all.data = all.data[grep("Node.*", all.data$Sample, invert = T),]
    }

    write.csv(all.data, file = paste0(output.name, "_abrsel-stats.csv"), row.names = F, quote = F)
    print(paste0("Finished ", hyphy.analysis, " data summary!"))

  }#end busted if
  #################################################################################





}#end function


