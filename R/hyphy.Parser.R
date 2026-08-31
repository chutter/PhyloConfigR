#' @title hyphy.Parser
#'
#' @description Parses HyPhy output JSON files (SLAC, BUSTED, or aBSREL) into
#'   summary tables. For SLAC, extracts by-branch dN/dS statistics. For BUSTED
#'   and aBSREL, extracts model fit and selection statistics.
#'
#' @param results.directory path to the directory of HyPhy results subfolders,
#'   each containing a JSON output file
#'
#' @param hyphy.analysis which HyPhy analysis to parse: "SLAC", "BUSTED", or "absrel"
#'
#' @param output.name prefix for the output CSV file name
#'
#' @param tips.only if TRUE, only tip branches are returned; internal nodes are excluded
#'
#' @param overwrite if TRUE overwrites existing output files
#'
#' @param quiet if TRUE suppresses progress messages
#'
#' @return writes a CSV summary table to file and returns it invisibly
#'
#' @examples
#'
#' hyphy.Parser(results.directory = "hyphy/SLAC",
#'              hyphy.analysis = "SLAC",
#'              output.name = "slac_results",
#'              tips.only = TRUE,
#'              overwrite = FALSE)
#'
#' @export

hyphy.Parser = function(results.directory = NULL,
                        hyphy.analysis = c("SLAC", "BUSTED", "absrel"),
                        output.name = NULL,
                        tips.only = TRUE,
                        overwrite = FALSE,
                        quiet = TRUE) {

  if (is.null(results.directory) == T){ stop("A directory of results is needed.") }
  if (is.null(output.name) == T){ stop("An output name is needed.") }
  if (length(hyphy.analysis) != 1){ stop("Please select exactly one hyphy.analysis type.") }

  results.files = list.dirs(results.directory, full.names = F)
  results.files = results.files[results.files != ""]

  #################################################################################
  ################################# SLAC ##########################################
  #################################################################################
  if (hyphy.analysis == "SLAC"){

    all.data = data.frame()
    for (i in 1:length(results.files)){

      json.file = paste0(results.directory, "/", results.files[i], "/SLAC-results.json")
      if (!file.exists(json.file)){ next }

      json.data = jsonlite::fromJSON(json.file)
      mle.data = json.data$MLE
      data.headers = mle.data$headers[,1]

      bybr.averaged = mle.data$content$`0`$`by-branch`$AVERAGED
      colnames(bybr.averaged) = data.headers
      taxa.names = mle.data$content$`0`$`by-branch`$NAMES

      data.headers = c("file", "sample", "ES", "EN", "S", "N", "PS",
                       "dS", "dN", "dN-dS", "P>1", "P<1", "branchLength")
      branch.data = data.frame(file = results.files[i], sample = taxa.names, bybr.averaged)
      colnames(branch.data) = data.headers
      branch.data$omega = branch.data$dN / branch.data$dS

      all.data = rbind(all.data, branch.data)
    }

    if (tips.only == TRUE){
      all.data = all.data[grep("Node.*", all.data$sample, invert = T),]
    }

    write.csv(all.data, file = paste0(output.name, "_SLAC_by-branch.csv"), row.names = F, quote = F)
    if (!quiet){ print(paste0("Finished ", hyphy.analysis, " data summary!")) }
    return(invisible(all.data))
  }

  #################################################################################
  ################################# BUSTED ########################################
  #################################################################################
  if (hyphy.analysis == "BUSTED"){

    all.data = data.frame()
    for (i in 1:length(results.files)){

      json.file = paste0(results.directory, "/", results.files[i], "/BUSTED-results.json")
      if (!file.exists(json.file)){ next }

      json.data = jsonlite::fromJSON(json.file)
      branch.att = unlist(json.data$`branch attributes`[[1]])

      con.stats   = branch.att[grep("\\.constrained", names(branch.att))]
      uncon.stats = branch.att[grep("\\.unconstrained", names(branch.att))]
      nucl.stats  = branch.att[grep("\\.Nucleotide", names(branch.att))]
      omega.stats = branch.att[grep("\\.MG94xREV", names(branch.att))]

      con.data   = data.frame(Sample = gsub("\\.constrained", "", names(con.stats)),
                              null_constrained = as.numeric(con.stats))
      uncon.data = data.frame(Sample = gsub("\\.unconstrained", "", names(uncon.stats)),
                              unconstrained = as.numeric(uncon.stats))
      nucl.data  = data.frame(Sample = gsub("\\.Nucleotide GTR", "", names(nucl.stats)),
                              gtr_model = as.numeric(nucl.stats))
      omega.data = data.frame(Sample = gsub("\\.MG94xREV.*", "", names(omega.stats)),
                              mg94xrev = as.numeric(omega.stats))

      if (nrow(con.data) != 0){
        a.data = merge(con.data, uncon.data, by = "Sample")
      } else {
        a.data = uncon.data
        a.data$null_constrained = NA
      }
      b.data    = merge(nucl.data, omega.data, by = "Sample")
      w.data    = merge(a.data, b.data, by = "Sample")
      save.data = cbind(Locus = results.files[i], w.data)
      all.data  = rbind(all.data, save.data)
    }

    if (tips.only == TRUE){
      all.data = all.data[grep("Node.*", all.data$Sample, invert = T),]
    }

    write.csv(all.data, file = paste0(output.name, "_BUSTED-stats.csv"), row.names = F, quote = F)
    if (!quiet){ print(paste0("Finished ", hyphy.analysis, " data summary!")) }
    return(invisible(all.data))
  }

  #################################################################################
  ################################# aBSREL ########################################
  #################################################################################
  if (hyphy.analysis == "absrel"){

    all.data = data.frame()
    for (i in 1:length(results.files)){

      json.file = paste0(results.directory, "/", results.files[i], "/aBSREL-results.json")
      if (!file.exists(json.file)){ next }

      json.data  = jsonlite::fromJSON(json.file)
      branch.att = unlist(json.data$`branch attributes`[[1]])

      full.stats  = branch.att[grep("\\.Full adaptive model$", names(branch.att))]
      fs.stats    = branch.att[grep("\\.Full adaptive model \\(synonymous subs/site\\)", names(branch.att))]
      fns.stats   = branch.att[grep("\\.Full adaptive model \\(non-synonymous subs/site\\)", names(branch.att))]
      mg94.stats  = branch.att[grep("\\.Baseline MG94xREV$", names(branch.att))]
      omega.stats = branch.att[grep("\\.Baseline MG94xREV omega ratio", names(branch.att))]
      upval.stats = branch.att[grep("\\.Uncorrected P-value", names(branch.att))]
      cpval.stats = branch.att[grep("\\.Corrected P-value", names(branch.att))]
      rate1.stats = branch.att[grep("\\.Rate Distributions1", names(branch.att))]
      rate2.stats = branch.att[grep("\\.Rate Distributions2", names(branch.att))]
      ngtr.stats  = branch.att[grep("\\.Nucleotide GTR", names(branch.att))]
      rcl.stats   = branch.att[grep("\\.Rate classes", names(branch.att))]
      lrt.stats   = branch.att[grep("\\.LRT", names(branch.att))]

      full.data  = data.frame(Sample = gsub("\\.Full adaptive model$", "", names(full.stats)), full_adaptive_model = as.numeric(full.stats))
      fs.data    = data.frame(Sample = gsub("\\.Full adaptive model \\(synonymous subs/site\\)", "", names(fs.stats)), full_adaptive_model_syn = as.numeric(fs.stats))
      fns.data   = data.frame(Sample = gsub("\\.Full adaptive model \\(non-synonymous subs/site\\)", "", names(fns.stats)), full_adaptive_model_nonsyn = as.numeric(fns.stats))
      mg94.data  = data.frame(Sample = gsub("\\.Baseline MG94xREV$", "", names(mg94.stats)), baseline_mg94xrev = as.numeric(mg94.stats))
      omega.data = data.frame(Sample = gsub("\\.Baseline MG94xREV omega ratio", "", names(omega.stats)), baseline_omega_ratio = as.numeric(omega.stats))
      upval.data = data.frame(Sample = gsub("\\.Uncorrected P-value", "", names(upval.stats)), uncorrected_pval = as.numeric(upval.stats))
      cpval.data = data.frame(Sample = gsub("\\.Corrected P-value", "", names(cpval.stats)), corrected_pval = as.numeric(cpval.stats))
      rate1.data = data.frame(Sample = gsub("\\.Rate Distributions1", "", names(rate1.stats)), rate_distributions1 = as.numeric(rate1.stats))
      rate2.data = data.frame(Sample = gsub("\\.Rate Distributions2", "", names(rate2.stats)), rate_distributions2 = as.numeric(rate2.stats))
      ngtr.data  = data.frame(Sample = gsub("\\.Nucleotide GTR", "", names(ngtr.stats)), nucleotide_gtr = as.numeric(ngtr.stats))
      rcl.data   = data.frame(Sample = gsub("\\.Rate classes", "", names(rcl.stats)), rate_classes = as.numeric(rcl.stats))
      lrt.data   = data.frame(Sample = gsub("\\.LRT", "", names(lrt.stats)), lrt = as.numeric(lrt.stats))

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
      all.data  = rbind(all.data, save.data)
    }

    if (tips.only == TRUE){
      all.data = all.data[grep("Node.*", all.data$Sample, invert = T),]
    }

    write.csv(all.data, file = paste0(output.name, "_aBSREL-stats.csv"), row.names = F, quote = F)
    if (!quiet){ print(paste0("Finished ", hyphy.analysis, " data summary!")) }
    return(invisible(all.data))
  }

}#end function
