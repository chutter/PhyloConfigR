#' @title alignmentsToStructure
#'
#' @description Function to trim a set of alignments to a provided target. The function operates across a directory of alignments that correspond to a single fasta file of capture targets
#'
#' @param alignment.directory path to a folder of sequence alignments
#'
#' @param alignment.format available input alignment formats: fasta or phylip
#'
#' @param output.directory new alignment directory where the trimmed output files are saved
#'
#' @param output.format available output formats: phylip
#'
#' @param target.file path to the fasta file with the target sequences. These should be the entire marker, not the probe.
#'
#' @param target.direction TRUE ensures output alignments are the same direction as the targets
#'
#' @param min.alignment.length minimum alignment length to save in bp (default: 100)
#'
#' @param min.taxa.alignment mininum number of taxa to save alignment (default: 4)
#'
#' @param threads number of CPU threads / processes
#'
#' @param memory memory in GB
#'
#' @param overwrite TRUE to overwrite output files with the same name
#'
#' @param mafft.path system path to the mafft program
#'
#' @return a new output directory with the trimmed alignments
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

alignmentsToStructure = function(alignment.directory = NULL,
                                alignment.format = "phylip",
                                population.file = NULL,
                                output.file = NULL,
                                min.count.taxa = 4,
                                min.percent.taxa = 50,
                                min.percent.column = 50,
                                threads = 1,
                                memory = 1,
                                overwrite = FALSE,
                                mafft.path = NULL) {


  # alignment.directory = "/Volumes/LaCie/Structure_maker/trimmed_genes"
  # alignment.format = "phylip"
  # output.file = "/Volumes/LaCie/Structure_maker/output_genes_filtered.str"
  # population.file = "/Volumes/LaCie/Structure_maker/Populations.csv"
  # threads = 1 #no multithreadhing yet
  # memory = 10
  # overwrite = TRUE

  # # Filters
  # min.count.taxa = 4
  # min.percent.taxa = 75 #whichever is more strict w/ min.count.taxa
  # min.percent.column = 75



  require(foreach)

  # if (dir.exists(output.dir) == FALSE) { dir.create(output.dir) }

  if (file.exists(output.file) == TRUE) {
    if (overwrite == TRUE) {
      system(paste0("rm ", output.file))
    }
  }

  #Gathers alignments
  align.files = list.files(alignment.directory)

  if (length(align.files) == 0) { stop("alignment files could not be found.") }

  # Sets up the empty file and doubles it up
  pop.file = read.csv(population.file)
  pop.file2 = pop.file
  pop.file$Sample = paste0(pop.file$Sample, "_1")
  pop.file2$Sample = paste0(pop.file2$Sample, "_2")

  final.locus = rbind(pop.file, pop.file2)
  final.locus = final.locus[order(final.locus$Sample), ]

  # Loops through each alignment
  for (i in seq_along(align.files)) {

    #Reads in phylip format
    if (alignment.format == "phylip") {
      align <- Biostrings::readAAMultipleAlignment(file = paste0(alignment.directory, "/", align.files[i]), format = "phylip")

      #  align = Biostrings::readDNAStringSet(file = paste0(alignment.dir, "/", align.files[i]), format = "phylip")
      #  align = readLines(paste0(alignment.dir, "/", align.files[i]))[-1]
      #  align = gsub(".*\\ ", "", align)
      #  char.count = nchar(align)

      align <- Biostrings::DNAStringSet(align)
      save.name <- gsub(".phy$", "", align.files[i])
      save.name <- gsub(".phylip$", "", save.name)
    } # end phylip

    #Reads in fasta format
    if (alignment.format == "fasta") {
      align <- Biostrings::readDNAStringSet(paste0(alignment.directory, "/", align.files[i]))
      save.name <- gsub(".fa$", "", align.files[i])
      save.name <- gsub(".fasta$", "", save.name)
    } # end phylip

    if (is.null(min.count.taxa) != TRUE) {
      if (length(align) <= min.count.taxa) {
        print(paste0(align.files[i], "alignment too few taxa, skip."))
        next
      }
    }

    if (is.null(min.percent.taxa) != TRUE) {
      if (length(align) / nrow(pop.file) <= min.percent.taxa / 100) {
        print(paste0(align.files[i], "alignment too few taxa, skip."))
        next
      }
    }

    # alignment length finds
    alignment = as.matrix(align)
    x.len = dim(alignment)[2]

    # empty data frame to fill later
    save.data = data.frame(
      Column = as.numeric(),
      Het_Count = as.numeric(),
      Miss_Count = as.numeric(),
      Per_Miss = as.numeric()
    )
    #Loops through each column to tablulate heterozygous columns and those with missing data
    for (k in 1:x.len){
      # Pulls out column
      align.col = alignment[, k]
      # Gets character count for column
      col.counts = table(align.col)[!names(table(align.col)) %in% c("N", "-")]
      #Counts het sizes and missing data
      het.count = stringr::str_count(as.character(paste0(align.col, collapse = "")), pattern = c("R", "Y", "K", "M", "S", "W", "B", "D", "H", "V"))
      miss.count = stringr::str_count(as.character(paste0(align.col, collapse = "")), pattern = c("N", "-", "\\?"))
      #Saves data
      temp.save = data.frame(
        Column = k,
        Het_Count = sum(het.count),
        Miss_Count = sum(miss.count),
        Per_Miss = (length(align.col) - sum(miss.count)) / nrow(pop.file)
      )
      save.data = rbind(save.data, temp.save)

    } # end k loop

    #removes columns with too much missing data
    if (is.null(min.percent.column) != TRUE) {
      save.data = save.data[save.data$Per_Miss >= min.percent.column / 100, ]
      if (nrow(save.data) < 10) {
        print(paste0(align.files[i], "alignment too few characters, skip."))
        next
      }
    }

    #Selects best column by maximizing heterozygosity and minmimizing missing data
    best.column = save.data[save.data$Het_Count == max(save.data$Het_Count), ]
    best.column = best.column[best.column$Miss_Count == min(best.column$Miss_Count), ]

    #Randomly selects columns if there ties
    if (nrow(best.column) != 1) {
      best.column = best.column[sample(1:nrow(best.column), 1), ]
    }

    #to do, better way to deal with ties
    #if (nrow(best.column) > 1){ stop("too many") }

    #Seelcts data from best column
    final.site = alignment[, best.column$Column]

    #Formats final data from alignment using best column data
    final.align = data.frame(names(final.site), final.site)
    rownames(final.align) = NULL
    colnames(final.align) = c("Sample", gsub(".phy$", "", align.files[i]))

    #Makes sure each row is unique or it all goes to hell
    final.align2 <- final.align
    final.align$Sample <- paste0(final.align$Sample, "_1")
    final.align2$Sample <- paste0(final.align2$Sample, "_2")

    final.align = rbind(final.align, final.align2)
    final.align = final.align[order(final.align$Sample),]

    for (k in 1:nrow(final.align)){

      if (final.align[k, 2] == "M") {
        final.align[k, 2] <- "A"
        final.align[k + 1, 2] <- "C"
        k <- k + 1
      }

      if (final.align[k,2] == "R"){
          final.align[k, 2] = "A"
          final.align[k + 1, 2] = "G"
          k = k + 1
      }

      if (final.align[k, 2] == "W") {
        final.align[k, 2] <- "A"
        final.align[k + 1, 2] <- "T"
        k <- k + 1
      }

      if (final.align[k, 2] == "S") {
        final.align[k, 2] <- "C"
        final.align[k + 1, 2] <- "G"
        k <- k + 1
      }

      if (final.align[k, 2] == "Y") {
        final.align[k, 2] <- "C"
        final.align[k + 1, 2] <- "T"
        k <- k + 1
      }

      if (final.align[k, 2] == "K") {
        final.align[k, 2] <- "G"
        final.align[k + 1, 2] <- "T"
        k <- k + 1
      }

    }#end k loop

    #Assigns numerical values
    final.align[final.align == "A"] <- 0
    final.align[final.align == "T"] <- 1
    final.align[final.align == "G"] <- 2
    final.align[final.align == "C"] <- 3
    final.align[final.align == "N"] <- -9
    final.align[final.align == "-"] <- -9
    final.align[final.align == "?"] <- -9

    # Data frame structure
    # Label  Locus1 Locus2
    # Samp1    0     0
    # Samp1    1     0
    # Samp2    0     1
    # Samp2    0     1

    # Adds to combined data
    final.locus = merge(final.locus, final.align, all = TRUE, by = "Sample")
    final.locus[is.na(final.locus) == TRUE] <- -9

    print(paste0("alignment number ", i, " finished!"))

  } # end i loop

  final.locus$Sample = gsub("_1$|_2$", "", final.locus$Sample)
  colnames(final.locus)[3:length(colnames(final.locus))] = seq(1:(length(colnames(final.locus)) - 2))
  #final.locus$Label = gsub(".*_", "", final.locus$Label)

  write.table(final.locus, output.file,
    row.names = FALSE, sep = " ", quote = FALSE, col.names = FALSE
  )

  #########################
  ###### END SCRIPT
  #########################

} #end function
