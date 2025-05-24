#' CircosVisualization
#'
#' @param completeness A completeness of genome. If NULL, all nucleotides are recognized as chromosome driven.
#' @param fna_count A completeness of genome. If NULL, all nucleotides are recognized as chromosome driven.
#' @param gb_dir A completeness of genome. If NULL, all nucleotides are recognized as chromosome driven.
#' @return A genome map.
#' @export
#'
# Function to generate genome maps for multiple chromosomes

ContigsVisualization <- function(fna_count = NULL, completeness = FALSE, gb_dir = NULL){
  # Check if contig_length exists, if not, generate it from GenBank files
  if (!exists("contig_length", envir = .GlobalEnv)) {
    message("Generating 'contig_length' from GenBank files...")

    # If gb_dir is NULL, use the current working directory
    if (is.null(gb_dir)) {
      gb_dir <- getwd()
    }
    # Get the list of .gbk, .gb, or .gbff files in the directory
    gb_files <- list.files(
      gb_dir,
      pattern = "\\.gbk$|\\.gb$|\\.gbff$",
      full.names = TRUE
    )
    gb_files <- gb_files[!grepl("^~\\$", gb_files)]  # Exclude temporary files

    if (length(gb_files) == 0) {
      stop("No GenBank files found in the specified directory.")
    }

    # Read the first GenBank file (assuming one file for simplicity)
    genbank_file <- gb_files[1]
    origin <- read.csv(genbank_file, header = FALSE, sep = "\n")

    # Calculate contig length
    contig_length <- as.data.frame(grep("source + 1..", origin$V1, ignore.case = TRUE, value = TRUE))
    colnames(contig_length) <- "contig_length"
    contig_length <- tidyr::separate(contig_length, col = `contig_length`, into = c("start", "end"), sep = "\\.\\.")
    contig_length <- as.data.frame(contig_length$end)
    colnames(contig_length) <- "length"
    contig_length$order <- seq_len(nrow(contig_length))
    contig_length <- contig_length[order(contig_length$order), ]
    rownames(contig_length) <- NULL
    contig_length <- as.data.frame(contig_length$length)
    colnames(contig_length) <- "length"

    # Assign contig_length to the global environment
    assign("contig_length", contig_length, envir = .GlobalEnv)
    message("'contig_length' has been generated and assigned to the global environment.")
  } else {
    contig_length <- get("contig_length", envir = .GlobalEnv)
  }

  # Set fna_count to the number of rows in contig_length if it is NULL
  if (is.null(fna_count)) {
    fna_count <- nrow(contig_length)
  }

  circos.clear()
  col_text <- "grey40"
  circos.par("start.degree" = 90, "track.height" = 0.8, gap.degree = 5, cell.padding = c(0, 0, 0, 0), track.margin = c(0.02, 0.02))

  # Check if gff exists, if not, create it using GenomeDrawer
  if (!exists("gff", envir = .GlobalEnv)) {
    message("Generating 'gff' using GenomeDrawer...")

    # If gb_dir is NULL, use the current working directory
    if (is.null(gb_dir)) {
      gb_dir <- getwd()
    }

    # Get the list of .gff files in the directory
    gff_files <- list.files(
      gb_dir,
      pattern = "\\.gff$",
      full.names = TRUE
    )
    gff_files <- gff_files[!grepl("^~\\$", gff_files)]  # Exclude temporary files

    if (length(gff_files) == 0) {
      stop("No GFF files found in the specified directory.")
    }

    # Read the first GFF file (assuming one file for simplicity)
    gff_file <- gff_files[1]
    gff <- ape::read.gff(file = gff_file)
    gff <- gff %>%
      dplyr::mutate(locus_tag = substr(.$attributes, stringr::str_locate(.$attributes, "(locus_tag)\\=.{1,15}?_(RS)?[0-9]{1,}")[, 1],
                                       stringr::str_locate(.$attributes, "(locus_tag)\\=.{1,15}?_(RS)?[0-9]{1,}")[, 2])) %>%
      dplyr::mutate(locus_tag = gsub("locus_tag\\=", "", .$locus_tag))
    # Create separate GFF data for each region
    for (i in 1:length(which(gff$source == "Local" | gff$type == "region"))) {
      assign(paste((gff %>% dplyr::filter(source == "Local" | type == "region") %>% dplyr::arrange(dplyr::desc(end)) %>% dplyr::select(seqid))[i, ] %>% as.character(), "gff", sep = "_"),
             gff %>% dplyr::group_by(seqid) %>% split(.$seqid) %>% magrittr::extract2(i))
    }
  }


  # Check if genbank_table exists in the global environment
  if (!exists("genbank_table", envir = .GlobalEnv)) {
    stop("Error: 'genbank_table' does not exist in the global environment.")
  }

  # Load genbank_table from the global environment
  genbank_table <- get("genbank_table", envir = .GlobalEnv)

  # Process GenBank data
  for (i in 1:length(which(gff$source == "Local" | gff$type == "region"))) {
    current_seqid <- (gff %>% dplyr::filter(source == "Local" | type == "region") %>% dplyr::select(seqid))[i, ] %>% as.character()
    assign(paste(current_seqid, "contig_final", sep = "_"),
           genbank_table %>%
             dplyr::mutate(
               seqid = dplyr::case_when(
                 locus_tag %in% (
                   get(paste(current_seqid, "gff", sep = "_")) %>%
                     dplyr::select(locus_tag) %>%
                     na.omit() %>%
                     .$locus_tag %>%
                     c()
                 ) ~ current_seqid,
                 TRUE ~ NA_character_
               )
             ) %>%
             dplyr::filter(seqid == current_seqid))
  }

  factors <- (gff %>% filter(source == "Local" | type == "region") %>%
                dplyr::arrange(dplyr::desc(end)) %>%
                dplyr::select(seqid)) %>%
    distinct() %>%
    pull() %>%
    as.character()

  # Create xlim matrix for all factors
  xlim <- matrix(nrow = length(factors), ncol = 2)
  for (i in 1:length(factors)) {
    seqid <- factors[i]
    xlim[i, ] <- c(as.numeric(get(paste(seqid, "length", sep = "_"))$start),
                   as.numeric(get(paste(seqid, "length", sep = "_"))$end))
  }

  pdf(file = paste0(paste(Genome_summary$File,"map"), ".pdf"))

  circos.initialize(factors = factors, xlim = xlim)

  # Plot genome base circle
  circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    chr <- CELL_META$sector.index
    xlim <- CELL_META$xlim
    ylim <- CELL_META$ylim
    circos.text(mean(xlim), mean(ylim), chr, cex = 0.5, col = col_text, facing = "bending.inside", niceFacing = TRUE)
  }, bg.col = "grey50", bg.border = FALSE, track.height = 0.03)

  # Create dynamic break points for axis labeling for each factor
  for (i in 1:length(factors)) {
    genome_length <- as.numeric(get(paste(factors[i], "length", sep = "_"))$end)
    step_size <- genome_length / 24  # Divide the genome into approximately 24 labels
    step_size <- ceiling(step_size / 1000) * 1000  # Round step_size to the nearest 1,000 for readability

    # Determine label unit (Mb for large genomes, kb for smaller genomes)
    if (genome_length >= 1e6) {
      label_unit <- "Mb"
      scale_factor <- 1e6
    } else {
      label_unit <- "kb"
      scale_factor <- 1e3
    }

    brk <- seq(0, genome_length, by = step_size)

    # Add genome x-axis with dynamic labels
    circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
      circos.axis(h = "top",
                  major.at = brk,
                  labels = paste(round(brk / scale_factor, 2), label_unit),
                  labels.cex = 0.8,
                  col = col_text,
                  labels.col = col_text,
                  lwd = 0.7,
                  labels.facing = "clockwise")
    }, bg.border = FALSE)
  }

  # Plot COG for 3' to 5' and 5' to 3' separately
  circos.track(ylim = c(0, 1), bg.col = "grey90", bg.border = FALSE, track.height = 0.06, panel.fun = function(x, y) {
    sector_index <- CELL_META$sector.index
    if (exists(paste(sector_index, "3to5", sep = "_"), inherits = TRUE)) {
      cog_data_3to5 <- get(paste(sector_index, "3to5", sep = "_"))
      if (!is.null(cog_data_3to5) && nrow(cog_data_3to5) > 0) {
        circos.genomicRect(
          region = matrix(c(as.numeric(cog_data_3to5$start), as.numeric(cog_data_3to5$end)), ncol = 2),
          value = cog_data_3to5$COG_category_for_plot,
          col = cog_data_3to5$COG_color,
          border = NA
        )
      }
    }
  })

  circos.track(ylim = c(0, 1), bg.col = "grey90", bg.border = FALSE, track.height = 0.06, panel.fun = function(x, y) {
    sector_index <- CELL_META$sector.index
    if (exists(paste(sector_index, "5to3", sep = "_"), inherits = TRUE)) {
      cog_data_5to3 <- get(paste(sector_index, "5to3", sep = "_"))
      if (!is.null(cog_data_5to3) && nrow(cog_data_5to3) > 0) {
        circos.genomicRect(
          region = matrix(c(as.numeric(cog_data_5to3$start), as.numeric(cog_data_5to3$end)), ncol = 2),
          value = cog_data_5to3$COG_category_for_plot,
          col = cog_data_5to3$COG_color,
          border = NA
        )
      }
    }
  })

  # Plot RNA for all factors in the same track
  circos.track(ylim = c(0, 1), bg.col = "grey90", bg.border = FALSE, track.height = 0.06, panel.fun = function(x, y) {
    sector_index <- CELL_META$sector.index
    if (!is.null(RNA) && any(RNA$seqid == sector_index)) {
      RNA_data <- RNA %>% filter(seqid == sector_index)
      if (nrow(RNA_data) > 0) {
        circos.genomicRect(
          region = matrix(c(as.numeric(RNA_data$start), as.numeric(RNA_data$end)), ncol = 2),
          value = RNA_data$type,
          col = RNA_data$color,
          border = NA
        )
      }
    }
  })

  # Function to calculate global ylim for a given type across all contigs
  calculate_global_ylim <- function(factors, type, scale_factor = 1.2) {  # Increased scale_factor
    max_range <- 0
    for (i in seq_along(factors)) {
      mapped_index <- paste0("window_contig_", i, "_seq")
      if (exists(mapped_index, inherits = TRUE)) {
        gc_data <- get(mapped_index, inherits = TRUE)
        if (!is.null(gc_data) && nrow(gc_data) > 0) {
          valid_values <- as.numeric(gc_data[[paste0(type, "_minus_average")]])
          if (any(!is.na(valid_values))) {
            max_range <- max(max_range, max(abs(valid_values), na.rm = TRUE))
          }
        }
      }
    }
    if (max_range == 0) {
      return(c(-0.1, 0.1))  # Reduced default range for tighter control
    }
    return(c(-max_range * scale_factor, max_range * scale_factor))
  }

  # Replace the GC skew and GC ratio plotting loop with this
  for (type in c("gc_skew", "gc_ratio")) {
    # Compute global ylim for this type
    ylim <- calculate_global_ylim(factors, type, scale_factor = 1.2)

    circos.track(bg.col = NA, bg.border = FALSE, track.height = 0.08, panel.fun = function(x, y) {  # Reduced track.height
      sector_index <- CELL_META$sector.index
      contig_index <- which(factors == sector_index)

      if (length(contig_index) > 0) {
        mapped_index <- paste0("window_contig_", contig_index, "_seq")

        if (exists(mapped_index, inherits = TRUE)) {
          gc_data <- get(mapped_index, inherits = TRUE)

          if (!is.null(gc_data) && nrow(gc_data) > 0 && any(!is.na(gc_data[[paste0(type, "_minus_average")]]))) {
            circos.barplot(
              value = as.numeric(gc_data[[paste0(type, "_minus_average")]]),
              pos = (gc_data$start + gc_data$end) / 2,
              col = gc_data[[paste0(type, "_color")]],
              border = NA,
              bar_width = 10000
            )
          } else {
            message(paste("No valid data for", type, "in", sector_index))
          }
        } else {
          message(paste("Data for", mapped_index, "not found."))
        }
      } else {
        message(paste("Sector index", sector_index, "not found in factors."))
      }
    }, ylim = ylim)
  }


  # Add genome name and total base pairs at the center
  text(0, 0, paste(Genome_summary$File,
                   paste(formatC(sum(as.numeric(contig_length$length)), format = "f", digits = 0, big.mark = ","), "bp"),
                   sep = "\n"), cex = 1.5)

  # Draw legend for COG
  COG <- genbank_table %>% dplyr::group_by(COG_legend) %>% dplyr::slice_head() %>% dplyr::select(COG_category_for_plot, COG_color, COG_legend) %>% na.omit()
  lgd_list <- ComplexHeatmap::Legend(
    labels = COG$COG_legend,
    title = "COG category",
    type = "points",
    legend_gp = grid::gpar(col = NA),
    ncol = 1,
    by_row = TRUE,
    direction = "vertical",
    background = COG$COG_color
  )
  pd <- ComplexHeatmap::packLegend(lgd_list)

  ComplexHeatmap::draw(pd, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))
  dev.off()
  circos.clear()
}
