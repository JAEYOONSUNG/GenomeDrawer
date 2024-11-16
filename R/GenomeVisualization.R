#' CircosVisualization
#'
#' @param i A number of contigs.
#' @param completeness A completeness of genome. If NULL, all nucleotides are recognized as chromosome driven.
#' @return A genome map.
#' @export
#'

# Function to generate genome maps for multiple chromosomes
GenomeVisualization <- function(fna_count = NULL, completeness = TRUE) {
  # Set fna_count to the number of rows in contig_length if it is NULL
  if (is.null(fna_count)) {
    fna_count <- nrow(contig_length)
  }

  for (i in 1:fna_count) {
    if (completeness) {
      pdf(file = paste0(paste(paste((gff %>% filter(source == "Local" | type == "region") %>% dplyr::select(seqid))[i, ] %>% as.character()), "map", sep="_"), ".pdf"))
    }

    # Clear previous plots and set up circos parameters
    circos.clear()
    col_text <- "grey40"
    circos.par("start.degree" = 90, "track.height" = 0.8, gap.degree = 0, cell.padding = c(0, 0, 0, 0))

    # Initialize factors for the chromosome i
    seqid <- paste((gff %>% filter(source == "Local"|type == "region") %>% dplyr::arrange(dplyr::desc(end)) %>% dplyr::select(seqid))[i, ] %>% as.character())
    chr_length <- get(paste(seqid, "length", sep = "_"))

    circos.initialize(factors = c(seqid), xlim = matrix(c(chr_length$start, chr_length$end), ncol = 2))

    # Plot genome base circle
    circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
      chr <- CELL_META$sector.index
      xlim <- CELL_META$xlim
      ylim <- CELL_META$ylim
      circos.text(mean(xlim), mean(ylim), chr, cex = 0.5, col = col_text, facing = "bending.inside", niceFacing = TRUE)
    }, bg.col = "grey50", bg.border = FALSE, track.height = 0.03)

    # Calculate genome length and determine step size and label units
    genome_length <- nchar(get(paste(names(nucleotide)[i], sep = "_")))
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

    # Create dynamic break points for axis labeling
    brk <- seq(0, genome_length, by = step_size)

    # Add genome x-axis with dynamic labels
    circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
      circos.axis(h = "top",
                  major.at = brk,  # Use dynamically generated break points
                  labels = paste(round(brk / scale_factor, 2), label_unit),  # Labels in kb or Mb
                  labels.cex = 0.8,
                  col = col_text,
                  labels.col = col_text,
                  lwd = 0.7,
                  labels.facing = "clockwise")
    }, bg.border = FALSE)

    # Plot COG in 3' to 5' and 5' to 3' directions if seqid exists and data is available
    plot_cog <- function(seqid_suffix) {
      cog_data <- get(paste(seqid, seqid_suffix, sep = "_"))
      if (!is.null(cog_data) && nrow(cog_data) > 0) {
        circos.track(factors = cog_data$seqid,
                     panel.fun = function(region, value, ...) {
                       circos.genomicRect(
                         region = matrix(c(cog_data$start %>% as.numeric(),
                                           cog_data$end %>% as.numeric()), ncol = 2),
                         value = cog_data$COG_category_for_plot,
                         col = cog_data$COG_color,
                         border = NA
                       )
                     }, ylim = c(0, 1), track.index = get.current.track.index() + 1, bg.col = "grey90", bg.border = FALSE, track.height = 0.06)
      } else {
        message(paste("Skipping COG track for", seqid, "suffix", seqid_suffix, "due to lack of data."))
      }
    }

    plot_cog("3to5")
    plot_cog("5to3")

    # Plot RNA annotations if seqid exists and data is available
    if (!is.null(RNA) && any(RNA$seqid == seqid)) {
      RNA_data <- RNA %>% filter(seqid == seqid)
      if (nrow(RNA_data) > 0) {
        circos.track(factors = RNA_data$seqid %>% as.character(),
                     panel.fun = function(region, value, ...) {
                       circos.genomicRect(
                         region = matrix(c(RNA_data$start %>% as.numeric(), RNA_data$end %>% as.numeric()), ncol = 2),
                         value = RNA_data$type %>% as.character(),
                         col = RNA_data$color,
                         border = NA
                       )
                     }, ylim = c(0, 1), track.index = get.current.track.index() + 1, bg.col = "grey90", bg.border = FALSE, track.height = 0.06)
      } else {
        message(paste("Skipping RNA track for", seqid, "due to lack of data."))
      }
    }

    # Plot GC skew and GC ratio tracks if data exists
    plot_gc_track <- function(type) {
      gc_data <- get(paste("window", names(nucleotide)[i], sep = "_"), inherits = TRUE)
      if (!is.null(gc_data) && nrow(gc_data) > 0 && all(gc_data[[paste0(type, "_minus_average")]] != 0)) {
        if (any(!is.na(gc_data[[paste0(type, "_minus_average")]]))) {
          circos.track(factors = gc_data$seqid,
                       panel.fun = function(value, pos, ...) {
                         circos.barplot(
                           value = gc_data[[paste0(type, "_minus_average")]] %>% as.numeric(),
                           pos = (gc_data$start + gc_data$end) / 2,
                           col = gc_data[[paste0(type, "_color")]],
                           border = NA,
                           bar_width = 10000
                         )
                       }, ylim = c(-max(abs(gc_data[[paste0(type, "_minus_average")]] %>% as.numeric()), na.rm = TRUE),
                                   max(abs(gc_data[[paste0(type, "_minus_average")]] %>% as.numeric()), na.rm = TRUE)),
                       track.index = get.current.track.index() + 1, bg.col = NA, bg.border = FALSE, track.height = 0.12)
        } else {
          message(paste("Skipping GC track for", seqid, "type", type, "due to lack of non-NA data."))
        }
      } else {
        message(paste("Skipping GC track for", seqid, "type", type, "due to lack of data or zero range."))
      }
    }

    # Check and plot GC skew and GC ratio within loop to prevent errors
    if (exists(paste("window", names(nucleotide)[i], sep = "_"))) {
      plot_gc_track("gc_skew")
      plot_gc_track("gc_ratio")
    } else {
      message(paste("Skipping GC skew and ratio tracks for", seqid, "due to missing data."))
    }

    # Add genome name and total base pairs at the center
    text(0, 0, paste(seqid,
                     paste(formatC(nchar(get(paste(names(nucleotide)[i], sep = "_"))),format = "f", digits = 0, big.mark = ","),"bp"),
                     sep = "\n"), cex = 1.5)

    # Draw legend for COG
    COG <- genbank_table %>% group_by(COG_legend) %>% slice_head() %>% dplyr::select(COG_category_for_plot, COG_color, COG_legend) %>% na.omit()
    lgd_list = ComplexHeatmap::Legend(
      labels = COG$COG_legend,
      title = "COG category",
      type = "points",
      legend_gp = gpar(col = NA),
      ncol = 1,
      by_row = TRUE,
      direction = "vertical",
      background = COG$COG_color
    )
    pd = ComplexHeatmap::packLegend(lgd_list)

    draw(pd, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))
    dev.off()
  }
  circos.clear()
}

