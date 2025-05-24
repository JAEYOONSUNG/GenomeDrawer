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
    fna_count <- nrow(contig_length) # Assuming contig_length is defined elsewhere
  }

  for (i in 1:fna_count) {
    if (completeness) {
      # Assuming gff is defined elsewhere
      grDevices::pdf(file = paste0(paste(paste((gff %>% dplyr::filter(source == "Local" | type == "region") %>% dplyr::select(seqid))[i, ] %>% as.character()), "map", sep="_"), ".pdf"))
    }

    # Clear previous plots and set up circos parameters
    circlize::circos.clear()
    col_text <- "grey40"
    circlize::circos.par("start.degree" = 90, "track.height" = 0.8, gap.degree = 0, cell.padding = c(0, 0, 0, 0))

    # Initialize factors for the chromosome i
    # Assuming gff is defined elsewhere
    seqid <- paste((gff %>% dplyr::filter(source == "Local"|type == "region") %>% dplyr::arrange(dplyr::desc(end)) %>% dplyr::select(seqid))[i, ] %>% as.character())
    chr_length <- get(paste(seqid, "length", sep = "_")) # Relies on object 'seqid_length' existing

    circlize::circos.initialize(factors = c(seqid), xlim = matrix(c(chr_length$start, chr_length$end), ncol = 2))

    # Plot genome base circle
    circlize::circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
      chr <- circlize::get.cell.meta.data("sector.index") # Corrected way to get CELL_META components
      xlim <- circlize::get.cell.meta.data("xlim")
      ylim <- circlize::get.cell.meta.data("ylim")
      circlize::circos.text(mean(xlim), mean(ylim), chr, cex = 0.5, col = col_text, facing = "bending.inside", niceFacing = TRUE)
    }, bg.col = "grey50", bg.border = FALSE, track.height = 0.03)

    # Calculate genome length and determine step size and label units
    # Assuming nucleotide is defined elsewhere
    genome_length <- nchar(get(paste(names(nucleotide)[i], sep = "_"))) # Relies on object 'nucleotideName_' existing
    step_size <- genome_length / 24 # Divide the genome into approximately 24 labels
    step_size <- ceiling(step_size / 1000) * 1000 # Round step_size to the nearest 1,000 for readability

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
    circlize::circos.track(track.index = circlize::get.current.track.index(), panel.fun = function(x, y) {
      circlize::circos.axis(h = "top",
                            major.at = brk, # Use dynamically generated break points
                            labels = paste(round(brk / scale_factor, 2), label_unit), # Labels in kb or Mb
                            labels.cex = 0.8,
                            col = col_text,
                            labels.col = col_text,
                            lwd = 0.7,
                            labels.facing = "clockwise")
    }, bg.border = FALSE)

    # Plot COG in 3' to 5' and 5' to 3' directions if seqid exists and data is available
    plot_cog <- function(seqid_suffix) {
      cog_data <- get(paste(seqid, seqid_suffix, sep = "_")) # Relies on 'seqid_3to5' or 'seqid_5to3' existing
      if (!is.null(cog_data) && nrow(cog_data) > 0) {
        circlize::circos.track(factors = cog_data$seqid,
                               panel.fun = function(region, value, ...) { # 'region' and 'value' are arguments for panel.fun, not direct data
                                 circlize::circos.genomicRect(
                                   region = matrix(c(cog_data$start %>% as.numeric(),
                                                     cog_data$end %>% as.numeric()), ncol = 2),
                                   value = cog_data$COG_category_for_plot, # This should be a vector of values for the regions
                                   col = cog_data$COG_color, # This should be a vector of colors
                                   border = NA
                                 )
                               }, ylim = c(0, 1), track.index = circlize::get.current.track.index() + 1, bg.col = "grey90", bg.border = FALSE, track.height = 0.06)
      } else {
        message(paste("Skipping COG track for", seqid, "suffix", seqid_suffix, "due to lack of data."))
      }
    }

    plot_cog("3to5")
    plot_cog("5to3")

    # Plot RNA annotations if seqid exists and data is available
    # Assuming RNA is defined elsewhere
    if (!is.null(RNA) && any(RNA$seqid == seqid)) {
      RNA_data <- RNA %>% dplyr::filter(seqid == seqid)
      if (nrow(RNA_data) > 0) {
        circlize::circos.track(factors = RNA_data$seqid %>% as.character(),
                               panel.fun = function(region, value, ...) { # 'region' and 'value'
                                 circlize::circos.genomicRect(
                                   region = matrix(c(RNA_data$start %>% as.numeric(), RNA_data$end %>% as.numeric()), ncol = 2),
                                   value = RNA_data$type %>% as.character(), # This should be a vector of values
                                   col = RNA_data$color, # This should be a vector of colors
                                   border = NA
                                 )
                               }, ylim = c(0, 1), track.index = circlize::get.current.track.index() + 1, bg.col = "grey90", bg.border = FALSE, track.height = 0.06)
      } else {
        message(paste("Skipping RNA track for", seqid, "due to lack of data."))
      }
    }

    # Plot GC skew and GC ratio tracks if data exists
    plot_gc_track <- function(type) {
      # Relies on 'window_nucleotideName_' existing
      gc_data <- get(paste("window", names(nucleotide)[i], sep = "_"), inherits = TRUE)
      if (!is.null(gc_data) && nrow(gc_data) > 0 && all(gc_data[[paste0(type, "_minus_average")]] != 0)) {
        if (any(!is.na(gc_data[[paste0(type, "_minus_average")]]))) {
          circlize::circos.track(factors = gc_data$seqid, # factors should match the sector
                                 panel.fun = function(x,y) { # x and y are implicit here for the current sector
                                   current_sector_data <- gc_data[gc_data$seqid == circlize::get.cell.meta.data("sector.index"),]
                                   if(nrow(current_sector_data) > 0){
                                     circlize::circos.barplot(
                                       value = current_sector_data[[paste0(type, "_minus_average")]] %>% as.numeric(),
                                       pos = (current_sector_data$start + current_sector_data$end) / 2,
                                       col = current_sector_data[[paste0(type, "_color")]],
                                       border = NA,
                                       bar_width = 10000 # bar_width might need adjustment based on data scale
                                     )
                                   }
                                 }, ylim = c(-max(abs(gc_data[[paste0(type, "_minus_average")]] %>% as.numeric()), na.rm = TRUE),
                                             max(abs(gc_data[[paste0(type, "_minus_average")]] %>% as.numeric()), na.rm = TRUE)),
                                 track.index = circlize::get.current.track.index() + 1, bg.col = NA, bg.border = FALSE, track.height = 0.12)
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
    graphics::text(0, 0, paste(seqid,
                               paste(formatC(nchar(get(paste(names(nucleotide)[i], sep = "_"))),format = "f", digits = 0, big.mark = ","),"bp"),
                               sep = "\n"), cex = 1.5)

    # Draw legend for COG
    # Assuming genbank_table is defined elsewhere
    COG <- genbank_table %>% dplyr::group_by(COG_legend) %>% dplyr::slice_head() %>% dplyr::select(COG_category_for_plot, COG_color, COG_legend) %>% stats::na.omit()
    lgd_list = ComplexHeatmap::Legend(
      labels = COG$COG_legend,
      title = "COG category",
      type = "points",
      legend_gp = grid::gpar(col = NA), # Using grid::gpar
      ncol = 1,
      by_row = TRUE,
      direction = "vertical",
      background = COG$COG_color
    )
    pd = ComplexHeatmap::packLegend(lgd_list)

    ComplexHeatmap::draw(pd, x = grid::unit(4, "mm"), y = grid::unit(4, "mm"), just = c("left", "bottom")) # Using grid::unit
    grDevices::dev.off()
  }
  circlize::circos.clear()
}
