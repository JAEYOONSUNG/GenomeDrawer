#' GenomeDrawer #
#'
#' @param dir A string specifying the directory where the gff output files are located. If NULL, the current working directory is used.
#' @return A list containing information for plotting genome map.
#' @export
#'

GenomeDrawer <- function(dir = NULL) {
  # If dir is NULL, use the current working directory
  if (is.null(dir)) {
    dir <- getwd()
  }

  # Get the list of .tsv files in the directory
  gff_files <- list.files(
    dir,
    pattern = "\\.gff$",
    full.names = FALSE
  )
  # Exclude temporary files starting with '~$'
  gff_files <- gff_files[!grepl("^~\\$", gff_files)]

  # Check if genbank_table exists in the global environment
  if (!exists("genbank_table", envir = .GlobalEnv)) {
    stop("Error: 'genbank_table' does not exist in the global environment.")
  }

  # Load genbank_table from the global environment
  genbank_table <- get("genbank_table", envir = .GlobalEnv)

  # Load and process GFF files

  gff <- ape::read.gff(file = gff_files)
  gff <- gff %>%
    dplyr::mutate(locus_tag = substr(.$attributes, stringr::str_locate(.$attributes, "(locus_tag)\\=.{1,15}?_(RS)?[0-9]{1,}")[, 1],
                                     stringr::str_locate(.$attributes, "(locus_tag)\\=.{1,15}?_(RS)?[0-9]{1,}")[, 2])) %>%
    dplyr::mutate(locus_tag = gsub("locus_tag\\=", "", .$locus_tag))

  # Create separate GFF data for each region
  for (i in 1:length(which(gff$source == "Local" | gff$type == "region"))) {
    assign(paste((gff %>% dplyr::filter(source == "Local" | type == "region") %>% dplyr::arrange(dplyr::desc(end)) %>% dplyr::select(seqid))[i, ] %>% as.character(), "gff", sep = "_"),
           gff %>% dplyr::group_by(seqid) %>% split(.$seqid) %>% magrittr::extract2(i))
  }


  # Process GenBank data
  for (i in 1:length(which(gff$source == "Local" | gff$type == "region"))) {
    current_seqid <- (gff %>% dplyr::filter(source == "Local" | type == "region") %>% dplyr::select(seqid))[i, ] %>% as.character()
    assign(paste(current_seqid, "contig_final", sep = "_"),
           genbank_table %>%
             dplyr::mutate(seqid = dplyr::case_when(locus_tag %in% (get(paste(current_seqid, "gff", sep = "_")) %>% dplyr::select(locus_tag) %>%na.omit() %>% .$locus_tag %>% c()) ~ current_seqid,
                 TRUE ~ NA_character_)) %>%dplyr::filter(seqid == current_seqid), envir = .GlobalEnv)
  }

  # Extract contig lengths
  for (i in 1:length(which(gff$source == "Local" | gff$type == "region"))) {
    assign(paste((gff %>% dplyr::filter(source == "Local" | type == "region") %>% dplyr::select(seqid))[i, ] %>% as.character(), "length", sep = "_"),
           get(paste((gff %>% dplyr::filter(source == "Local" | type == "region") %>% dplyr::select(seqid))[i, ] %>% as.character(), "gff", sep = "_")) %>% dplyr::slice_head(),
           envir = .GlobalEnv)
  }

  # Remove the first entry in each GFF data
  for (i in 1:length(which(gff$source == "Local" | gff$type == "region"))) {
    assign(paste((gff %>% filter(source == "Local" | type == "region") %>% dplyr::select(seqid))[i, ] %>% as.character(), "gff", sep = "_"),
           get(paste((gff %>% filter(source == "Local" | type == "region") %>% dplyr::select(seqid))[i, ] %>% as.character(), "gff", sep = "_")) %>% dplyr::slice(-1),
           envir = .GlobalEnv)
  }

  # Extract + strand data
  for (i in 1:length(which(gff$source == "Local" | gff$type == "region"))) {
    assign(paste((gff %>% filter(source == "Local" | type == "region") %>% dplyr::select(seqid))[i, ] %>% as.character(), "5to3", sep = "_"),
           get(paste((gff %>% filter(source == "Local" | type == "region") %>% dplyr::select(seqid))[i, ] %>% as.character(), "contig_final", sep = "_")) %>%
             dplyr::filter(direction == "+") %>%
             dplyr::select(seqid, locus_tag, start, end, COG_category_for_plot, COG_color),
            envir = .GlobalEnv)
  }

  # Extract - strand data
  for (i in 1:length(which(gff$source == "Local" | gff$type == "region"))) {
    assign(paste((gff %>% filter(source == "Local" | type == "region") %>% dplyr::select(seqid))[i, ] %>% as.character(), "3to5", sep = "_"),
           get(paste((gff %>% filter(source == "Local" | type == "region") %>% dplyr::select(seqid))[i, ] %>% as.character(), "contig_final", sep = "_")) %>%
             dplyr::filter(direction == "-") %>%
             dplyr::select(seqid, locus_tag, start, end, COG_category_for_plot, COG_color),
           envir = .GlobalEnv)
  }

  # Process RNA data
  RNA <- gff %>% filter(grepl("tRNA|rRNA", type)) %>%
    dplyr::mutate(locus_tag = substr(.$attributes, stringr::str_locate(.$attributes, "(locus_tag)\\=.{5,15}?_(RS)?[0-9]{1,}")[, 1], stringr::str_locate(.$attributes, "(locus_tag)\\=.{5,15}?_(RS)?[0-9]{1,}")[, 2])) %>%
    dplyr::mutate(locus_tag = gsub("locus_tag\\=", "", .$locus_tag)) %>%
    dplyr::select(seqid, locus_tag, type, start, end)

  RNA <- RNA %>% dplyr::mutate(color = dplyr::case_when(
    type %in% c("tRNA") ~ 'blue',
    type %in% c("rRNA") ~ 'red'
  ))

  assign("RNA", RNA, envir = .GlobalEnv)

  # Initialize an empty list named nucleotide
  nucleotide <- list()
  for (i in 1:nrow(contig_length)) {
    nucleotide[[i]] <- get(paste("contig", i, "seq", sep = "_"))
  }
  names(nucleotide) <- paste("contig", 1:nrow(contig_length), "seq", sep = "_")

  for (i in 1:length(names(nucleotide))) {
  # Generate sequence ID
  seqid <- paste((gff %>% filter(source == "Local" | type == "region") %>% dplyr::select(seqid))[i, ] %>% as.character())
  seq_length <- nchar(get(paste(names(nucleotide)[i], sep = "_")))

  # Check for sequences shorter than 10,000 bp
  if (seq_length <= 10000) {
    # Single range for sequences shorter than or equal to 10,000 bp
    start_vals <- 1
    end_vals <- seq_length
  } else {
    # Generate ranges for sequences longer than 10,000 bp
    start_vals <- c(10000 * (1:ceiling(seq_length / 10000)) - 9999)
    end_vals <- c(10000 * (1:floor(seq_length / 10000)), seq_length)
  }

  # Remove duplicates and invalid ranges
  window_data <- data.table::data.table(
    seqid = seqid,
    start = start_vals,
    end = end_vals
  ) %>%
    dplyr::filter(end > start) %>%     # Remove invalid ranges where end <= start
    dplyr::distinct(start, end, .keep_all = TRUE)  # Remove duplicate ranges

  # Assign the cleaned data to the global environment
  assign(paste("window", names(nucleotide)[i], sep = "_"), window_data, envir = .GlobalEnv)
}

  for (i in 1:length(names(nucleotide))) {
    assign(paste("window", names(nucleotide)[i], sep = "_"),
           get(paste("window", names(nucleotide)[i], sep = "_")) %>% dplyr::mutate(
             "window_nt" = stringr::str_sub(get(paste(names(nucleotide)[i], sep = "_")),
                                            start = start,
                                            end = end)),
           envir = .GlobalEnv)
  }

  for (i in 1:length(names(nucleotide))) {
    assign(paste(names(nucleotide)[i], "RNA", sep = "_"),
           get(paste("window", names(nucleotide)[i], sep = "_")) %>%
             dplyr::mutate(seqid = (gff %>% dplyr::filter(source == "Local" | type == "region") %>% dplyr::select(seqid))[i, ] %>% as.character()),
           envir = .GlobalEnv)
  }

  # Initialize an empty list named nucleotide
  nucleotide <- list()
  for (i in 1:nrow(contig_length)) {
    nucleotide[[i]] <- get(paste("contig", i, "seq", sep = "_"))
  }
  names(nucleotide) <- paste("contig", 1:nrow(contig_length), "seq", sep = "_")
  assign("nucleotide", nucleotide, envir = .GlobalEnv)

  for (i in 1:length(names(nucleotide))) {
    assign(paste("window", names(nucleotide)[i], sep = "_"),
           data.table::data.table(
             seqid = paste((gff %>% filter(source == "Local" | type == "region") %>% dplyr::select(seqid))[i, ] %>% as.character()),
             start = c(10000 * (1:ceiling(nchar(get(paste(names(nucleotide)[i], sep = "_")))/10000)) - 9999),
             end = c(10000 * (1:floor(nchar(get(paste(names(nucleotide)[i], sep = "_")))/10000)),
                     nchar(get(paste(names(nucleotide)[i], sep = "_"))))),
           envir = .GlobalEnv)
  }

 for (i in 1:length(names(nucleotide))) {
    assign(paste("window", names(nucleotide)[i], sep = "_"),
           get(paste("window", names(nucleotide)[i], sep = "_")) %>% dplyr::mutate(
             "window_nt" = stringr::str_sub(get(paste(names(nucleotide)[i], sep = "_")),
                                            start = start,
                                            end = end)),
           envir = .GlobalEnv)
  }

  for (i in 1:length(names(nucleotide))) {
    assign(paste(names(nucleotide)[i], "RNA", sep = "_"),
           get(paste("window", names(nucleotide)[i], sep = "_")) %>%
             dplyr::mutate(seqid = (gff %>% dplyr::filter(source == "Local" | type == "region") %>% dplyr::select(seqid))[i, ] %>% as.character()),
           envir = .GlobalEnv)
  }

  for (i in 1:length(names(nucleotide))) {
    assign(paste("window", names(nucleotide)[i], "gc_skew", sep = "_"),
           data.table::data.table(gc_skew = gc_skew(get(paste("window", names(nucleotide)[i], sep = "_"))$window_nt[1])), # sometimes error occur
           envir = .GlobalEnv)

    for (j in (1:nrow(get(paste("window", names(nucleotide)[i], sep = "_"))))[-c(1)]) {
      addition <- setNames(as.data.frame(gc_skew(get(paste("window", names(nucleotide)[i], sep = "_"))$window_nt[j])), "gc_skew")
      assign(paste("window", names(nucleotide)[i], "gc_skew", sep = "_"), rbind(get(paste("window", names(nucleotide)[i], "gc_skew", sep = "_")), addition),
             envir = .GlobalEnv)
      rm(addition)
    }
  }

  for (i in 1:length(names(nucleotide))) {
    assign(paste("window", names(nucleotide)[i], "gc_ratio", sep = "_"),
           data.table::data.table(gc_ratio = gc_ratio(get(paste("window", names(nucleotide)[i], sep = "_"))$window_nt[1])),
           envir = .GlobalEnv)

    for (j in (1:nrow(get(paste("window", names(nucleotide)[i], sep = "_"))))[-c(1)]) {
      addition <- setNames(as.data.frame(gc_ratio(get(paste("window", names(nucleotide)[i], sep = "_"))$window_nt[j])), "gc_ratio")
      assign(paste("window", names(nucleotide)[i], "gc_ratio", sep = "_"), rbind(get(paste("window", names(nucleotide)[i], "gc_ratio", sep = "_")), addition),
             envir = .GlobalEnv)
      rm(addition)
    }
  }
  for (i in 1:length(names(nucleotide))) {
    assign(paste("window", names(nucleotide)[i], sep = "_"),
           cbind(get(paste("window", names(nucleotide)[i], sep = "_")),
                 get(paste("window", names(nucleotide)[i], "gc_skew", sep = "_")),
                 get(paste("window", names(nucleotide)[i], "gc_ratio", sep = "_"))),
           envir = .GlobalEnv)
  }

  for (i in 1:length(names(nucleotide))) {
    assign(paste("window", names(nucleotide)[i], sep = "_"),
           get(paste("window", names(nucleotide)[i], sep = "_")) %>%
             dplyr::mutate(gc_ratio_minus_average =
                             (.$gc_ratio) - (gc_ratio(get(paste(names(nucleotide)[i], sep = "_"))))),
           envir = .GlobalEnv)
  }

  for (i in 1:length(names(nucleotide))) {
    assign(paste("window", names(nucleotide)[i], sep = "_"),
           get(paste("window", names(nucleotide)[i], sep = "_")) %>%
             dplyr::mutate(gc_skew_minus_average =
                             (.$gc_skew) - (gc_skew(get(paste(names(nucleotide)[i], sep = "_"))))),
           envir = .GlobalEnv)
  }

  for (i in 1:length(names(nucleotide))) {
    assign(paste("window", names(nucleotide)[i], sep = "_"),
           get(paste("window", names(nucleotide)[i], sep = "_")) %>%
             dplyr::mutate(gc_skew_color = dplyr::case_when(gc_skew_minus_average > 0 ~ "green",
                                                            gc_skew_minus_average < 0 ~ "red")),
           envir = .GlobalEnv)
  }

  for (i in 1:length(names(nucleotide))) {
    assign(paste("window", names(nucleotide)[i], sep = "_"),
           get(paste("window", names(nucleotide)[i], sep = "_")) %>%
             dplyr::mutate(gc_ratio_color = dplyr::case_when(gc_ratio_minus_average > 0 ~ "blue",
                                                             gc_ratio_minus_average < 0 ~ "yellow")),
           envir = .GlobalEnv)
  }
  assign("gff", gff, envir = .GlobalEnv)
}


