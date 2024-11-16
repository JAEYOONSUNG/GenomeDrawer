#' RunDrawer
#'
#' @param completness A genbank completeness of data.
#' @return A genome map.
#' @export
#'

run_drawer <- function(completeness = NULL) {
  if (is.null(completeness)) {
    completeness <- FALSE
  }

  # Call Function Genbank_organizer
  Genbank_organizer()
  message("Step1. genbank_organizer function executed")

  # Call Genbank fna extractor Function
  genbank_fna_extractor()
  message("Step2. Genbank fna extractor function executed")

  # Step 3: Call EggNOG_annotations and merge with genbank_table
  # Check if 'emapper.annotations' exists in the current working directory
  if (any(grepl("\\.emapper\\.annotations\\.", list.files(getwd())))) {
    # Step 3: Call EggNOG_annotations and merge with genbank_table
    EggNOG_annotations()
    message("Step3. EggNOG_annotations function executed")

    # Check if EggNOG_table is NULL and exists in the environment
    if (is.null(EggNOG_table)) {
      if (exists("EggNOG_table", envir = .GlobalEnv)) {
        EggNOG_table <- get("EggNOG_table", envir = .GlobalEnv)
        message("EggNOG_table fetched from the global environment.")
      } else {
        message("EggNOG_table not found in the environment. Skipping this step.")
        return(invisible(NULL))  # Skip this step and continue
      }
    }

    # Check if the necessary columns exist in the tables before merging
    if (!"locus_tag" %in% colnames(genbank_table)) {
      message("Error: 'locus_tag' not found in genbank_table. Skipping the merge step.")
      return(invisible(NULL))
    }
    if (!"query" %in% colnames(EggNOG_table)) {
      message("Error: 'query' not found in EggNOG_table. Skipping the merge step.")
      return(invisible(NULL))
    }

    # Merge the genbank_table with EggNOG_table
    genbank_table <- merge(genbank_table, EggNOG_table, by.x = "locus_tag", by.y = "query", all.x = TRUE)

    # Assign the result to 'genbank_table' in the global environment
    assign("genbank_table", genbank_table, envir = .GlobalEnv)
    message("The result has been saved to the R environment variable 'genbank_table'.")
  } else {
    message("No 'emapper.annotations' file found in the working directory. Skipping EggNOG_annotations.")
  }


  # Merge the genbank_table with EggNOG_table
  gb_info()
  message("Genbank info has been saved to the R environment variable 'Genome_summary.")

  # Assign the result to 'genbank_table' in the global environment
  assign("genbank_table", genbank_table, envir = .GlobalEnv)
  # Call Genbank fna extractor Function
  GenomeDrawer(dir = NULL)
  message("Step4. factor extractor from gff")

  # Conditional step based on 'completeness'
  if (completeness == TRUE) {
    GenomeVisualization()
    message("Step5. Genome map visualization for completeness.")
  } else {
    ContigsVisualization()
    message("Step5. Contigs visualization for non-completeness.")
  }
}
