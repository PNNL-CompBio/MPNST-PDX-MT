#' Utilities for importing MSigDB gene sets into leapR
#'
#' @description
#' Functions for listing MSigDB collections, retrieving a selected collection,
#' and converting gene sets into a format readable by leapR.
#'
#' @details
#' \itemize{
#'   \item [list_msigdb_genesets()] lists collections available for a database
#'     species.
#'   \item [get_genesets_from_msigdb()] retrieves a selected MSigDB collection
#'     and converts it to a leapR gene-set object.
#'   \item [.transfer_to_leapr()] is an internal helper that converts a named
#'     gene list to a leapR gene-set object using a temporary GMT file.
#' }
#'
#' @seealso
#' [list_msigdb_genesets()], [get_genesets_from_msigdb()],
#' [.transfer_to_leapr()]
#'
#' @name msigdbr_to_leapr

#' List available MSigDB gene set collections
#'
#' Wrapper for [msigdbr::msigdbr_collections()] that returns the available
#' MSigDB collections for the specified database species.
#'
#' @param species Database species code passed to
#'   [msigdbr::msigdbr_collections()]. Defaults to `"HS"`.
#'
#' @return A data frame of available MSigDB collections.
#' @export
list_msigdb_genesets <- function(
  species = "HS"
) {
  msigdbr::msigdbr_collections(db_species = species)
}

#' Retrieve an MSigDB gene set collection
#'
#' Retrieves gene sets from MSigDB using [msigdbr::msigdbr()] and converts them
#' to a gene-set object readable by leapR.
#'
#' @param species Species name passed to [msigdbr::msigdbr()].
#'   Defaults to `"Homo sapiens"`.
#' @param msigdb_subcollection MSigDB subcollection passed to
#'   [msigdbr::msigdbr()]. Available subcollections can be found with
#'   [list_msigdb_genesets()]. Defaults to `"CP:PID"`.
#'
#' @return A leapR gene-set object containing the retrieved gene sets.
#'
#' @examples
#' \dontrun{
#' pid_genesets <- get_genesets_from_msigdb()
#' }
get_genesets_from_msigdb <- function(
  species = "Homo sapiens",
  msigdb_subcollection = "CP:PID"
) {
  msigdbr_df <- msigdbr::msigdbr(
    species = species,
    subcollection = msigdb_subcollection
  ) %>%
    select(gs_name, gene_symbol) %>%
    distinct()
  msigdbr_list <- split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)

  geneset <- .transfer_to_leapr(msigdbr_list)
  return(geneset)
}

#' Convert a named gene list to a leapR gene-set object
#'
#' Writes the gene sets to a temporary GMT file, reads them with
#' [leapR::read_gene_sets()], and transposes the resulting matrix.
#'
#' @param gene_list A named list of character vectors. Each list name is a gene
#'   set name, and each vector contains that set's gene symbols.
#'
#' @return A leapR gene-set object containing the supplied gene sets.
#'
#' @examples
#' \dontrun{
#' gene_list <- list(example_set = c("TP53", "EGFR"))
#' geneset <- .transfer_to_leapr(gene_list)
#' }
#' @keywords internal
.transfer_to_leapr <- function(gene_list) {
  tmp_file <- tempfile(fileext = ".gmt")
  con <- file(tmp_file, "w")

  for (set_name in names(gene_list)) {
    genes <- gene_list[[set_name]]
    # GMT format requires: Name, Description (or NA), followed by genes
    row_data <- c(set_name, "NA", genes)
    # Write as a tab-separated line
    writeLines(paste(row_data, collapse = "\t"), con)
  }
  close(con)

  geneset <- leapR::read_gene_sets(tmp_file)
  unlink(tmp_file)

  geneset$matrix <- t(geneset$matrix)
  return(geneset)
}
