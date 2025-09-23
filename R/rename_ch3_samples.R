#' Rename sample names in a table
#'
#' Update \code{sample_name} values in a specified DuckDB table using a mapping
#' from old names to new names. Accepts either a named character vector
#' (\code{c(old="new", ...)}), or a two-column data.frame with columns
#' \code{old} and \code{new}.
#'
#' @param ch3_db Path to a \code{.ch3.db} file or a \code{"ch3_db"} object.
#' @param table Character scalar: the table to modify (e.g. \code{"positions"},
#'   \code{"windows"}, \code{"regions"}, or \code{"calls"}).
#' @param samples_map Either a named character vector (\code{names = old, values = new})
#'   or a data.frame with columns \code{old} and \code{new}.
#' @param strict Logical; if \code{TRUE} (default) stop on issues (missing old
#'   names, empty/NA new names, or duplicate mappings). If \code{FALSE}, tries
#'   to proceed after dropping invalid rows with a warning.
#' @param preview Logical; if \code{TRUE}, prints a small before/after summary
#'   of distinct sample names.
#'
#' @return (Invisibly) the updated \code{"ch3_db"} object (connection closed on return).
#'
#' @examples
#' \dontrun{
#' # Named character vector
#' rename_ch3_samples("my_db.ch3.db", table = "positions",
#'                    samples_map = c("Astrocytes" = "Astro",
#'                                    "Cortical_Neurons" = "Cortical"))
#'
#' # data.frame mapping
#' m <- data.frame(old = c("Ctrl1","Ctrl2"), new = c("Control_1","Control_2"))
#' rename_ch3_samples("my_db.ch3.db", "windows", m)
#' }
#'
#' @importFrom DBI dbExistsTable dbGetQuery dbWriteTable dbQuoteIdentifier dbExecute
#' @export
rename_ch3_samples <- function(ch3_db,
                               table,
                               samples_map,
                               strict = TRUE,
                               preview = TRUE)
{
  # ---- normalize mapping --------------------------------------------------------
  if (is.character(samples_map) && !is.null(names(samples_map))) {
    old <- names(samples_map)
    new <- unname(samples_map)
    map <- data.frame(old = old, new = new, stringsAsFactors = FALSE)
  } else if (is.data.frame(samples_map) && all(c("old","new") %in% names(samples_map))) {
    map <- data.frame(old = as.character(samples_map$old),
                      new = as.character(samples_map$new),
                      stringsAsFactors = FALSE)
  } else {
    stop("samples_map must be a named character vector OR a data.frame with columns 'old' and 'new'.")
  }
  
  # Basic hygiene
  map$old[is.na(map$old)] <- ""
  map$new[is.na(map$new)] <- ""
  bad <- (map$old == "") | (map$new == "")
  if (any(bad)) {
    msg <- sprintf("Found %d mapping rows with empty old/new.", sum(bad))
    if (strict) stop(msg)
    warning(msg, " Dropping those rows.")
    map <- map[!bad, , drop = FALSE]
  }
  
  # Deduplicate identical old→new pairs; check for duplicate 'old'
  map <- unique(map)
  dup_old <- duplicated(map$old)
  if (any(dup_old)) {
    msg <- paste0("Duplicate 'old' names in mapping: ",
                  paste(unique(map$old[dup_old]), collapse = ", "))
    if (strict) stop(msg)
    warning(msg, " Keeping the first occurrence.")
    map <- map[!duplicated(map$old), , drop = FALSE]
  }
  
  # Optional: warn if some 'new' names collide among themselves (not fatal)
  if (any(duplicated(map$new))) {
    warning("Multiple old names map to the same new name: ",
            paste(unique(map$new[duplicated(map$new)]), collapse = ", "),
            ". Resulting table will have merged samples under that name.")
  }
  
  # ---- open DB, checks ----------------------------------------------------------
  ch3_db <- MethylSeqR:::.ch3helper_connectDB(ch3_db)
  on.exit({ MethylSeqR:::.ch3helper_closeDB(ch3_db) }, add = TRUE)
  
  if (!DBI::dbExistsTable(ch3_db$con, table)) {
    stop("Table '", table, "' does not exist in the database.")
  }
  
  tbl_id <- as.character(DBI::dbQuoteIdentifier(ch3_db$con, table))
  
  # Ensure table has a sample_name column
  cols <- DBI::dbGetQuery(ch3_db$con, sprintf("PRAGMA table_info(%s);", tbl_id))$name
  if (!"sample_name" %in% cols) {
    stop("Table '", table, "' has no 'sample_name' column.")
  }
  
  # Check old names exist
  existing <- DBI::dbGetQuery(ch3_db$con,
                              sprintf("SELECT DISTINCT sample_name FROM %s;", tbl_id))$sample_name
  missing_old <- setdiff(map$old, existing)
  if (length(missing_old)) {
    msg <- paste0("These 'old' sample names were not found in '", table, "': ",
                  paste(missing_old, collapse = ", "))
    if (strict) stop(msg)
    warning(msg, " They will be ignored.")
    map <- map[!map$old %in% missing_old, , drop = FALSE]
  }
  
  if (!nrow(map)) {
    message("No valid mappings to apply; nothing changed.")
    ch3_db$current_table <- table
    return(invisible(ch3_db))
  }
  
  # ---- preview before (single stream) ------------------------------------------
  message("")
  if (preview) {
    message("Before rename: distinct sample names in '", table, "':")
    if (length(existing)) {
      show <- sort(existing)[1:min(20, length(existing))]
      message(paste(show, collapse = "\n"))
      if (length(existing) > 20) message("... (showing first 20)")
    } else {
      message("(none)")
    }
  }
  
  # ---- write temp mapping + update ---------------------------------------------
  DBI::dbExecute(ch3_db$con, "DROP TABLE IF EXISTS temp_sample_map;")
  DBI::dbWriteTable(ch3_db$con, "temp_sample_map", map, temporary = TRUE, overwrite = TRUE)
  
  # Perform UPDATE with join (DuckDB supports UPDATE ... FROM)
  sql_update <- sprintf("
    UPDATE %s AS t
    SET sample_name = m.new
    FROM temp_sample_map AS m
    WHERE t.sample_name = m.old;
  ", tbl_id)
  
  n_changed <- DBI::dbExecute(ch3_db$con, sql_update)
  DBI::dbExecute(ch3_db$con, "DROP TABLE IF EXISTS temp_sample_map;")
  
  # ---- preview after (single stream) -------------------------------------------
  if (preview) {
    after <- DBI::dbGetQuery(ch3_db$con,
                             sprintf("SELECT DISTINCT sample_name FROM %s;", tbl_id))$sample_name
    message("\nAfter rename: distinct sample names in '", table, "':")
    if (length(after)) {
      show <- sort(after)[1:min(20, length(after))]
      message(paste(show, collapse = "\n"))
      if (length(after) > 20) message("... (showing first 20)")
    } else {
      message("(none)")
    }
    message("\nRows updated: ", n_changed)
  }
  
  ch3_db$current_table <- table
  invisible(ch3_db)
}
