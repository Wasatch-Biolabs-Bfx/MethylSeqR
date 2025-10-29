# -------------------------------------------------------------
# Helpers
# -------------------------------------------------------------

# Auto-detect resource caps: all-but-one cores and ~80% of RAM
.auto_duckdb_resource_caps <- function(memory_fraction = 0.80) {
  cores <- tryCatch(parallel::detectCores(), error = function(e) 2)
  threads <- max(1, cores - 1)
  
  total_bytes <- NA_real_
  sysname <- Sys.info()[["sysname"]]
  if (file.exists("/proc/meminfo")) {
    ln <- tryCatch(readLines("/proc/meminfo"), error = function(e) "")
    mt <- sub(".*MemTotal:\\s*([0-9]+) kB.*", "\\1", grep("^MemTotal:", ln, value = TRUE))
    if (nzchar(mt)) total_bytes <- as.numeric(mt) * 1024
  } else if (identical(sysname, "Darwin")) {
    out <- suppressWarnings(system("sysctl -n hw.memsize", intern = TRUE))
    if (length(out)) total_bytes <- as.numeric(out)
  } else if (.Platform$OS.type == "windows") {
    out <- suppressWarnings(system("wmic computersystem get TotalPhysicalMemory /value", intern = TRUE))
    if (length(out)) {
      val <- sub(".*TotalPhysicalMemory=([0-9]+).*", "\\1", paste(out, collapse = ""))
      if (grepl("^[0-9]+$", val)) total_bytes <- as.numeric(val)
    }
    if (!is.finite(total_bytes)) {
      ps <- suppressWarnings(system('powershell -Command "(Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory"', intern = TRUE))
      if (length(ps) && grepl("^[0-9]+$", ps[1])) total_bytes <- as.numeric(ps[1])
    }
  }
  if (!is.finite(total_bytes)) total_bytes <- 8 * 1024^3  # fallback 8GB
  
  limit_bytes <- floor(memory_fraction * total_bytes)
  limit_mb    <- max(1024, floor(limit_bytes / 1024^2))   # >= 1GB
  mem_str     <- paste0(limit_mb, "MB")
  list(threads = threads, memory_limit = mem_str)
}


# Parse mod_code like c("m","h","m + h","a"); returns data.frame(label, codes(list))
.parse_mod_specs <- function(mod_code) {
  if (length(mod_code) == 0) return(data.frame(label=character(), codes=I(list())))
  specs <- lapply(mod_code, function(s) {
    s2    <- gsub("\\s+", "", s)
    parts <- strsplit(s2, "\\+", fixed = FALSE)[[1]]
    label <- paste0(parts, collapse = "")
    list(label = label, codes = parts)
  })
  data.frame(
    label = vapply(specs, `[[`, character(1), "label"),
    codes = I(lapply(specs, `[[`, "codes")),
    stringsAsFactors = FALSE
  )
}




# # Build SQL snippets for per-position counts for unmodified and mods/combinations
# .build_pos_count_sql <- function(unmod_code, unmod_label, specs_df) {
#   esc <- function(x) gsub("'", "''", x)
#   pieces <- c(sprintf("SUM(CASE WHEN call_code = '%s' THEN 1 ELSE 0 END) AS %s_counts",
#                       esc(unmod_code), unmod_label))
#   for (i in seq_len(nrow(specs_df))) {
#     lab <- specs_df$label[i]
#     codes <- specs_df$codes[[i]]
#     if (length(codes) == 1) {
#       pieces <- c(pieces, sprintf(
#         "SUM(CASE WHEN call_code = '%s' THEN 1 ELSE 0 END) AS %s_counts",
#         esc(codes), lab))
#     } else {
#       pieces <- c(pieces, sprintf(
#         "SUM(CASE WHEN call_code IN ('%s') THEN 1 ELSE 0 END) AS %s_counts",
#         paste(esc(codes), collapse = "','"), lab))
#     }
#   }
#   labels_all <- c(unmod_label, specs_df$label)
#   list(select_counts_pos = paste(pieces, collapse = ",\n          "),
#        labels_all = labels_all)
# }

# Build SQL snippets for per-position counts for unmodified and mods/combinations
.build_pos_count_sql <- function(unmod_code, unmod_label, specs_df) {
  esc <- function(x) gsub("'", "''", x)
  
  # --- helpers to make identifiers safe only if numeric ---
  .is_numeric_code <- function(x) grepl("^[0-9]+$", x)
  .make_safe_label <- function(x) if (.is_numeric_code(x)) paste0("m_", x) else x
  
  pieces <- character()
  
  # Unmodified base (make it safe too just in case)
  unmod_lab_raw  <- unmod_label
  unmod_lab_safe <- .make_safe_label(unmod_lab_raw)
  pieces <- c(pieces,
              sprintf("SUM(CASE WHEN call_code = '%s' THEN 1 ELSE 0 END) AS %s_counts",
                      esc(unmod_code), unmod_lab_safe)
  )
  
  # Modified codes and combinations
  for (i in seq_len(nrow(specs_df))) {
    lab_raw   <- specs_df$label[i]
    lab_safe  <- .make_safe_label(lab_raw)
    codes_vec <- specs_df$codes[[i]]
    
    if (length(codes_vec) == 1) {
      pieces <- c(pieces, sprintf(
        "SUM(CASE WHEN call_code = '%s' THEN 1 ELSE 0 END) AS %s_counts",
        esc(codes_vec), lab_safe))
    } else {
      pieces <- c(pieces, sprintf(
        "SUM(CASE WHEN call_code IN ('%s') THEN 1 ELSE 0 END) AS %s_counts",
        paste(esc(codes_vec), collapse = "','"), lab_safe))
    }
  }
  
  # Return both the SQL and the safe labels
  labels_raw  <- c(unmod_lab_raw,  specs_df$label)
  labels_safe <- c(unmod_lab_safe, vapply(specs_df$label, .make_safe_label, character(1)))
  
  list(
    select_counts_pos = paste(pieces, collapse = ",\n          "),
    labels_all        = labels_safe,                         # safe labels for downstream use
    label_map         = setNames(labels_safe, labels_raw)     # optional: map raw -> safe
  )
}






# Build a chrom filter clause if chrs provided (character vector) else ""
.chrom_filter_sql <- function(chrs) {
  if (is.null(chrs) || length(chrs) == 0) return("")
  vals <- paste(sprintf("'%s'", gsub("'", "''", chrs)), collapse = ", ")
  paste0(" AND chrom IN (", vals, ")")
}
