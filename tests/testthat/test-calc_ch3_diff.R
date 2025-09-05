test_that("calc_ch3_diff() creates mod_diff_positions with expected fields and values", {
  skip_if_not_installed("duckdb")
  skip_if_not_installed("DBI")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tibble")
  
  library(DBI)
  library(dplyr)
  library(tibble)
  
  # ---- Minimal 2-locus, 2v2 design
  df <- tibble(
    sample_name = c("case1","case2","ctrl1","ctrl2","case1","case2","ctrl1","ctrl2"),
    chrom       = c("chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1"),
    start       = c(100,100,100,100,200,200,200,200),
    end         = c(101,101,101,101,201,201,201,201),
    num_calls   = c(10, 12, 11,  9,  8, 10, 12,  7),
    mh_counts   = c( 7,  8,  3,  4,  2,  2,  6,  3)
  )
  # sanity: locus 100 → meth_diff > 0; locus 200 → meth_diff < 0
  
  tmpdir <- withr::local_tempdir()
  dbfile <- file.path(tmpdir, "tiny.ch3.db")
  con    <- DBI::dbConnect(duckdb::duckdb(dbfile, read_only = FALSE))
  DBI::dbWriteTable(con, "positions", df)
  
  # Minimal ch3_db object per your structure
  ch3_db_obj <- list(db_file = dbfile, current_table = NULL, con = con)
  class(ch3_db_obj) <- "ch3_db"
  
  testthat::with_mocked_bindings(
    # Use the open connection we already made; don't defer close.
    .ch3helper_connectDB = function(x) x,
    # Do not purge/vacuum/close during test; return as-is so we can inspect results.
    .ch3helper_cleanup   = function(x) x,
    {
      # ---------- fast_fisher
      msg <- capture.output(
        res_fast <- calc_ch3_diff(
          ch3_db    = ch3_db_obj,
          call_type = "positions",
          cases     = c("case1","case2"),
          controls  = c("ctrl1","ctrl2"),
          mod_type  = "mh",
          calc_type = "fast_fisher"
        )
      )
      
      expect_s3_class(res_fast, "ch3_db")
      expect_true(DBI::dbExistsTable(con, "mod_diff_positions"))
      
      out <- dplyr::tbl(con, "mod_diff_positions") |> arrange(chrom, start) |> collect()
      
      # core columns
      expect_true(all(c("p_val", "p_adjust", "meth_diff") %in% names(out)))
      # mod_type rename happened: "mod_*" -> "mh_*"
      expect_true(any(grepl("^mh_counts_", names(out))), info = "mod_type rename to 'mh_' should occur")
      
      # BH-adjusted p-values in [0,1]
      expect_true(all(out$p_adjust >= 0 & out$p_adjust <= 1, na.rm = TRUE))
      
      # meth_diff directions by locus
      expect_gt(out$meth_diff[out$start == 100], 0)
      expect_lt(out$meth_diff[out$start == 200], 0)
      
      # ---------- r_fisher should closely match fast_fisher p-values
      DBI::dbRemoveTable(con, "mod_diff_positions")
      invisible(calc_ch3_diff(
        ch3_db    = ch3_db_obj,
        call_type = "positions",
        cases     = c("case1","case2"),
        controls  = c("ctrl1","ctrl2"),
        mod_type  = "mh",
        calc_type = "r_fisher"
      ))
      out_r <- dplyr::tbl(con, "mod_diff_positions") |> arrange(chrom, start) |> collect()
      
      DBI::dbRemoveTable(con, "mod_diff_positions")
      invisible(calc_ch3_diff(
        ch3_db    = ch3_db_obj,
        call_type = "positions",
        cases     = c("case1","case2"),
        controls  = c("ctrl1","ctrl2"),
        mod_type  = "mh",
        calc_type = "fast_fisher"
      ))
      out_fast <- dplyr::tbl(con, "mod_diff_positions") |> arrange(chrom, start) |> collect()
      
      expect_equal(order(out_fast$p_val), order(out_r$p_val))
    }
  )
  
  DBI::dbDisconnect(con, shutdown = TRUE)
})

test_that("calc_ch3_diff() supports different mod_type renames", {
  skip_if_not_installed("duckdb"); skip_if_not_installed("DBI")
  skip_if_not_installed("dplyr");  skip_if_not_installed("tibble")
  
  tmpdir <- withr::local_tempdir()
  dbfile <- file.path(tmpdir, "tiny2.ch3.db")
  con    <- DBI::dbConnect(duckdb::duckdb(dbfile, read_only = FALSE))
  
  df <- tibble::tibble(
    sample_name = c("case","ctrl"),
    chrom="chr1", start=1L, end=2L,
    num_calls=10L, mh_counts=6L
  )
  DBI::dbWriteTable(con, "positions", df)
  ch3_db_obj <- list(db_file = dbfile, current_table = NULL, con = con)
  class(ch3_db_obj) <- "ch3_db"
  
  testthat::with_mocked_bindings(
    .ch3helper_connectDB = function(x) x,
    .ch3helper_cleanup   = function(x) x,
    {
      invisible(calc_ch3_diff(ch3_db_obj, "positions",
                              cases="case", controls="ctrl",
                              mod_type = "h", calc_type="fast_fisher"))
      cols <- colnames(dplyr::tbl(con, "mod_diff_positions"))
      # "mod_*" → "h_*"
      expect_true(any(grepl("^h_counts_", cols)))
      expect_false(any(grepl("^mod_counts_", cols)))
    }
  )
  
  DBI::dbDisconnect(con, shutdown = TRUE)
})

test_that("calc_ch3_diff() errors clearly for missing table or sample names", {
  skip_if_not_installed("duckdb"); skip_if_not_installed("DBI"); skip_if_not_installed("tibble")
  
  tmpdir <- withr::local_tempdir()
  dbfile <- file.path(tmpdir, "tiny3.ch3.db")
  con    <- DBI::dbConnect(duckdb::duckdb(dbfile, read_only = FALSE))
  ch3_db_obj <- list(db_file = dbfile, current_table = NULL, con = con)
  class(ch3_db_obj) <- "ch3_db"
  
  testthat::with_mocked_bindings(
    .ch3helper_connectDB = function(x) x,
    .ch3helper_cleanup   = function(x) x,
    {
      # No positions table
      expect_error(
        calc_ch3_diff(ch3_db_obj, "positions", cases="a", controls="b"),
        regexp = "table does not exist", fixed = FALSE
      )
      
      # Only controls present
      DBI::dbWriteTable(con, "positions",
                        tibble::tibble(sample_name=c("ctrl1","ctrl2"),
                                       chrom="chr1", start=1L, end=2L,
                                       num_calls=10L, mh_counts=5L)
      )
      expect_error(
        calc_ch3_diff(ch3_db_obj, "positions", cases="case1", controls=c("ctrl1","ctrl2")),
        regexp = "^Check case names - some case samples are missing", fixed = FALSE
      )
      
      DBI::dbRemoveTable(con, "positions")
      DBI::dbWriteTable(con, "positions",
                        tibble::tibble(sample_name=c("case1","case2"),
                                       chrom="chr1", start=1L, end=2L,
                                       num_calls=10L, mh_counts=5L)
      )
      expect_error(
        calc_ch3_diff(ch3_db_obj, "positions", cases=c("case1","case2"), controls="ctrl1"),
        regexp = "^Check control names- some control samples are not in the data", fixed = FALSE
      )
    }
  )
  
  DBI::dbDisconnect(con, shutdown = TRUE)
})

test_that(".fast_fisher matches fisher.test() on simple 2x2s", {
  skip_if_not_installed("stats")
  
  # Access internal via ::: since it's not exported
  ff <- MethylSeqR:::`.fast_fisher`
  rf <- function(a,b,c,d) fisher.test(matrix(c(a,b,c,d), nrow = 2))$p.value
  
  mats <- list(
    c(7, 3, 3, 7),   # symmetric
    c(8, 2, 4, 6),
    c(10,0,0,10),
    c(12,5,8,15)
  )
  for (m in mats) {
    # Map to q,m,n,k per your implementation:
    # q = mod_case, m = mod_case + mod_ctrl, n = unmod_case + unmod_ctrl, k = num_case
    a <- m[1]; b <- m[2]; c <- m[3]; d <- m[4]
    q <- b
    mpar <- a + b
    npar <- c + d
    kpar <- b + d
    expect_equal(ff(q, mpar, npar, kpar), rf(a,b,c,d), tolerance = 1e-8)
  }
})
