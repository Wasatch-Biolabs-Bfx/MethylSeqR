qc_modseq <- function(modseq_dat, plot = TRUE) {
  print("calculating coverage stats...")
  get_cov_stats(modseq_dat = modseq_dat, plot = plot)
  print("calculating mod stats...")
  get_mod_stats(modseq_dat = modseq_dat, plot = plot)
  print("calculating correlations...")
  cor_modseq(modseq_dat = modseq_dat, plot = plot)
  print("calculating cluster...")
  cluster_modseq(modseq_dat = modseq_dat)
  print("running pca...")
  pca_modseq(modseq_dat = modseq_dat)
  print("calculating methylation density...")
  density_modseq(modseq_dat = modseq_dat)
}
