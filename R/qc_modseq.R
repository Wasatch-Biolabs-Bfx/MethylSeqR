qc_modseq <- function(modseq_dat, plot = TRUE) {
  get_cov_stats(modseq_dat = modseq_dat, plot = plot)
  get_mod_stats(modseq_dat = modseq_dat, plot = plot)
  cor_modseq(modseq_dat = modseq_dat, plot = plot)
  cluster_modseq(modseq_dat = modseq_dat, plot = plot)
  pca_modseq(modseq_dat = modseq_dat)
  density_modseq(modseq_dat = modseq_dat)
}
