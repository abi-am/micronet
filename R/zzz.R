###############################################################
#########   fixing errors in ANCOMBC and Phylosmith   #########
###############################################################

.onLoad <- function(libname, pkgname) {
  # Override co_occurrence with your modified version
  utils::assignInNamespace(".bias_em", .bias_em_patched, ns = asNamespace("ANCOMBC"))
  utils::assignInNamespace("co_occurrence", modified_co_occurrence, ns = asNamespace("phylosmith"))
}
