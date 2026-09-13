# Shared spatial workflow validation and default signatures.
.biomed_vis_score_dataframe <- function(x) {
  if (inherits(x, "Seurat")) return(x[[]])
  if (is.data.frame(x) || is.matrix(x)) return(as.data.frame(x, check.names = FALSE))
  stop("Supply a Seurat object, score matrix or metadata data frame.")
}

.biomed_vis_positive_integer <- function(x, name, minimum = 1L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < minimum || x != floor(x))
    stop(name, " must be an integer >= ", minimum, ".")
  invisible(x)
}

.biomed_vis_neural_pathways <- function() list(
  Adrenergic = c("ADRA1A","ADRA1B","ADRA1D","ADRA2A","ADRA2B","ADRA2C","ADRB1","ADRB2","ADRB3"),
  CGRP = c("CALCRL","RAMP1","RAMP2","RAMP3","CALCR","GNAS","ADCY3","ADCY6","PRKACA","CREB1"),
  Tachykinin = c("TACR1","TACR2","TACR3","TAC1","TAC3","TAC4","GNAQ","GNA11","PLCB1","PRKCA"),
  NPY = c("NPY1R","NPY2R","NPY4R","NPY5R","NPY","PYY","PPY","GNAI1","GNAI2","GNAI3"),
  Cholinergic = c("CHRNA7","CHRNA5","CHRNA9","CHRNA10","CHRM1","CHRM2","CHRM3","CHRM4","CHRM5","ACHE","SLC5A7"),
  Neurotrophin = c("NTRK1","NTRK2","NTRK3","NGFR","NGF","BDNF","NTF3","NTF4","SORT1","NGFRAP1"),
  Glutamatergic = c("GRIA1","GRIA2","GRIA3","GRIA4","GRIN1","GRIN2A","GRIN2B","GRM1","GRM3","GRM5"),
  GABAergic = c("GABRA1","GABRA2","GABRA3","GABRB1","GABRB2","GABRB3","GABRG1","GABRG2","GABBR1","GABBR2"),
  Serotonergic = c("HTR1A","HTR1B","HTR1D","HTR2A","HTR2B","HTR2C","HTR3A","HTR4","HTR6","HTR7"),
  Dopaminergic = c("DRD1","DRD2","DRD3","DRD4","DRD5","SLC6A3","DDC","COMT","GNAS","GNAI2"),
  Neural_Remodeling = c("UCHL1","TUBB3","PRPH","NEFL","NEFM","GAP43","SNAP25","SYP","NRXN1","NCAM1"),
  Axon_Guidance = c("SEMA3A","SEMA3C","SEMA3F","NRP1","NRP2","SLIT2","ROBO1","EFNA1","EPHA2","NGF","BDNF")
)
