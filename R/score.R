#' Score target sequences against a weight matrix
#'
#' Implements the paper's linear regression-based scoring: the final score for
#' a sequence is the sum of the per-position weights of the bases it carries,
#' over the positions retained in the weight matrix.
#'
#' Higher scores indicate a greater predicted likelihood of being functional,
#' and therefore *lower* expected target gene expression.
#'
#' @param sequences Character vector of target sequences. Each must be the same
#'   length as the sequences used to build `wm` (50 nucleotides for the bundled
#'   training set: a 20-nucleotide targeting site flanked by 15 nucleotides on
#'   each side). `T` is accepted and read as `U`.
#' @param wm An `sdrna_wm` weight matrix from [sdrna_weight_matrix()].
#'
#' @return A numeric vector of scores, named with `sequences` when those are
#'   named.
#'
#' @examples
#' wm <- sdrna_weight_matrix(seed = 1)
#' scores <- sdrna_score(sdrna356$sequence, wm)
#' # Higher score should track lower target expression.
#' cor(scores, sdrna356$expression)
#' @export
sdrna_score <- function(sequences, wm) {
  if (!inherits(wm, "sdrna_wm")) {
    stop("`wm` must be an sdrna_wm object from sdrna_weight_matrix().",
         call. = FALSE)
  }
  m <- sdrna_split(sequences)
  if (ncol(m) != ncol(wm)) {
    stop("Sequences are ", ncol(m), " nucleotides but the weight matrix ",
         "covers ", ncol(wm), ".", call. = FALSE)
  }
  positions <- attr(wm, "positions")
  idx <- match(m[, positions, drop = FALSE], SDRNA_BASES)
  idx <- matrix(idx, nrow = nrow(m))
  # Column j of the subset maps to weight-matrix column positions[j].
  w <- vapply(seq_along(positions),
              function(j) wm[idx[, j], positions[j]],
              numeric(nrow(m)))
  if (is.null(dim(w))) w <- matrix(w, nrow = nrow(m))
  out <- rowSums(w)
  names(out) <- names(sequences)
  out
}
