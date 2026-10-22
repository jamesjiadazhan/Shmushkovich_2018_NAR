#' Positive predictive power versus sensitivity
#'
#' Reproduces the performance curves of the paper's Figure 2B and Figure 3.
#' Sequences are ranked by score and the threshold is swept from the
#' highest-scoring compound downwards. At each threshold:
#'
#' * positive predictive power (PPP) is the percent of selected compounds that
#'   are genuinely functional;
#' * sensitivity is the percent of all functional compounds that have been
#'   selected.
#'
#' @param scores Numeric vector of scores, as returned by [sdrna_score()].
#' @param expression Numeric vector of measured target gene expression, percent
#'   of control, the same length as `scores`.
#' @param functional_cutoff Expression value below which a compound counts as
#'   functional.
#'
#' @return A data frame with one row per threshold and columns `threshold`,
#'   `n_selected`, `n_true`, `ppp` and `sensitivity`, the last two as percents.
#'
#' @examples
#' wm <- sdrna_weight_matrix(seed = 1)
#' s <- sdrna_score(sdrna356$sequence, wm)
#' curve <- sdrna_ppp_sensitivity(s, sdrna356$expression, 24)
#' head(curve)
#' @export
sdrna_ppp_sensitivity <- function(scores, expression, functional_cutoff = 24) {
  stopifnot(length(scores) == length(expression))
  ord <- order(scores, decreasing = TRUE)
  truth <- expression[ord] < functional_cutoff
  n_func <- sum(expression < functional_cutoff)
  k <- seq_along(scores)
  tp <- cumsum(truth)
  data.frame(
    threshold = scores[ord],
    n_selected = k,
    n_true = tp,
    ppp = 100 * tp / k,
    sensitivity = if (n_func > 0) 100 * tp / n_func else rep(NA_real_, length(k))
  )
}

#' Positive predictive power at a target sensitivity
#'
#' Reads a single operating point off the curve from
#' [sdrna_ppp_sensitivity()], the way the paper quotes figures such as
#' "85 percent accuracy with 40 percent sensitivity".
#'
#' @param curve A data frame from [sdrna_ppp_sensitivity()].
#' @param sensitivity Target sensitivity, as a percent.
#'
#' @return A one-row data frame: the first threshold reaching at least
#'   `sensitivity`.
#' @export
sdrna_ppp_at <- function(curve, sensitivity = 25) {
  hit <- which(curve$sensitivity >= sensitivity)
  if (!length(hit)) {
    return(curve[0, , drop = FALSE])
  }
  curve[hit[1L], , drop = FALSE]
}

#' Area under the ROC curve
#'
#' Computes AUC by the rank (Mann-Whitney) identity, with ties averaged, as
#' used for the paper's Supplementary Figure 2.
#'
#' @param scores Numeric vector of scores.
#' @param labels Logical vector, `TRUE` for functional compounds.
#'
#' @return A single numeric AUC value, or `NA_real_` when either class is
#'   empty.
#'
#' @examples
#' wm <- sdrna_weight_matrix(seed = 1)
#' s <- sdrna_score(sdrna356$sequence, wm)
#' sdrna_auc(s, sdrna356$expression < 24)
#' @export
sdrna_auc <- function(scores, labels) {
  stopifnot(length(scores) == length(labels))
  labels <- as.logical(labels)
  pos <- scores[labels]
  neg <- scores[!labels]
  if (!length(pos) || !length(neg)) return(NA_real_)
  r <- rank(c(pos, neg))
  (sum(r[seq_along(pos)]) - length(pos) * (length(pos) + 1) / 2) /
    (length(pos) * length(neg))
}

#' ROC curve points
#'
#' @param scores Numeric vector of scores.
#' @param labels Logical vector, `TRUE` for functional compounds.
#'
#' @return A data frame with columns `threshold`, `fpr` and `tpr`.
#' @export
sdrna_roc <- function(scores, labels) {
  stopifnot(length(scores) == length(labels))
  labels <- as.logical(labels)
  ord <- order(scores, decreasing = TRUE)
  truth <- labels[ord]
  tpr <- cumsum(truth) / sum(truth)
  fpr <- cumsum(!truth) / sum(!truth)
  data.frame(threshold = scores[ord],
             fpr = c(0, fpr),
             tpr = c(0, tpr))[seq_len(length(scores) + 1L), ]
}
