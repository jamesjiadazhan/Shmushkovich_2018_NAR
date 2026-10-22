#' Cross-validated performance of the scoring algorithm
#'
#' The performance figures quoted in the paper are resubstitution estimates:
#' the weight matrix is built from all 356 training compounds and then scored
#' on those same compounds. Because every significant cell of the matrix is
#' fitted to the data it is then evaluated on, such estimates are optimistic,
#' and the more stringent the functional cutoff, the fewer compounds define the
#' matrix and the larger that optimism becomes.
#'
#' This function gives an out-of-sample estimate instead. The compounds are
#' split into `k` folds; for each fold the matrix is rebuilt from the other
#' folds only and used to score the held-out compounds. Scores from all folds
#' are pooled and evaluated once.
#'
#' Note that folds are drawn over compounds, which are not fully independent:
#' the panel targets only 17 genes, so compounds against the same transcript
#' share sequence context. Grouping folds by gene would be stricter still, but
#' the gene identity of each compound is not published in Supplementary Table 1.
#'
#' @param data Training data, as in [sdrna_weight_matrix()].
#' @param functional_cutoff,nonfunctional_cutoff Functionality bins, as in
#'   [sdrna_weight_matrix()].
#' @param k Number of folds.
#' @param seed Optional integer seed controlling both the fold assignment and
#'   the per-fold simulation.
#' @param ... Further arguments passed to [sdrna_weight_matrix()].
#'
#' @return A list with elements `auc` (pooled out-of-sample AUC), `cor` (pooled
#'   Pearson correlation between score and expression), `scores` (the pooled
#'   held-out scores, in the order of `data`), `fold` (fold assignment) and
#'   `k`.
#'
#' @examples
#' cv <- sdrna_cv(k = 5, seed = 1)
#' cv$auc
#' @export
sdrna_cv <- function(data = sdrnapred::sdrna356,
                     functional_cutoff = 24,
                     nonfunctional_cutoff = 44,
                     k = 10L,
                     seed = NULL,
                     ...) {
  stopifnot(is.data.frame(data), k >= 2L)
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(data)
  fold <- sample(rep_len(seq_len(k), n))
  scores <- rep(NA_real_, n)

  for (i in seq_len(k)) {
    train <- data[fold != i, , drop = FALSE]
    test <- data[fold == i, , drop = FALSE]
    n_f <- sum(train$expression < functional_cutoff)
    n_nf <- sum(train$expression > nonfunctional_cutoff)
    if (n_f < 2L || n_nf < 2L) {
      stop("Fold ", i, " leaves too few compounds in a subset; ",
           "use a smaller k or a less stringent cutoff.", call. = FALSE)
    }
    wm <- sdrna_weight_matrix(data = train,
                              functional_cutoff = functional_cutoff,
                              nonfunctional_cutoff = nonfunctional_cutoff,
                              ...)
    scores[fold == i] <- sdrna_score(test$sequence, wm)
  }

  labels <- data$expression < functional_cutoff
  list(auc = sdrna_auc(scores, labels),
       cor = stats::cor(scores, data$expression),
       scores = scores,
       fold = fold,
       k = k)
}
