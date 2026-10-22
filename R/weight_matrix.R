#' Bases used throughout the package, in matrix row order
#' @keywords internal
SDRNA_BASES <- c("A", "C", "G", "U")

#' Split sequences into a character matrix
#'
#' @param sequences Character vector of equal-length RNA sequences.
#' @return A character matrix with one row per sequence and one column per
#'   position.
#' @keywords internal
sdrna_split <- function(sequences) {
  sequences <- toupper(gsub("T", "U", sequences, fixed = TRUE))
  widths <- nchar(sequences)
  if (length(unique(widths)) != 1L) {
    stop("All sequences must have the same length; found lengths ",
         paste(sort(unique(widths)), collapse = ", "), ".", call. = FALSE)
  }
  m <- do.call(rbind, strsplit(sequences, "", fixed = TRUE))
  bad <- setdiff(unique(as.vector(m)), SDRNA_BASES)
  if (length(bad)) {
    stop("Sequences contain characters outside A/C/G/U: ",
         paste(bad, collapse = ", "), ".", call. = FALSE)
  }
  m
}

#' Per-position base counts
#'
#' @param m Character matrix from [sdrna_split()].
#' @return A 4 x L integer matrix of counts, rows named A/C/G/U.
#' @keywords internal
sdrna_counts <- function(m) {
  L <- ncol(m)
  out <- matrix(0L, nrow = 4L, ncol = L,
                dimnames = list(SDRNA_BASES, seq_len(L)))
  for (j in seq_len(L)) {
    tab <- table(factor(m[, j], levels = SDRNA_BASES))
    out[, j] <- as.integer(tab)
  }
  out
}

#' Build a positional base-preference weight matrix
#'
#' Reimplements the weight-matrix construction of Shmushkovich, Monopoli et al.
#' (2018), section "Generating the weight matrix".
#'
#' The training compounds are split into a functional subset (target expression
#' strictly below `functional_cutoff`, as a percent of control) and a
#' non-functional subset (target expression strictly above
#' `nonfunctional_cutoff`). For every base at every retained position the
#' weight is the functional per-position base frequency minus the
#' non-functional frequency, but only where that difference is significant;
#' non-significant cells are set to zero.
#'
#' Significance follows the paper's Equation 2. Two collections of random
#' sequence sets are simulated, one sized like the functional subset (`q`
#' sequences) and one like the non-functional subset (`z` sequences), both
#' drawn position by position from the base composition of the whole training
#' set. `n_replicates` sets of each size are generated, giving a simulated
#' spread of per-position base counts, whose standard deviations enter the
#' denominator of the test statistic. A two-sided one-sample t-test with
#' Bonferroni correction across all tested cells yields the p-values.
#'
#' @param data Training data frame with columns `sequence` (equal-length RNA
#'   strings) and `expression` (target gene expression as a percent of
#'   control). Defaults to the bundled [sdrna356] dataset.
#' @param functional_cutoff Numeric. Compounds below this expression value are
#'   functional. The paper considers 17, 24 and 35; it selected 24.
#' @param nonfunctional_cutoff Numeric. Compounds above this value are
#'   non-functional. 44 in the paper.
#' @param positions Integer vector of sequence positions to retain. The default
#'   `2:49` reproduces the paper's primary algorithm, which excludes positions
#'   1 and 50 to avoid artefacts from the cloning of multiple target sites. Use
#'   `16:35` for the 20-nucleotide targeting region alone, as in the paper's
#'   Figure 4.
#' @param alpha Significance threshold applied after Bonferroni correction.
#'   0.001 in the paper.
#' @param n_replicates Number of simulated random sets per subset size. The
#'   paper generates `N = 356` of each.
#' @param background Optional 4 x L matrix of per-position base probabilities
#'   used to simulate the random sets, rows ordered A/C/G/U. Defaults to the
#'   composition of `data` itself. Supplying an external background (for
#'   example RefSeq frequencies) reproduces the paper's control algorithm; see
#'   [sdrna_control_matrix()].
#' @param reference Optional 4 x L matrix of per-position base probabilities
#'   that *replaces the non-functional subset* on the comparison side of the
#'   weight calculation, rows ordered A/C/G/U. This is how the paper builds its
#'   control algorithm; see [sdrna_control_matrix()]. When `NULL` (the
#'   default) the non-functional subset of `data` is used, which is the main
#'   algorithm.
#' @param scale Multiplier applied to the returned weights. The default 100
#'   expresses them in percentage points, matching the integer values printed
#'   in the paper's Figure 2A and Figure 4A.
#' @param seed Optional integer seed for the simulation, for reproducibility.
#'
#' @return An object of class `sdrna_wm`: a 4 x L numeric matrix of weights
#'   (zero at positions not in `positions`), carrying attributes `positions`,
#'   `p_values`, `significant`, `n_functional`, `n_nonfunctional`, `cutoffs`
#'   and `scale`.
#'
#' @references Shmushkovich T, Monopoli KR, Homsy D, Leyfer D,
#'   Betancur-Boissel M, Khvorova A, Wolfson AD (2018). Functional features
#'   defining the efficacy of cholesterol-conjugated, self-deliverable,
#'   chemically modified siRNAs. *Nucleic Acids Research* 46(20):10905-10916.
#'
#' @examples
#' wm <- sdrna_weight_matrix(functional_cutoff = 24, seed = 1)
#' round(wm[, 16:20])
#' @export
sdrna_weight_matrix <- function(data = sdrnapred::sdrna356,
                                functional_cutoff = 24,
                                nonfunctional_cutoff = 44,
                                positions = 2:49,
                                alpha = 0.001,
                                n_replicates = 356L,
                                background = NULL,
                                reference = NULL,
                                scale = 100,
                                seed = NULL) {
  stopifnot(is.data.frame(data),
            all(c("sequence", "expression") %in% names(data)))
  if (functional_cutoff >= nonfunctional_cutoff) {
    stop("functional_cutoff must be below nonfunctional_cutoff.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)

  m <- sdrna_split(data$sequence)
  L <- ncol(m)
  positions <- sort(unique(as.integer(positions)))
  if (any(positions < 1L | positions > L)) {
    stop("positions must lie within 1:", L, ".", call. = FALSE)
  }

  expr <- as.numeric(data$expression)
  f <- which(expr < functional_cutoff)
  nf <- which(expr > nonfunctional_cutoff)
  q <- length(f)
  z <- length(nf)
  if (q < 2L || z < 2L) {
    stop("Both subsets need at least two compounds; got q = ", q,
         " and z = ", z, ".", call. = FALSE)
  }

  count_f <- sdrna_counts(m[f, , drop = FALSE])
  if (is.null(reference)) {
    count_nf <- sdrna_counts(m[nf, , drop = FALSE])
  } else {
    stopifnot(is.matrix(reference), nrow(reference) == 4L,
              ncol(reference) == L)
    reference <- sweep(reference, 2L, colSums(reference), "/")
    # The reference composition stands in for the non-functional subset,
    # expressed as the counts it would produce in a subset of that size.
    count_nf <- reference * z
  }

  if (is.null(background)) {
    background <- sdrna_counts(m)
    background <- sweep(background, 2L, colSums(background), "/")
  } else {
    stopifnot(is.matrix(background), nrow(background) == 4L,
              ncol(background) == L)
    background <- sweep(background, 2L, colSums(background), "/")
  }

  # Simulated spread of per-position base counts at each subset size.
  sd_f <- matrix(NA_real_, 4L, L, dimnames = dimnames(count_f))
  sd_nf <- sd_f
  for (j in seq_len(L)) {
    p <- background[, j]
    sd_f[, j] <- apply(stats::rmultinom(n_replicates, q, p), 1L, stats::sd)
    sd_nf[, j] <- apply(stats::rmultinom(n_replicates, z, p), 1L, stats::sd)
  }

  # Equation 2.
  numerator <- count_nf / z - count_f / q
  denominator <- sqrt((sd_nf^2 / z^2) / z + (sd_f^2 / q^2) / q)
  tstat <- numerator / denominator
  tstat[!is.finite(tstat)] <- 0

  df <- n_replicates - 1L
  praw <- 2 * stats::pt(-abs(tstat), df = df)
  n_tests <- 4L * length(positions)
  padj <- pmin(praw * n_tests, 1)
  padj[, setdiff(seq_len(L), positions)] <- NA_real_

  sig <- !is.na(padj) & padj < alpha
  weights <- (count_f / q - count_nf / z) * scale
  weights[!sig] <- 0

  structure(weights,
            class = c("sdrna_wm", "matrix", "array"),
            positions = positions,
            p_values = padj,
            significant = sig,
            n_functional = q,
            n_nonfunctional = z,
            cutoffs = c(functional = functional_cutoff,
                        nonfunctional = nonfunctional_cutoff),
            scale = scale)
}

#' Build the paper's control weight matrix
#'
#' The paper derives a control algorithm by replacing the training set's own
#' per-position base composition with average RefSeq database frequencies,
#' adjusted for the initial bias of the training set, and otherwise repeating
#' the construction unchanged. Supply those frequencies as `refseq`; when they
#' are omitted a uniform 25 percent composition is used, which is a neutral
#' stand-in rather than the paper's own control.
#'
#' @param data Training data, as in [sdrna_weight_matrix()].
#' @param refseq Optional numeric vector of four base probabilities (A, C, G,
#'   U) or a 4 x L matrix of per-position probabilities.
#' @param ... Further arguments passed to [sdrna_weight_matrix()].
#'
#' @return An `sdrna_wm` object, as returned by [sdrna_weight_matrix()].
#' @export
sdrna_control_matrix <- function(data = sdrnapred::sdrna356,
                                 refseq = NULL, ...) {
  m <- sdrna_split(data$sequence)
  L <- ncol(m)
  if (is.null(refseq)) {
    bg <- matrix(0.25, 4L, L, dimnames = list(SDRNA_BASES, seq_len(L)))
  } else if (is.matrix(refseq)) {
    bg <- refseq
  } else {
    stopifnot(length(refseq) == 4L)
    bg <- matrix(refseq, nrow = 4L, ncol = L,
                 dimnames = list(SDRNA_BASES, seq_len(L)))
  }
  sdrna_weight_matrix(data = data, background = bg, reference = bg, ...)
}

#' @export
print.sdrna_wm <- function(x, ...) {
  cut <- attr(x, "cutoffs")
  cat("sdRNA positional base-preference weight matrix\n")
  cat("  functional < ", cut[["functional"]], "% expression (n = ",
      attr(x, "n_functional"), ")\n", sep = "")
  cat("  non-functional > ", cut[["nonfunctional"]], "% expression (n = ",
      attr(x, "n_nonfunctional"), ")\n", sep = "")
  cat("  positions retained: ", min(attr(x, "positions")), "-",
      max(attr(x, "positions")), " (", length(attr(x, "positions")),
      " of ", ncol(x), ")\n", sep = "")
  cat("  non-zero weights: ", sum(x != 0), " of ",
      4L * length(attr(x, "positions")), " tested cells\n", sep = "")
  invisible(x)
}
