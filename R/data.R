#' Training panel of 356 chemically modified self-delivering siRNAs
#'
#' The complete screening panel reported in Supplementary Table 1 of
#' Shmushkovich, Monopoli et al. (2018): 356 cholesterol-conjugated,
#' chemically modified, self-delivering siRNAs (sdRNAs) targeting 17 genes,
#' with the efficacy measured for each in a dual-luciferase reporter assay in
#' HeLa cells at 1 micromolar passive uptake.
#'
#' Each `sequence` is the 50-nucleotide mRNA target region: the 20-nucleotide
#' siRNA targeting site flanked by 15 nucleotides on each side. In the printed
#' table the targeting site is set in bold, which plain text cannot carry; it
#' is therefore characters 16 to 35 of every sequence. One row is an exception:
#' ID 1284 is printed with no bold at all in the published PDF, so its
#' targeting site is not marked in the source either. It is reproduced as
#' printed rather than inferred.
#'
#' `expression` is target gene expression remaining as a percent of
#' untreated control, so *lower* values mean *more* potent silencing. `sd` is
#' the standard deviation over three replicates. Two expression values depart
#' from the table's usual two-decimal convention and are reproduced as printed:
#' ID 1095 is 25.405 and ID 1107 is 29.2.
#'
#' The paper's functionality bins are expression below 17, 24 or 35 percent
#' (41, 91 and 138 compounds) and non-functional above 44 percent (157
#' compounds).
#'
#' @format A data frame with 356 rows and 4 columns:
#' \describe{
#'   \item{id}{Integer compound identifier, 1000 to 1355, contiguous.}
#'   \item{sequence}{Character. The 50-nucleotide target region, A/C/G/U.}
#'   \item{expression}{Numeric. Target gene expression, percent of control.}
#'   \item{sd}{Numeric. Standard deviation across three replicates.}
#' }
#'
#' @source Supplementary Table 1 of Shmushkovich T, Monopoli KR, Homsy D,
#'   Leyfer D, Betancur-Boissel M, Khvorova A, Wolfson AD (2018). Functional
#'   features defining the efficacy of cholesterol-conjugated,
#'   self-deliverable, chemically modified siRNAs. *Nucleic Acids Research*
#'   46(20):10905-10916. \doi{10.1093/nar/gky745}. Published under CC BY 4.0
#'   and redistributed here under that licence. Extracted from the publisher's
#'   PDF and independently verified cell by cell against a second parser.
"sdrna356"
