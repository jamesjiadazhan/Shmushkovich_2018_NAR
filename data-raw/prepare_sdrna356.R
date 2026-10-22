## Builds data/sdrna356.rda from the verified Paper2Agent reading package.
## Source: Supplementary Table 1 of doi:10.1093/nar/gky745 (CC BY 4.0),
## delivered as nine CSV parts by the paper skill.
skill <- "../monopoli-2018-hsirna-efficacy-paper"
parts <- sprintf("%s/assets/supp_table/supplementary-table-1-part-%d.csv", skill, 1:9)
stopifnot(all(file.exists(parts)))

raw <- do.call(rbind, lapply(parts, read.csv, check.names = FALSE,
                             colClasses = "character"))
stopifnot(nrow(raw) == 356L)

sdrna356 <- data.frame(
  id         = as.integer(raw[[1]]),
  sequence   = toupper(trimws(raw[[2]])),
  expression = as.numeric(raw[[3]]),
  sd         = as.numeric(raw[[4]]),
  stringsAsFactors = FALSE
)

## Invariants asserted by the source review, re-checked here.
stopifnot(
  identical(sdrna356$id, 1000:1355),
  all(nchar(sdrna356$sequence) == 50L),
  !any(grepl("[^ACGU]", sdrna356$sequence)),
  all(is.finite(sdrna356$expression)), all(is.finite(sdrna356$sd))
)

save(sdrna356, file = "data/sdrna356.rda", compress = "xz", version = 2)
cat("wrote", nrow(sdrna356), "rows\n")
