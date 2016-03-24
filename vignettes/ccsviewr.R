## ---- echo=TRUE, results='asis'------------------------------------------
library(ccsviewr)
# For this example, we will load some sample data
# distributed with the package.
ccs_name = system.file("extdata", "sample.aligned.ccs.bam", package="ccsviewr")
subreads_name = system.file("extdata", "sample.subreads.bam", package="ccsviewr")
fasta_name = system.file("extdata", "sample.fna", package="ccsviewr")

# Now let's collect data for one ZMW
hole = 9 # This is the ZMW we want data for
alns = getAlignments(hole, ccs_name, subreads_name, fasta_name)

## ---- echo=TRUE----------------------------------------------------------
trimNamePrefix = function(x) sub("m160311_183010_42237_c100992522550000001823224507191630_s1_p0/9/", "", x)
scores = data.frame(name = sapply(alns, function(x) trimNamePrefix(x$id)),
                    score = sapply(alns, function(x) x$score))

knitr::kable(scores)

## ---- echo=TRUE----------------------------------------------------------
bad_alns = which(scores$score < -500)
clean_alns = alns[-bad_alns]

## ---- echo=TRUE----------------------------------------------------------
  df = AlnsToDataFrame(clean_alns)


## ---- echo=TRUE, eval=FALSE----------------------------------------------
#  plotMSA(df, "FullAlignment.pdf")

## ---- echo=TRUE, eval=FALSE----------------------------------------------
#  # Clean up the read names
#  rnames = as.character(df$id)
#  df$id = factor(sapply(rnames, trimNamePrefix))
#  # Now make a cleaner looking PDF
#  plotMSA(df, "AlignmentWindow.pdf", start = 40, end=60, showPositions = FALSE)

