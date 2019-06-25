ECHO Downloading the files
R CMD BATCH synapse.r

ECHO Calculating the differential gene expression
ECHO BM10
R CMD BATCH DiffExpr_BM10.R
ECHO BM22
R CMD BATCH DiffExpr_BM22.R
ECHO BM36
R CMD BATCH DiffExpr_BM36.R
ECHO BM44
R CMD BATCH DiffExpr_BM44.R
