#!/usr/bin/env bash

echo Downloading the files
R CMD BATCH synapse.r

echo Calculating the differential gene expression
echo BM10
R CMD BATCH DiffExpr_BM10.R
echo BM22
R CMD BATCH DiffExpr_BM22.R
echo BM36
R CMD BATCH DiffExpr_BM36.R
echo BM44
R CMD BATCH DiffExpr_BM44.R
