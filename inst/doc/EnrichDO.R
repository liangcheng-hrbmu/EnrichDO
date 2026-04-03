## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----install_chunk, eval=FALSE------------------------------------------------
# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# 
# BiocManager::install("EnrichDO")

## ----install_chunk2, eval=FALSE-----------------------------------------------
# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# 
# library(devtools)
# devtools::install_github("liangcheng-hrbmu/EnrichDO")

## ----setup,results='hide'-----------------------------------------------------
library(EnrichDO)

## -----------------------------------------------------------------------------
Alzheimer <- read.delim(file.path(
  system.file("extdata", package = "EnrichDO"),
  "Alzheimer_curated.csv"
), header = FALSE)
input_example <- Alzheimer[, 1]

## -----------------------------------------------------------------------------
demo.data <- c(1636, 351, 102, 2932, 3077, 348, 4137, 54209, 5663, 5328, 23621, 3416, 3553)

## ----eval=TRUE,results='hide',cache=TRUE,message=FALSE------------------------
weighted_demo_result <- doEnrich(interestGenes = demo.data)

## ----eval=FALSE---------------------------------------------------------------
# weighted_demo_result <- doEnrich(
#   interestGenes = demo.data,
#   test = "fisherTest",
#   method = "holm",
#   m = 1,
#   minGsize = 10,
#   maxGsize = 2000,
#   delta = 0.05,
#   penalize = TRUE
# )

## ----echo=FALSE, results='hold'-----------------------------------------------
weighted_demo_result <- doEnrich(
  interestGenes = demo.data, test = "fisherTest", method = "holm", m = 1, minGsize = 10,
  maxGsize = 2000, delta = 0.05, penalize = TRUE
)

## -----------------------------------------------------------------------------
Tradition_demo_result <- doEnrich(demo.data, traditional = TRUE)

## ----eval=FALSE---------------------------------------------------------------
# show(weighted_demo_result)

## ----echo=FALSE, results='hold'-----------------------------------------------
show(weighted_demo_result)

## ----detail-result------------------------------------------------------------
# Extract the full results
detail_res <- viewDetailResult(weighted_demo_result)

## ----summary-result-----------------------------------------------------------
# Extract the significant results
summary_res <- viewSummaryResult(weighted_demo_result)

## -----------------------------------------------------------------------------
writeResult(EnrichResult = weighted_demo_result, file = file.path(tempdir(), "result.txt"), Q = 1, P = 1)

## ----fig.cap="bar plot",fig.align='center',fig.width=8,fig.height=6-----------
drawBarGraph(EnrichResult = weighted_demo_result, n = 10, delta = 0.05)

## ----fig.cap="point plot",fig.align='center',fig.width=8,fig.height=6---------
drawPointGraph(EnrichResult = weighted_demo_result, n = 10, delta = 0.05)

## ----fig.cap="tree plot",fig.align='center',fig.width=8,fig.height=6----------
drawGraphViz(EnrichResult = weighted_demo_result, n = 10, numview = FALSE, pview = FALSE, labelfontsize = 17)

## ----fig.cap="heatmap",fig.align='center',fig.width=8,fig.height=6------------
drawHeatmap(
  interestGenes = demo.data,
  EnrichResult = weighted_demo_result,
  gene_n = 10,
  fontsize_row = 8
)

## ----fig.width=8,fig.height=6-------------------------------------------------
# Firstly, read the wrireResult output file,using the following two lines
data <- read.delim(file.path(system.file("examples", package = "EnrichDO"), "result.txt"))
enrich <- convDraw(resultDO = data)

# then, Use the drawing function you need
drawGraphViz(enrich = enrich) # Tree diagram

drawPointGraph(enrich = enrich, delta = 0.05) # Bubble diagram

drawBarGraph(enrich = enrich, delta = 0.05) # Bar plot

## ----session-info,cache = FALSE,echo=TRUE,message=TRUE,warning=FALSE----------
sessionInfo()

