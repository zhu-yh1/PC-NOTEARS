# PC-NOTEARS
A hybrid method for causal effect estimation

### Installation
```
devtools::install_github("zhu-yh1/PC-NOTEARS/PCnt")
```

### Example

#### Load libraries
``` r
library(PCnt)
library(pcalg)
library(Rcpp)
library(RcppArmadillo)
library(ggnetwork)
library(igraph)
library(intergraph)
library(ggpubr)
```

#### Load dataset
```r
# read in example data from package
data =  read.table(file.path(file_path, "example_input.tsv"), sep="\t")
data = as.matrix(data)

# take 1000 samples a small exmaple for shorter runtime
data = data[1:1000,]

truth = read.table(file.path(file_path, "example_truth.tsv"), sep="\t")
truth = as.matrix(truth)
```
#### Perform PC
```r
suffStat = list(C = cor(data), n = nrow(data))
pc.fit = pc(suffStat,
            indepTest = gaussCItest,
            labels = colnames(data),
            alpha = 0.05)
pcres = summary(object = pc.fit)
PC_no_edge = 1-pcres
```
#### Run NOTEARS with PC output as edge constraint
```r
```
#### Enforce dag constraint on output
```r
```
#### Look at the metrics
```r
```
#### Plot output graph alone
```r
```
#### Plot truth and output with same layout
```r
```
