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

truth = read.table(file.path(file_path, "example_truth.tsv"), sep="\t")
truth = as.matrix(truth)

# take 1000 samples a small exmaple for shorter runtime
data = data[1:1000,]
```
#### Perform PC
```r
# Run PCnt
# set lambda1 for NOTEARS l1 penalty and alpha for confidence level for PC
output = PC_nt(data, lambda=0.01, alpha=0.05)

# default lambda1=0.01m alpha=0.05
output = PC_nt(data)
```

#### Enforce dag constraint on output
```r
output_dag = adj2dag(output)
```

#### Look at the metrics
```r
# F1 for orientation
getF1(output_dag, truth)
#F1 for adjacency
getAdjF1(output_dag, truth)
#Orient accuracy
getOrientAccuracy(output_dag, truth)
#SHD
myShd(output_dag, truth)
```

#### Plot output graph alone
```r
gvarType = c(rep("Causal gene",5), rep("Altered gene", 25))
gvarShape = rep("Gene", 30)
manual_colors = c("Causal gene" = "#a2d2f1", "Altered gene" = "#c8d3d5")
network_output = network_visualize(as.matrix(output_dag), gvarType, gvarShape)

network_output$p + scale_color_manual(values = manual_colors)
```
![image text](https://github.com/zhu-yh1/PC-NOTEARS/blob/main/exmaples/output_network.png)

#### Plot truth and output with same layout
```r
network_compare(output_dag, truth, gvarType, gvarShape, manual_colors, seed = 2)
```
![image text](https://github.com/zhu-yh1/PC-NOTEARS/blob/main/exmaples/network_compare.png)
