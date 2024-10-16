# PC-NOTEARS
A hybrid method for causal effect estimation

### Paper
[A hybrid constrained continuous optimization approach for optimal causal discovery from biological data](https://academic.oup.com/bioinformatics/article/40/Supplement_2/ii87/7749067)
```
@article{Zhu2024Sep,
	author = {Zhu, Yuehua and Benos, Panayiotis V. and Chikina, Maria},
	title = {{A hybrid constrained continuous optimization approach for optimal causal discovery from biological data}},
	journal = {Bioinformatics},
	volume = {40},
	number = {Supplement_2},
	pages = {ii87--ii97},
	year = {2024},
	month = sep,
	issn = {1367-4811},
	publisher = {Oxford Academic},
	doi = {10.1093/bioinformatics/btae411}
}
```

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
file_path <- system.file(package = "PCnt")
print(file_path)

# input data
data =  read.table(file.path(file_path, "example_input.tsv"), sep="\t")
data = as.matrix(data)

# matched ground truth
truth = read.table(file.path(file_path, "example_truth.tsv"), sep="\t")
truth = as.matrix(truth)

# take 1000 samples a small exmaple for shorter runtime
data = data[1:1000,]
```
#### Perform PCnt
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
# in our example ground truth causal graph there are 5 causal genes and 25 altered genes. 
gvarType = c(rep("Causal gene",5), rep("Altered gene", 25))
# This is a subgroup
gvarShape = rep("Gene", 30)
# set manual colors for nodes
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
