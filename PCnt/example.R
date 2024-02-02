library(PCnt)
library(pcalg)
library(Rcpp)
library(RcppArmadillo)
library(ggnetwork)
library(igraph)
library(intergraph)
library(ggpubr)

file_path <- system.file(package = "PCnt")
print(file_path)

# read in example data from package
data =  read.table(file.path(file_path, "example_input.tsv"), sep="\t")
data = as.matrix(data)

# take 1000 samples a small exmaple for shorter runtime
data = data[1:1000,]

truth = read.table(file.path(file_path, "example_truth.tsv"), sep="\t")
truth = as.matrix(truth)

# Perform PC
suffStat = list(C = cor(data), n = nrow(data))
pc.fit = pc(suffStat,
            indepTest = gaussCItest,
            labels = colnames(data),
            alpha = 0.05)
pcres = summary(object = pc.fit)
PC_no_edge = 1-pcres

# run NOTEARS with PC output as edge constraint
ntres = notearsInterceptMultiLoss(data, lambda1 = 0.01, no_edge = PC_no_edge)

# get DAG
output_dag = adj2dag(ntres$graph)

# look at the metrics
# F1 for orientation
getF1(output_dag, truth)
#F1 for adjacency
getAdjF1(output_dag, truth)
#Orient accuracy
getOrientAccuracy(output_dag, truth)
#SHD
myShd(output_dag, truth)


# plot output graph alone
gvarType = c(rep("Causal gene",5), rep("Altered gene", 25))
gvarShape = rep("Gene", 30)
manual_colors = c("Causal gene" = "#a2d2f1", "Altered gene" = "#c8d3d5")
network_output = network_visualize(as.matrix(output_dag), gvarType, gvarShape)

pdf("test/output_network.pdf", width=10,height=8)
network_output$p + scale_color_manual(values = manual_colors)
dev.off()

# plot truth and output with truth layout

pdf("test/network_compare.pdf", width=20,height=8)
network_compare(output_dag, truth, gvarType, gvarShape, manual_colors, seed = 2)
dev.off()

