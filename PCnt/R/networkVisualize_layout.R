library(ggnetwork)
library(igraph)
library(intergraph)
library(ggpubr)

getLabel = function(graph, truth) {
	n = ncol(truth)
	m = ncol(graph)
	res = matrix(0, nrow=m, ncol=m)
	rownames(res) = rownames(graph)
	colnames(res) = rownames(graph)
	for (i in 1:n) {
		for (j in 1:n) {
			if (graph[i,j] != 0) {
				if (truth[i,j] != 0) {
					res[i,j] = "truePositive"
				} else if (truth[j,i] != 0) {
					res[i,j] = "wrongDirection"
				} else {
					res[i,j] = "falsePositive"
				}
			}
		}
	}
	res
}

network2coordinate <- function(network, gene_list) {
	coord_df = network[c("x","y","vertex.names")]
	coord_df = unique(coord_df)
	rownames(coord_df) = coord_df$vertex.names
	coord_df = coord_df[gene_list,]
	coord_matrix = as.matrix(coord_df[c("x","y")])
}

#' Visualize network based on adjacency matrix
#'
#' network visualization function.
#'
#' @param graph adjacency matrix for graph.
#' @param gvarType variable group1. annotated in different colors.
#' @param gvarShape variable group2. annotated in different shapes.
#' @param truth truth matrix for edge labeling.
#' @param mcoord coordinates for plot embedding.
#' @param topE number of top edges to be plotted.
#' 
#' @return network graph and node coordinates.
#'
#' @export
network_visualize <- function(graph, gvarType, gvarShape, truth = NULL, mcoord = NULL, topE=50, layout.name = "fruchtermanreingold", alpha.range = c(0.05,1), edge.weight = T, use.alpha=F, sym = F, label.color="white", seed=42, nodeLabel=T){
	set.seed(seed)
	if (sym) {
		for (i in 1:nrow(graph))
			for (j in 1:i) {
				graph[i,j] = 0
			}
	}

	edgeLabel = NULL

	if (!is.null(truth)) {
		edgeLabel=getLabel(graph, truth)
	}
	
	cutoff=sort(abs(as.vector(graph)),T)[topE+1]
	
	graph[abs(graph)<cutoff]=0
	zero.vars=which(rowSums(abs(graph))==0&colSums(abs(graph))==0)
	if (!is.null(edgeLabel))
        edgeLabel[which(abs(graph)<cutoff)]=0

	# if(length(zero.vars)>0){
	# 	graph=graph[-zero.vars, -zero.vars]
	# 	gvarType=gvarType[-zero.vars]
	# 	gvarShape=gvarShape[-zero.vars]
	# }  

	#convert the graph to data.frame
	#######################################################################
	graph.binary <- matrix(0, nrow = nrow(graph), ncol = ncol(graph))
	graph.binary[graph!=0] <- 1
	rownames(graph.binary) <- rownames(graph)
	colnames(graph.binary) <- colnames(graph)
	net <- graph_from_adjacency_matrix(graph.binary, mode = "directed", diag = F)
	
	if (!is.null(mcoord)) {
		net.dat <- ggnetwork(asNetwork(net), layout = mcoord, params = list(exp = 0))
	} else {
		net.dat <- ggnetwork(asNetwork(net), layout = layout.name, params = list(exp = 0))
	}
	
	
	
	#add gvartype, node color, node shape
	#######################################################################
	num_of_edges <- sum(graph.binary==1)
	num_of_nodes <- nrow(graph)
	
	#add var
	var <- rep(NA, nrow(net.dat))
	var[(num_of_edges+1):nrow(net.dat)] <- gvarType
	net.dat$var <- var
	
	shape <- rep(NA, nrow(net.dat))
	shape[(num_of_edges+1):nrow(net.dat)] <- gvarShape
	net.dat$shape <- shape
	
	
	#add edge weight option
	#using [0,1] to rescale the weights
	#######################################################################
	weights <- rep(NA, nrow(net.dat))
    edgeOpt <- rep(NA, nrow(net.dat))
	node.ref <- net.dat[(num_of_edges+1):nrow(net.dat),]
	graph.verify <- matrix(0, nrow = nrow(graph), ncol = ncol(graph))
	rownames(graph.verify) <- rownames(graph)
	colnames(graph.verify) <- colnames(graph)
	
	for (i in 1:num_of_edges){
		row.name <- as.character(node.ref$vertex.names)[which(node.ref[,1]==net.dat[i,1]& node.ref[,2]== net.dat[i,2])]
		col.name <- as.character(node.ref$vertex.names)[which.min(abs(net.dat[i,4]-node.ref[,1])+abs(net.dat[i,5]-node.ref[,2]))]
		
		graph.verify[row.name, col.name] <- graph[row.name, col.name]
		#go back to graph and find weights
		weights[i] <-   graph[row.name, col.name]
		if(!is.null(edgeLabel))
			edgeOpt[i] <-   edgeLabel[row.name, col.name]
		if (graph[row.name, col.name]==0){
			print("warning")
			print(i)
		}#if
	}#for i
	
	#verify the parsing
	#all(graph==graph.verify)
	
	
	net.dat$weights <- weights
	net.dat$edgeOpt <- edgeOpt

	#may add gap to make the edge look better
	#gap already exists
	proj <- function(x, x.min = alpha.range[1], x.max = alpha.range[2]){
		if(max(x)==min(x)) {
			x.frac=1
		}
		else {
			x.frac <- (x-min(x))/(max(x)-min(x))
		}
		res <- x.frac*(x.max-x.min)+x.min
		return(res)
	}#proj
	
	
	weights.inuse <- weights
	weights.inuse[!is.na(weights.inuse)] <- abs(weights.inuse[!is.na(weights.inuse)])
	weights.inuse[!is.na(weights.inuse)] <- proj(weights.inuse[!is.na(weights.inuse)])
	
	net.dat$weights.inuse <- weights.inuse
	#https://www.andrea-rau.com/post/alpha-transparency-in-ggplot2/, doesn't work
	#net.dat <- net.dat%>% mutate(weights.inuse = I(weights.inuse))
	
	
	#fix edge weights
	######################
	if (!edge.weight){
		if (!sym) {
			p <- ggplot(net.dat, aes(x=x, y = y, xend = xend, yend = yend), arrow.gap = 1)+
				theme_blank()+
				geom_edges( curvature = 0, size = 0.3)+
				geom_edges(linetype = "blank",arrow = arrow(length = unit(5, "pt"), type = "closed"), curvature = 0)+
				geom_nodelabel_repel(aes(label = vertex.names, color = var, alpha=0),box.padding = unit(1, "lines"), size = 4)+
				labs(shape="", color ="", linetype="")
		}
		else {
			p <- ggplot(net.dat, aes(x=x, y = y, xend = xend, yend = yend), arrow.gap = 1)+
				theme_blank()+
				geom_edges( curvature = 0, size = 0.3)+
				geom_edges(linetype = "blank",arrow = arrow(length = unit(5, "pt"), type = "closed"), curvature = 0)+
				geom_nodelabel_repel(aes(label = vertex.names, color = var, alpha=0),box.padding = unit(1, "lines"), size = 4)+
				labs(shape="", color ="", linetype="")
		}
		
	}
	else{
		#proportion weights
		######################
		p <- ggplot(net.dat, aes(x =x, y = y, xend = xend, yend = yend), arrow.gap = 1) +theme_blank()
		
		for (i in 1:num_of_edges){
			#p <- p+geom_edges(linetype = "blank",arrow = arrow(length = unit(7, "pt"), type = "closed"), curvature = 0, alpha = net.dat[i,"weights.inuse"])
			if(use.alpha){
				if (!sym) {
					p <- p+geom_edges(data = net.dat[i,,drop = F], aes(x=x, y = y, xend = xend, yend = yend), linetype = "blank",arrow = arrow(length = unit(9, "pt"), type = "closed"), alpha = net.dat[i,"weights.inuse"], curvature = 0) +
						geom_nodelabel_repel(data = net.dat[(num_of_edges+1):nrow(net.dat),],
											 aes(x = x, y = y, label = vertex.names, fill = var, alpha=0),
											 box.padding = unit(1, "lines"),
											 size = 5, color = "white")+
						
						guides(fill = guide_legend(title = "Type", override.aes = aes(label="")))
				}
				else {
					p <- p+geom_edges(data = net.dat[i,,drop = F], aes(x=x, y = y, xend = xend, yend = yend), linetype = "blank", alpha = net.dat[i,"weights.inuse"], curvature = 0) +
						geom_nodelabel_repel(data = net.dat[(num_of_edges+1):nrow(net.dat),],
											 aes(x = x, y = y, label = vertex.names, fill = var, alpha=0),
											 box.padding = unit(1, "lines"),
											 size = 5, color = "white")+
						guides(fill = guide_legend(title = "Type", override.aes = aes(label="")))
				}
				if (net.dat[i,"weights"] < 0){
					p <- p+geom_edges(data = net.dat[i,,drop = F], aes(x=x, y = y, xend = xend, yend = yend), alpha = net.dat[i,"weights.inuse"], curvature = 0, linetype = "dashed")  
				}
				else{
					p <- p+geom_edges(data = net.dat[i,,drop = F], aes(x=x, y = y, xend = xend, yend = yend), alpha = net.dat[i,"weights.inuse"], curvature = 0, linetype = "solid")  
				}#else
			}
			else{
				if (!sym) {
                    if (!is.null(edgeLabel)) {
                        if (net.dat[i,"edgeOpt"] == "falsePositive") {
                            p <- p+geom_edges(data = net.dat[i,,drop = F], aes(x=x, y = y, xend = xend, yend = yend), size = net.dat[i,"weights.inuse"]+0.3, curvature = 0, color="grey",arrow = arrow(length = unit(9, "pt"), type = "closed"))  
                        } else if (net.dat[i,"edgeOpt"] == "truePositive") {
                            p <- p+geom_edges(data = net.dat[i,,drop = F], aes(x=x, y = y, xend = xend, yend = yend), size = net.dat[i,"weights.inuse"]+0.3, curvature = 0, color="#009300",arrow = arrow(length = unit(9, "pt"), type = "closed"))
                        } else if (net.dat[i,"edgeOpt"] == "wrongDirection") {
                            p <- p+geom_edges(data = net.dat[i,,drop = F], aes(x=x, y = y, xend = xend, yend = yend), size = net.dat[i,"weights.inuse"]+0.3, curvature = 0, color="#e10000",arrow = arrow(length = unit(9, "pt"), type = "closed"))
                        }
                    } else {
                        if (net.dat[i,"weights"] > 0){
                            p <- p+geom_edges(data = net.dat[i,,drop = F], aes(x=x, y = y, xend = xend, yend = yend), size = net.dat[i,"weights.inuse"]+0.3, curvature = 0, color="black",arrow = arrow(length = unit(9, "pt"), type = "closed"))  
                        }
                        else{
                            p <- p+geom_edges(data = net.dat[i,,drop = F], aes(x=x, y = y, xend = xend, yend = yend), size = net.dat[i,"weights.inuse"]+0.3, curvature = 0, color="grey",arrow = arrow(length = unit(9, "pt"), type = "closed"))  
                        }#else
                    }
				}
				else {
					if (!is.null(edgeLabel)) {
                        if (net.dat[i,"edgeOpt"] == "falsePositive") {
                            p <- p+geom_edges(data = net.dat[i,,drop = F], aes(x=x, y = y, xend = xend, yend = yend), size = net.dat[i,"weights.inuse"]+0.3, curvature = 0, color="grey")  
                        } else if (net.dat[i,"edgeOpt"] == "truePositive") {
                            p <- p+geom_edges(data = net.dat[i,,drop = F], aes(x=x, y = y, xend = xend, yend = yend), size = net.dat[i,"weights.inuse"]+0.3, curvature = 0, color="#009300")
                        } else if (net.dat[i,"edgeOpt"] == "wrongDirection") {
                            p <- p+geom_edges(data = net.dat[i,,drop = F], aes(x=x, y = y, xend = xend, yend = yend), size = net.dat[i,"weights.inuse"]+0.3, curvature = 0, color="#009300")
                        }
                    }
					else {
						if (net.dat[i,"weights"] > 0){
							p <- p+geom_edges(data = net.dat[i,,drop = F], aes(x=x, y = y, xend = xend, yend = yend), size = net.dat[i,"weights.inuse"]+0.3, curvature = 0, color="black")
						}
						else{
							p <- p+geom_edges(data = net.dat[i,,drop = F], aes(x=x, y = y, xend = xend, yend = yend), size = net.dat[i,"weights.inuse"]+0.3, curvature = 0, color="grey")  
						}#else
					}
					
				}
			}
		}#for i
		
		if (nodeLabel) {
			p <- p + geom_nodes(data = net.dat[(num_of_edges+1):nrow(net.dat),], aes(x=x, y=y, shape = shape, color = var), size = 10) +
			geom_nodetext(data = net.dat[(num_of_edges+1):nrow(net.dat),],
						  aes(x = x, y = y, label = vertex.names),
						  size = 4, color = "black", fontface = "bold")+
			xlim(c(-0.1,1.2)) +
			guides(shape = guide_legend(title = "Feature Type", override.aes = aes(label="", size = 5),)) +
			guides(color = guide_legend(title = "Feature Subtype", override.aes = aes(label="", size = 5))) +
			theme(legend.title =element_text(size=12),  legend.text=element_text(size=10))
		}
		else {
			p <- p + geom_nodes(data = net.dat[(num_of_edges+1):nrow(net.dat),], aes(x=x, y=y, shape = shape, color = var), size = 10) +
			xlim(c(-0.1,1.2)) +
			guides(shape = guide_legend(title = "Feature Type", override.aes = aes(label="", size = 5),)) +
			guides(color = guide_legend(title = "Feature Subtype", override.aes = aes(label="", size = 5))) +
			theme(legend.title =element_text(size=12),  legend.text=element_text(size=10))
		}

	}#else
	
	# return(p)
	return(list("p" = p, "net" = net.dat))
}#network_visualize

