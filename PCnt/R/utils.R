require(pcalg)
#' Remove edges to satisfy dag constrains on graph
#'
#' This function does something amazing.
#'
#' @param adj_matrix adjacency matrix
#' @param topE number of top edges selected. default 150.
#'
#' @export
adj2dag = function(adj_matrix, topE=150) {
    dag = convertToDAG(adj_matrix, topE)
    rownames(dag) = rownames(adj_matrix)
    colnames(dag) = colnames(adj_matrix)
    dag
}

#' Compare two networks with same nodes embedding
#'
#' @param output_dag adjacency matrix.
#' @param truth number of top edges selected. default 150.
#' @param gvarType description.
#' @param gvarShape description.
#' @param manual_colors node colors.
#' @param topE default 200.
#' @param seed default 1.
#'
#' @export
network_compare = function(output_dag, truth, gvarType, gvarShape, manual_colors=NULL, topE = 200, seed = 1) {
  output_dag = as.matrix(output_dag)
  truth = as.matrix(truth)
  
  # remove nodes without any connections in truth graph
  zerovars = which(rowSums(abs(truth))==0 & colSums(abs(truth))==0)
  
  gvarType = gvarType[-zerovars]
  gvarShape = gvarShape[-zerovars]
  truth = truth[-zerovars, -zerovars]
  output_dag = output_dag[-zerovars, -zerovars]
  
  # plot truth graph
  truth_output = network_visualize(truth, gvarType, gvarShape, truth, topE=topE, seed=seed)
  my_coordinate = network2coordinate(truth_output$net, rownames(truth))
  
  # use truth graph layout to plot result graph
  output = network_visualize(output_dag, gvarType, gvarShape, truth, mcoord = my_coordinate, topE=topE, seed=seed)
  
  truthp = truth_output$p + scale_color_manual(values = manual_colors)
  outputp = output$p + scale_color_manual(values = manual_colors)
  
  p = ggarrange(truthp, outputp, ncol=2, nrow=1, common.legend = T)
  p
}

#' PC+NOTEARS
#'
#' @param X a numeric data matrix with samples as rows and variables as columns.
#' @param lambda1 lambda for l1 penalty, default 0.01.
#' @param losses "logistic" for binary variables and "l2" for continuous variables.
#' @param no_parent constraint that certain variables do not have parent nodes.
#' @param no_edge a 0,1 matrix where 1 species no_edge.
#' @param max.iter default 100.
#' @param h.tol default 1e-6.
#' @param rho.max default 1e+6.
#' @param w.threshold default 0.01.
#' @param alpha alpha for PC algorithm, default 0.05.
#' 
#' @return Adjacency matrix.
#'
#' @export
PC_nt = function(data, lambda1=0.01, losses=NULL, no_parents=NULL, no_edge=NULL,
                max.iter=100, h.tol=1e-6, rho.max=1e+6, w.threshold=0.01, 
                seed=1, alpha=0.05) {
  
  data = as.matrix(data)
  # Perform PC
  suffStat = list(C = cor(data), n = nrow(data))
  pc.fit = pc(suffStat,
              indepTest = gaussCItest,
              labels = colnames(data),
              alpha = alpha)
  
  pcresult = as(pc.fit@graph, "matrix")

  PC_no_edge = 1-pcresult

  # run NOTEARS with PC output as edge constraint
  ntres = notearsInterceptMultiLoss(data, no_edge = PC_no_edge, lambda1=lambda1, alpha=alpha)

  ntres$graph
}
