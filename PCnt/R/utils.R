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