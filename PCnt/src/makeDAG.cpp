#include <RcppArmadillo.h>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <bits/stdc++.h>

// Using the Rcpp namespace
// using namespace Rcpp;

// Constants
const double oo = 1e18;
const double eps = 1e-10;
const int MAXN = 1e3, MAXM = 1e6, MAXC = 1e6;

// Edge structure
struct Edge {
	int from, to, next;
	double weight;
	
	Edge() {}
	Edge(int a, int b, double v) {
		from = a;
		to = b;
		weight = v;
	}
	
	void clear () {
		from = 0;
		to = 0;
		weight = 0;
	}
	
	void print() {
		Rcpp::Rcout << "Edge(" << from << "," << to << ") = " << weight << std::endl;
	}
};
Edge edges[MAXM];

// Cycle structure
struct Cycle {
	int len;
	std::vector<int> node_index;
	std::vector<int> edge_index;
	
	Cycle () {}
	
	void clear() {
		len = 0;
		node_index.clear();
		edge_index.clear();
	}
	
	void print() {
		Rcpp::Rcout << "new cycle" << std::endl << "cycle len: " << len << std::endl;
		Rcpp::Rcout << "Nodes:";
		for (int i = 0; i < node_index.size(); ++i)
			Rcpp::Rcout << node_index[i] << " ";
		Rcpp::Rcout << std::endl;
		Rcpp::Rcout << "Edges:" << std::endl;
		for (int i = 0; i < edge_index.size(); ++i)
			Rcpp::Rcout << "edges["<< edge_index[i] << "] = (" << edges[edge_index[i]].from << ", " << edges[edge_index[i]].to << ")" << std::endl;
		Rcpp::Rcout << std::endl << std::endl;
	}
};
Cycle cycles[MAXC];


// Global variables
int n = 0, cnt_e = 0, cnt_c = 0, cnt_scc = 0, tot = 0, first[MAXN] = {-1},  dfn[MAXN] = {-1}, low[MAXN] = {-1}, color[MAXN] = {0};
double adj_matrix[MAXN][MAXN] = {0};
bool inStack[MAXN] = {0}, visited[MAXN] = {0}, blocked[MAXN] = {0}, removed[MAXN][MAXN] = {0};
std::vector <int> B[MAXN];

std::stack <int> s;
std::stack <int> scc[MAXN]; // It is ok not to save this information
void tarjan(int now) {
	dfn[now] = tot;
	low[now] = tot;
	tot++;
	//	Rcpp::Rcout << "tarjan:" << now << " " << tot << " " << dfn[now] << " " << low[now] << std::endl;
	s.push(now);
	inStack[now] = 1;
	for (int e = first[now]; e != -1; e = edges[e].next) {
		int to = edges[e].to;
		if (dfn[to] == -1) {
			tarjan(to);
			low[now] = std::min(low[now], low[to]);
		}
		else if (inStack[to])
			low[now] = std::min(low[now], low[to]);
	}
	
	// save scc with > 1 ndoes in stack
	if (dfn[now] == low[now]) {
		int tmp = s.top();
		s.pop();
		inStack[tmp] = 0;
		
		if (tmp == now) return;
		
		scc[cnt_scc].push(tmp);
		color[tmp] = cnt_scc;
		
		while (tmp != now) {
			tmp = s.top();
			inStack[tmp] = 0;
			s.pop();
			
			scc[cnt_scc].push(tmp);
			color[tmp] = cnt_scc;
		}
		cnt_scc++;
	}
	return;
}

// unblock now-related nodes
void unblock(int now) {
	blocked[now] = 0;
	for (int i = 0; i < B[now].size(); ++i)
		if (blocked[B[now][i]])
			unblock(B[now][i]);
		B[now].clear();
}

// johnson
std::vector <int> johnNodes;
std::vector <int> johnEdges;
bool johnsonFindCycles(int start, int now) {
	bool find = 0;
	johnNodes.push_back(now);
	blocked[now] = 1;
	//	Rcpp::Rcout << start << " " << now << std::endl;
	for (int e = first[now]; e != -1; e = edges[e].next) {
		int to = edges[e].to;
		// continue if node 'to' is not in the same scc as node 'i'
		if (color[now] != color[to])
			continue;
		if (!visited[to]) {
			if (to == start) {
				// save new cycle
				cycles[cnt_c].node_index = johnNodes;
				cycles[cnt_c].edge_index = johnEdges;
				cycles[cnt_c].edge_index.push_back(e);
				cycles[cnt_c].len = johnNodes.size();
				cnt_c++;
				find = 1;
				
				if (cnt_c >= MAXC) {
					Rcpp::stop(" More than 1,000,000 cycles found. Please reduce topE.");
				}
			}
			else if (!blocked[to]) {
				johnEdges.push_back(e);
				find = johnsonFindCycles(start, to) || find;
				johnEdges.pop_back();
			}
		}
	}
	
	if (find)	unblock(now);
	else {
		for (int e = first[now]; e != -1; e = edges[e].next) {
			int to = edges[e].to;
			if (color[now] != color[to])
				continue;
			if (!visited[to])
				B[to].push_back(now);
		}
	}
	johnNodes.pop_back();
	return find;
}

// Remove all cycles
void removeAllCycles() {
	for (int i = 0; i < cnt_c; ++i) {
		double minE = oo;
		int min_index = -1;
		
		for (int j = 0; j < cycles[i].len; ++j) {
			int e = cycles[i].edge_index[j];
			int to = edges[e].to, from = edges[e].from;
			
			if (removed[from][to]) {
				min_index = -1;
				break;
			}
			
			if (edges[e].weight < minE) {
				minE = edges[e].weight;
				min_index = e;
			}
		}
		
		if (min_index != -1) {
			removed[edges[min_index].from][edges[min_index].to] = 1;
		}
	}
}

// Compare two cycles
bool cmp(Cycle x, Cycle y) {
	return x.len <= y.len;
}

void init() {
	// initialization
	cnt_c = 0;
	cnt_e = 0;
	cnt_scc = 0;
	memset(first, -1, sizeof(first));
	memset(removed, 0, sizeof(removed));
	
	//	tarjan arrays
	tot = 0;
	memset(inStack, 0, sizeof(inStack));
	memset(dfn, -1, sizeof(dfn));
	memset(low, -1, sizeof(low));
	memset(color, -1, sizeof(color));
	
	// johnson arrays
	memset(visited, 0, sizeof(visited));
}

void johnsonClear() {
	memset(blocked, 0, sizeof(blocked));
	johnNodes.clear();
	johnEdges.clear();
	for (int i = 0; i < n; ++i)
		B[i].clear();
}

void preprocess() {
	// remove self-cycle
	for (int i = 0; i < n; ++i)
		adj_matrix[i][i] = 0;
	
	// remove A-B-A cycles
	for (int i = 0; i < n; ++i)
		for (int j = i+1; j < n; ++j) {
			if (std::abs(adj_matrix[i][j]) >= std::abs(adj_matrix[j][i]))	adj_matrix[j][i] = 0;
			else	adj_matrix[i][j] = 0;
		}
}

double sortedEdges[MAXM] = {0};
void getTopE(int topE) {
	int cnt = 0;
	double thres = 0;
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
		{
			if(std::abs(adj_matrix[i][j]) > eps)
				sortedEdges[cnt++] = std::abs(adj_matrix[i][j]);
		}
		if (cnt < topE) {
			thres = 0;
			return;
		}
		
		std::sort(sortedEdges, sortedEdges + cnt);
		
		thres = sortedEdges[cnt - topE];
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j) {
				if (std::abs(adj_matrix[i][j]) < thres)
					adj_matrix[i][j] = 0;
			}
}

void clear() {
	for (int i = 0; i < cnt_e; ++i) {
		edges[i].clear();
	}
	for (int i = 0; i < cnt_c; ++i) {
		cycles[i].clear();
	}
}

// Convert adjacency matrix to DAG

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat convertToDAG(arma::mat adj, int topE) {
	n = adj.n_rows;
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			adj_matrix[i][j] = adj(i, j);

// initialization
	init();
	
// remove self cycle (and A-B-A cycles)
	preprocess();

// get topE edges in adj matrix
	getTopE(topE);
	
// construct directed graph
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j) {
			if (std::abs(adj_matrix[i][j]) > eps) {
				edges[cnt_e] = Edge(i, j, std::abs(adj_matrix[i][j]));
				edges[cnt_e].next = first[i];
				first[i] = cnt_e;
				cnt_e++;
			}
		}
	
// find strongly connected components (scc) with tarjan, O(n+e)
	for (int i = 0; i < n; ++i)
		if (dfn[i] == -1)
			tarjan(i);
		
// print sccs with >1 nodes
///*
	Rcpp::Rcout << "Found " << cnt_scc << " strongly connected components." << std::endl;
/*
		for (int i = 0; i < cnt_scc; ++i)  {
		Rcpp::Rcout << "scc[" << i << "].size = " << scc[i].size() << std::endl;
		while (!scc[i].empty()) {
			Rcpp::Rcout << scc[i].top() << " ";
			scc[i].pop();
		}
		Rcpp::Rcout << std::endl;

	}
*/
//*/
	
///*
// johnson algorithm to find all elementary circuits of a directed graph, O((n+e)(c+1))
// https://epubs.siam.org/doi/10.1137/0204007
// search for new cycles starting from i (contains i)
	bool find = 0;
	for (int i = 0; i < n; ++i) {
		if (color[i] != -1) {
			johnsonClear();
			find = johnsonFindCycles(i,i);
		}
		visited[i] = 1;
	}
	Rcpp::Rcout << "Found " << cnt_c << " elementary cycles" << std::endl;
//	for (int i = 0; i < cnt_c; ++i)
//		cycles[i].print();
//*/
	removeAllCycles();
	
	arma::mat result(n, n);
	
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
//			Rcpp::Rcout << adj_matrix[i][j] << " ";
			if (removed[i][j]) {
				adj_matrix[i][j] = 0;
			}
			result(i, j) = adj_matrix[i][j];
		}
//		Rcpp::Rcout << std::endl;
	}

// clear all the edges and cycles
	clear();
	
	return result;
}