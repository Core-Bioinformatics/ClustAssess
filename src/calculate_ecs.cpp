#include <Rcpp.h>
//[[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix myContTable(IntegerVector &a, IntegerVector &b, int minim_mb1, int minim_mb2) {
	int maxim_mb1 = max(a), maxim_mb2 = max(b);
    int nClustersA = maxim_mb1 - minim_mb1 + 1;
	int nClustersB = maxim_mb2 - minim_mb2 + 1, n = a.size();

    IntegerMatrix result(nClustersA, nClustersB);

    for(int i = 0; i < n; i++) {
        result(a(i) - minim_mb1, b(i) - minim_mb2)++;
    }

	return(result);
}

// [[Rcpp::export]]
NumericVector disjointECS(IntegerVector mb1, IntegerVector mb2) {
	int minim_mb1 = min(mb1), minim_mb2 = min(mb2);
	int cluster1, cluster2;
	
	IntegerMatrix contTable = myContTable(mb1, mb2, minim_mb1, minim_mb2);
    int n = mb1.size(), nClusters1 = contTable.nrow(), nClusters2 = contTable.ncol();

	NumericMatrix uniqueECSvals(nClusters1, nClusters2);
	NumericVector ecsScore(n);
	IntegerVector clustersSizes1 = rowSums(contTable), clustersSizes2 = colSums(contTable);

	for (int i = 0; i < nClusters1; i++) {
		for (int j = nClusters2 - 1; j >= 0; j--) {
			uniqueECSvals(i, j) = (double) contTable(i, j) / std::max(clustersSizes1(i), clustersSizes2(j));
		}
	}

	for (int i = 0; i < n; i++) {
		cluster1 = mb1(i) - minim_mb1;
		cluster2 = mb2(i) - minim_mb2;
		ecsScore(i) = uniqueECSvals(cluster1, cluster2);
	}

	return(ecsScore);
}

// [[Rcpp::export]]
double disjointECSaverage(IntegerVector mb1, IntegerVector mb2) {
	// NOTE Eigen doesn't seem to improve performance
	int minim_mb1 = min(mb1), minim_mb2 = min(mb2);
	int cluster1, cluster2;
	
	IntegerMatrix contTable = myContTable(mb1, mb2, minim_mb1, minim_mb2);
    int nClusters1 = contTable.nrow(), nClusters2 = contTable.ncol();
	int n = sum(contTable);

	double avgECS = 0, intersectSize;

	IntegerVector clustersSizes1 = rowSums(contTable), clustersSizes2 = colSums(contTable);

	for (int i = 0; i < nClusters1; i++) {
		for (int j = nClusters2 - 1; j >= 0; j--) {
			intersectSize = contTable(i, j);
			avgECS += intersectSize * intersectSize / (n * std::max(clustersSizes1(i), clustersSizes2(j)));
		}
	}

	return(avgECS);
}
