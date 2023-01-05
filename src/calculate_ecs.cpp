#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix myContTable(IntegerVector &a, IntegerVector &b) {
    int nClustersA = max(a) + 1, nClustersB = max(b) + 1, n = a.size();
    IntegerMatrix result(nClustersA, nClustersB);
    
    for(int i = 0; i < n; i++) {
        result(a(i), b(i))++;
    }

	return(result);
}


// [[Rcpp::export]]
NumericVector disjointECS(IntegerVector mb1, IntegerVector mb2) {
	int minim_mb1 = min(mb1), minim_mb2 = min(mb2);
	IntegerVector copy_mb1 = mb1;
	if (minim_mb1) {
		copy_mb1 = copy_mb1 - minim_mb1;
	}
	IntegerVector copy_mb2 = mb2;
	if (minim_mb2) {
		copy_mb2 = copy_mb2 - minim_mb2;
	}
	IntegerMatrix contTable = myContTable(copy_mb1, copy_mb2);
    int n = copy_mb1.size(), nClusters1 = contTable.nrow(), nClusters2 = contTable.ncol();

	NumericMatrix uniqueECSvals(nClusters1, nClusters2);
	NumericVector ecsScore(n);
	IntegerVector clustersSizes1 = rowSums(contTable), clusterSizes2 = colSums(contTable);
	uniqueECSvals.fill(-1);

	for (int i = 0; i < n; i++) {
		int cluster1 = copy_mb1(i), cluster2 = copy_mb2(i);

		if (uniqueECSvals(cluster1, cluster2) != -1) {
			ecsScore(i) = uniqueECSvals(cluster1, cluster2);
			continue;
		}

		double c1 = clustersSizes1(cluster1), c2 = clusterSizes2(cluster2), intersectSize = contTable(cluster1, cluster2);
		double pointECS = (1 / c1 - 1 / c2) * intersectSize;
		if (pointECS < 0) {
			pointECS *= -1;
		}

		ecsScore(i) = 1 - 0.5 * (pointECS + (c1 - intersectSize) / c1 + (c2 - intersectSize) / c2);
		uniqueECSvals(cluster1, cluster2) = ecsScore(i);
	}

	return(ecsScore);
}
