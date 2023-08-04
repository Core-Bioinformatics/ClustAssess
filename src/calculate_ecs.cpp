#include <Rcpp.h>

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
	double c1, c2, pointECS, intersectSize;
	
	IntegerMatrix contTable = myContTable(mb1, mb2, minim_mb1, minim_mb2);
    int n = mb1.size(), nClusters1 = contTable.nrow(), nClusters2 = contTable.ncol();

	NumericMatrix uniqueECSvals(nClusters1, nClusters2);
	NumericVector ecsScore(n);
	IntegerVector clustersSizes1 = rowSums(contTable), clusterSizes2 = colSums(contTable);
	uniqueECSvals.fill(-1);

	for (int i = 0; i < n; i++) {
		cluster1 = mb1(i) - minim_mb1;
		cluster2 = mb2(i) - minim_mb2;

		if (uniqueECSvals(cluster1, cluster2) != -1) {
			ecsScore(i) = uniqueECSvals(cluster1, cluster2);
			continue;
		}

		c1 = clustersSizes1(cluster1);
		c2 = clusterSizes2(cluster2);
		intersectSize = contTable(cluster1, cluster2);
		pointECS = (1 / c1 - 1 / c2) * intersectSize;
		if (pointECS < 0) {
			pointECS *= -1;
		}

		ecsScore(i) = 1 - 0.5 * (pointECS + (c1 - intersectSize) / c1 + (c2 - intersectSize) / c2);
		uniqueECSvals(cluster1, cluster2) = ecsScore(i);
	}

	return(ecsScore);
}
