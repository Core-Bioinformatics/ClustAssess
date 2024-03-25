#include <Rcpp.h>
#include <RcppEigen.h>
//[[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

// NOTE this function is not used, that has potential to further optimize the SNN calculation for varying k values
// [[Rcpp::export(rng = false)]]
List filterNNmatrix(Eigen::SparseMatrix<double> oldNN, Eigen::MatrixXd nnRanked, int oldK, int newK, double prune = 0) {
	int nRows = nnRanked.rows();

	for (int j = newK; j < oldK; j++) {
		for (int i = 0; i < nRows; i++) {
			oldNN.coeffRef(i, nnRanked(i, j) - 1) = 0;
			// oldNN.insert(i, nnRanked(i, j) - 1) = 1;
		}
	}

	oldNN.prune(0.0);

	Eigen::SparseMatrix<double> SNN = oldNN * (oldNN.transpose());
	for (int i = 0; i < SNN.outerSize(); i++) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(SNN, i); it; ++it) {
			it.valueRef() = it.value()/(newK + (newK - it.value()));
			if (it.value() < prune) {
				it.valueRef() = 0;
			}
		}
	}

	if (prune > 0) {
		SNN.prune(0.0); // actually remove pruned values
	}

	return List::create(Named("nn") = oldNN, Named("snn") = SNN);
}


// [[Rcpp::export(rng = false)]]
Eigen::SparseMatrix<double> pruneSNN(Eigen::SparseMatrix<double> snnMatrix, double prune = 0) {
	if (prune <= 0) {
		return snnMatrix;
	}

	for (int i = 0; i < snnMatrix.outerSize(); i++) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(snnMatrix, i); it; ++it) {
			if (it.value() <= prune) {
				it.valueRef() = 0;
			}
		}
	}

	snnMatrix.prune(0.0); // actually remove pruned values

	return snnMatrix;
}

// [[Rcpp::export(rng = false)]]
Eigen::SparseMatrix<double> computeSNN(Eigen::SparseMatrix<double> &nnMatrix, int k, double prune = 0) {
	Eigen::SparseMatrix<double> SNN = nnMatrix * (nnMatrix.transpose());

	for (int i = 0; i < SNN.outerSize(); i++) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(SNN, i); it; ++it) {
			it.valueRef() = it.value()/(k + (k - it.value()));
			if (prune <= 0) {
				continue;
			}

			if (it.value() <= prune) {
				it.valueRef() = 0;
			}
		}
	}

	if (prune > 0) {
 	 	SNN.prune(0.0); // actually remove pruned values
	}

	return(SNN);
}


//' Computes the NN adjacency matrix given the neighbours
//'
//' @param nnRanked A matrix with the lists of the nearest neighbours for each point
//' @param k The number of neighbours to consider. Defaults to `-1`, which means all neighbours.
//' @param start The index of the first neighbour to consider. Defaults to `0`.
//' @param prune The threshold to prune the SNN matrix. Defaults to `0`.
//'
//' @export
//'
//' @return A list with the NN and SNN adjacency matrices.
// [[Rcpp::export(rng = false)]]
List getNNmatrix(Eigen::MatrixXd nnRanked, int k = -1, int start = 0, double prune = 0) {
	int nRows = nnRanked.rows(), nCols = nnRanked.cols();

	if (k == -1 || k > nCols) {
		k = nCols;
	}

	if (start + k > nCols) {
		k = nCols - start;
	}
  
	// std::cout << nRows << ' ' << nCols << ' ' << k << '\n';
  	std::vector<Eigen::Triplet<double>> tripletList;
  	tripletList.reserve(nRows * k);

  	for (int j = start; j < start + k; j++) {
		for (int i = 0; i < nRows; i++) {
			tripletList.push_back(Eigen::Triplet<double>(i, nnRanked(i, j) - 1, 1));
    	}
  	}

  	Eigen::SparseMatrix<double> NN(nRows, nRows);
	NN.setFromTriplets(tripletList.begin(), tripletList.end());

	if (prune < 0) {
		return List::create(Named("nn") = NN);
	}
	

	return List::create(Named("nn") = NN, Named("snn") = computeSNN(NN, k, prune));
}
