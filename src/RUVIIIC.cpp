#include "Rcpp.h"
#include "RcppEigen.h"
#include <limits>
#include <SymEigs.h>

// [[Rcpp::export]]
Rcpp::NumericMatrix RUVIIIC(Rcpp::NumericMatrix input, int k, Rcpp::NumericMatrix M, Rcpp::CharacterVector controls, Rcpp::CharacterVector toCorrect, bool withW)
{
	//Result matrix
	Rcpp::NumericMatrix results(input.nrow(), toCorrect.size());
	Rcpp::colnames(results) = toCorrect;
	Rcpp::rownames(results) = Rcpp::rownames(input);
	//Eigen view of the results matrix
	Eigen::Map<Eigen::MatrixXd> resultsAsEigen(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(results));
	//Names of columns in the input
	std::vector<std::string> dataMatrixColumnNames = Rcpp::as<std::vector<std::string> >(colnames(input));
	//Number of peptides to correct
	int nCorrections = (int)toCorrect.size();
	//Number of samples in the data matrix
	int nRows = (int)input.nrow();
	int nColumns = (int)input.ncol();

	std::vector<std::string> controlNames = Rcpp::as<std::vector<std::string> >(controls);
	//Indices of control variables within the input matrix
	std::vector<int> controlIndices;
	for(std::string control : controlNames)
	{
		controlIndices.push_back((int)std::distance(dataMatrixColumnNames.begin(), std::find(dataMatrixColumnNames.begin(), dataMatrixColumnNames.end(), control)));
	}
	std::vector<std::string> toCorrectNames;
	for(int i = 0; i < nCorrections; i++)
	{
		toCorrectNames.push_back(Rcpp::as<std::string>(toCorrect[i]));
	}

	//Create an Eigen view of the matrix
	Eigen::Map<Eigen::MatrixXd> inputAsEigen(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(input));

	//Create an Eigen view of the M matrix
	Eigen::Map<Eigen::MatrixXd> MAsEigen(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(M));

	//Convert to row major, because we're going to be copying rows, so contiguous rows are better.
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> inputAsRowMajorImputed = inputAsEigen;

	//And we're also going to replace NAs with 0s in this row major copy
	inputAsRowMajorImputed = inputAsRowMajorImputed.unaryExpr([](double v){return v == v ? v : 0.0; });

	//Was there an error found during the program? We need a flag, because we can't throw exceptions from inside a parallel loop
	bool foundError = false;
	std::string error;
	#pragma omp parallel
	{
		//We copy the submatrix of the data matrix into here
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> submatrixData(nRows, nColumns);
		//Further reduce the submatrix to just the controls
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> submatrixDataControls(nRows, controlIndices.size());
		//We copy the submatrix of the design matrix into here
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> submatrixM(nRows, M.ncol());
		//And then subsetting down to the columns that contain non-zero values
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> submatrixMReducedColumns(nRows, M.ncol());
		//Matrix used to extract the effects on the control variables
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> ac(nRows, controlIndices.size());

		//The number of rows copied into the submatrix
		int nSubmatrixRows;
		#pragma omp for
		for(int i = 0; i < nCorrections; i++)
		{
			//If there's been an error somewhere, just NOP out the rest of the loop iterations
			if (foundError) continue;
			//The name of the column we are currently trying to correct
			std::string currentToCorrect = toCorrectNames[i];
			//The index of that column within the data matrix
			int columnIndexWithinInput = std::distance(dataMatrixColumnNames.begin(), std::find(dataMatrixColumnNames.begin(), dataMatrixColumnNames.end(), currentToCorrect));
			//This is going to count how many rows of the data matrix we're keeping
			nSubmatrixRows = 0;
			for(int row = 0; row < nRows; row++)
			{
				if(input(row, columnIndexWithinInput) == input(row, columnIndexWithinInput))
				{
					//Copy row of data matrix
					submatrixData.row(nSubmatrixRows) = inputAsRowMajorImputed.row(row);
					//copy row of M
					submatrixM.row(nSubmatrixRows) = MAsEigen.row(row);
					nSubmatrixRows++;
				}
			}
			//Now we need to work out how many columns of submatrixM contain a non-zero value, and keep only those columns
			auto submatrixMView = submatrixM.block(0, 0, nSubmatrixRows, M.ncol());
			Eigen::VectorXd columnSums = Eigen::VectorXd::Ones(nSubmatrixRows).transpose() * submatrixMView;
			//This is the number of columns in the column-reduced version of M (which has already been row-reduced)
			int columnsMReduced = 0;
			for(int columnCounterM = 0; columnCounterM < M.ncol(); columnCounterM++)
			{
				if(columnSums(columnCounterM) > 0)
				{
					submatrixMReducedColumns.col(columnsMReduced) = submatrixMView.col(columnCounterM);
					columnsMReduced++;
				}
			}
			//Create a block view of submatrixMReducedColumns, accounting for the fact that we're not using some of the columns and rows.
			auto submatrixMReducedColumnsView = submatrixMReducedColumns.block(0, 0, nSubmatrixRows, columnsMReduced);
			//Check if we even have enough dimensions in the matrix to remove k factors
			if(std::min(nSubmatrixRows - columnsMReduced, (int)controls.size()) < k)
			{
				resultsAsEigen.col(i) = Eigen::VectorXd::Constant(nRows, std::numeric_limits<double>::quiet_NaN());
				continue;
			}
			//Create a block view of submatrixData, accounting for th efact that we're not using some of the rows. 
			auto submatrixDataView = submatrixData.block(0, 0, nSubmatrixRows, nColumns);
			//Take projection of submatrixDataView on the subspace orthogonal to submatrixMReducedColumnsView
			auto orthogonalProjection = submatrixDataView - submatrixMReducedColumnsView * (submatrixMReducedColumnsView.transpose() * submatrixMReducedColumnsView).inverse() * submatrixMReducedColumnsView.transpose() * submatrixDataView;
			//Symmetrise. It seems that if you don't explicitly set this type, things go screwy. The automatic type would be some kind of expression template, and I guess that when that goes directly into the Spectra::DenseSymMatProd below, something unexpected happens?
			Eigen::MatrixXd symmetrisedOrthogonal = orthogonalProjection * orthogonalProjection.transpose();
			//Tell Spectra the type used for matrix operations. 
			Spectra::DenseSymMatProd<double> spectraMatrix(symmetrisedOrthogonal);
			//Parameter relating to convergence. Using the default recommended value here.
			int nConvergence = std::min(columnsMReduced, std::max(2 * k + 1, 20));
			Spectra::SymEigsSolver<double, Spectra::WHICH_LM, Spectra::DenseSymMatProd<double> > eigs(&spectraMatrix, k, nConvergence);
			eigs.init();
			int nFoundEigenValues = eigs.compute(1000, 1e-10);
			//Check that we found the right number of eigenvalues
			if(nFoundEigenValues < k)
			{
				#pragma omp critical
				{
					if(!foundError)
					{
						foundError = true;
						error = "Not enough eigenvalues converged";
					}
				}
			}
			else
			{
				//get out the eigenvectors
				Eigen::MatrixXd eigenVectors = eigs.eigenvectors();
				//Get out fullalpha - effects on all the variables
				Eigen::MatrixXd alpha = eigenVectors.transpose() * submatrixDataView;
				//Reduce fullalpha to effects of just k factors on all the variables
				//auto alpha = fullalpha.block(0, 0, std::min(k, (int)fullalpha.rows()), fullalpha.cols());
				//We need a view, to account for the fact that we're only using part of ac
				auto acView = ac.block(0, 0, alpha.rows(), controlIndices.size());
				//We need a view, to account for the fact that we're only using part of submatrixDataControls
				auto submatrixDataControlsView = submatrixDataControls.block(0, 0, nSubmatrixRows, controlIndices.size());
				//Further subset fullalpha, looking at just the columns corresponding to the control variables
				for(int controlCounter = 0; controlCounter < (int)controlIndices.size(); controlCounter++)
				{
					acView.col(controlCounter) = alpha.col(controlIndices[controlCounter]);
					submatrixDataControlsView.col(controlCounter) = submatrixDataView.col(controlIndices[controlCounter]);
				}
				auto currentPeptideW = submatrixDataControlsView * acView.transpose() * (acView * acView.transpose()).inverse();
				//Now store the results. We're putting the results back into the correct places in the results matrix.
				auto corrected = submatrixDataView.col(columnIndexWithinInput) - currentPeptideW * alpha.col(columnIndexWithinInput);
				//Now store the results. We're putting the results back into the correct places in the results matrix.
				int correctedCounter = 0;
				for(int row = 0; row < nRows; row++)
				{
					if(input(row, columnIndexWithinInput) != input(row, columnIndexWithinInput))
					{
						resultsAsEigen(row, i) = std::numeric_limits<double>::quiet_NaN();
					}
					else 
					{
						resultsAsEigen(row, i) = corrected(correctedCounter);
						correctedCounter++;
					}
				}
			}
		}
	}
	return results;
}
