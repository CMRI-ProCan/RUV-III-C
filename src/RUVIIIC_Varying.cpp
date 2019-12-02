#include "Rcpp.h"
#include "RcppEigen.h"
#include <limits>
#include <SymEigs.h>

/// @brief Apply RUV-III-C, with a varying set of control variables
///
/// @param input The input data matrix to correct
/// @param k The number of factors of unwanted variation to remove
/// @param M The replicate matrix
/// @param controls The names of the negative control variables
/// @param toCorrect The names of the variables to correct
/// @param withW Should we return the values of the W matrices for every corrected variable? Currently not implemented.
/// 
/// @details Written with reference to the R function RUVIIIC::RUVIII_C_Varying, which is itself based on ruv::RUVIII. See the following papers for further detail on the algorithm and parameters:
/// 	Molania, R., Gagnon-Bartsch, J. A., Dobrovic, A., and Speed, T. P. (2019). A new normalization for Nanostring nCounter gene expression data. Nucleic Acids Research, 47(12), 6073–6083.
///	Gagnon-Bartsch, J. A. and Speed, T. P. (2012). Using control genes to correct for unwanted variation in microarray data. Biostatistics, 13(3), 539–552.
///	Gagnon-Bartsch, J. A., Jacob, L., and Speed, T. P. (2013). Removing unwanted variation from high dimensional data with negative controls.
// [[Rcpp::export]]
Rcpp::RObject RUVIIIC_Varying(Rcpp::NumericMatrix input, int k, Rcpp::NumericMatrix M, Rcpp::CharacterVector controls, Rcpp::CharacterVector toCorrect, bool withW)
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
	int nCorrections = static_cast<int>(toCorrect.size());
	//Number of samples in the data matrix
	int nRows = static_cast<int>(input.nrow());
	int nColumns = static_cast<int>(input.ncol());

	//Vector of W results, if required.
	std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> > WValues(nCorrections);
	std::vector<int> residualDimensions(nCorrections, -1);

	std::vector<std::string> controlNames = Rcpp::as<std::vector<std::string> >(controls);
	//Indices of control variables within the input matrix
	std::vector<int> controlIndices;
	for(std::string control : controlNames)
	{
		controlIndices.push_back(static_cast<int>(std::distance(dataMatrixColumnNames.begin(), std::find(dataMatrixColumnNames.begin(), dataMatrixColumnNames.end(), control))));
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
	inputAsRowMajorImputed = inputAsRowMajorImputed.unaryExpr([](double v){return std::isnan(v) ? 0.0 : v; });

	Eigen::MatrixXd inputSymmetrised = inputAsRowMajorImputed * inputAsRowMajorImputed.transpose();

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
		//Matrix that selects out the relevant rows of the overall data matrix. E.g. submatrixData = selectRowsFromInput * inputAsRowMajorImputed
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> selectRowsFromInput(nRows, nRows);

		//The indices of the negative control variables used, to normalize the current variable
		std::vector<int> controlIndicesThisVariable;
		//The indices of the rows for which the target variable has non-missing values
		std::vector<int> nonMissingIndices;

		//The number of rows copied into the submatrix
		int nSubmatrixRows;
		#pragma omp for schedule(dynamic)
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
			nonMissingIndices.clear();
			for(int row = 0; row < nRows; row++)
			{
				if(!std::isnan(inputAsEigen(row, columnIndexWithinInput)))
				{
					//Copy row of data matrix
					submatrixData.row(nSubmatrixRows) = inputAsRowMajorImputed.row(row);
					//copy row of M
					submatrixM.row(nSubmatrixRows) = MAsEigen.row(row);
					nSubmatrixRows++;
					nonMissingIndices.push_back(row);
				}
			}
			//Create a block view of submatrixData, accounting for th efact that we're not using some of the rows. 
			auto submatrixDataView = submatrixData.block(0, 0, nSubmatrixRows, nColumns);
			//Clear out the list of control variables being used
			controlIndicesThisVariable.clear();
			for(int potentialControlIndex : controlIndices)
			{
				bool useThisControl = true;
				for(int row = 0; row < nSubmatrixRows; row++)
				{
					if(inputAsEigen(nonMissingIndices[row], potentialControlIndex) != inputAsEigen(nonMissingIndices[row], potentialControlIndex))
					{
						useThisControl = false;
						break;
					}
				}
				if(useThisControl) controlIndicesThisVariable.push_back(potentialControlIndex);
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
			int currentResidualDimensions = nSubmatrixRows - columnsMReduced;
			//Check if we even have enough dimensions in the matrix to remove k factors
			if(std::min(currentResidualDimensions, static_cast<int>(controlIndicesThisVariable.size())) < k)
			{
				resultsAsEigen.col(i) = Eigen::VectorXd::Constant(nRows, std::numeric_limits<double>::quiet_NaN());
				continue;
			}
			//In the special case that we're trying to remove the maximum possible number of factors, we can use the GLS formulation. 
			else if (currentResidualDimensions == k)
			{
				//We need a view, to account for the fact that we're only using part of submatrixDataControls
				auto submatrixDataControlsView = submatrixDataControls.block(0, 0, nSubmatrixRows, controlIndicesThisVariable.size());
				//Further subset fullalpha, looking at just the columns corresponding to the control variables
				for(int controlCounter = 0; controlCounter < static_cast<int>(controlIndicesThisVariable.size()); controlCounter++)
				{
					submatrixDataControlsView.col(controlCounter) = submatrixDataView.col(controlIndicesThisVariable[controlCounter]);
				}
				auto productControls = submatrixDataControlsView * submatrixDataControlsView.transpose();
				auto invProductControls = productControls.inverse();
				auto corrected = submatrixMReducedColumnsView * (submatrixMReducedColumnsView.transpose() * invProductControls * submatrixMReducedColumnsView).inverse() * submatrixMReducedColumnsView.transpose() * invProductControls * submatrixDataView.col(columnIndexWithinInput);
				//Now store the results. We're putting the results back into the correct places in the results matrix.
				int correctedCounter = 0;
				for(int row = 0; row < nRows; row++)
				{
					if(std::isnan(input(row, columnIndexWithinInput)))
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
			else
			{
				//Create the view which has the correct number of rows
				auto selectRowsFromInputView = selectRowsFromInput.block(0, 0, nSubmatrixRows, nRows);
				selectRowsFromInputView = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Constant(nSubmatrixRows, nRows, 0);
				for(int i = 0; i < nonMissingIndices.size(); i++)
				{
					selectRowsFromInputView(i, nonMissingIndices[i]) = 1;
				}
				//Take projection of submatrixDataView on the subspace orthogonal to submatrixMReducedColumnsView
				auto orthogonalProjection = Eigen::MatrixXd::Identity(nSubmatrixRows, nSubmatrixRows) - submatrixMReducedColumnsView * (submatrixMReducedColumnsView.transpose() * submatrixMReducedColumnsView).inverse() * submatrixMReducedColumnsView.transpose();
				//Symmetrise. It seems that if you don't explicitly set this type, things go screwy. The automatic type would be some kind of expression template, and I guess that when that goes directly into the Spectra::DenseSymMatProd below, something unexpected happens?
				Eigen::MatrixXd symmetrisedOrthogonal = orthogonalProjection * selectRowsFromInputView * inputSymmetrised * selectRowsFromInputView.transpose() * orthogonalProjection.transpose();
				//Tell Spectra the type used for matrix operations. 
				Spectra::DenseSymMatProd<double> spectraMatrix(symmetrisedOrthogonal);
				//Parameter relating to convergence. Using the default recommended value here.
				int nConvergence = std::min(static_cast<int>(symmetrisedOrthogonal.rows()), std::max(2 * k + 1, 20));
				try
				{
					Spectra::SymEigsSolver<double, Spectra::WHICH_LM, Spectra::DenseSymMatProd<double> > eigs(&spectraMatrix, k, nConvergence);
					eigs.init();
					int nFoundEigenValues = eigs.compute(1000, 1e-10);
					//Check that we found the right number of eigenvalues
					if(nFoundEigenValues < k)
					{
						#pragma omp critical
						{
							std::cout << "Not enough eigenvalues converged for variable " << currentToCorrect << ", this variable was not normalised." << std::endl;
							for(int row = 0; row < nRows; row++)
							{
								resultsAsEigen(row, i) = std::numeric_limits<double>::quiet_NaN();
							}
						}
					}
					else
					{
						//get out the eigenvectors
						Eigen::MatrixXd eigenVectors = eigs.eigenvectors();
						//Get out alpha - effects on all the variables
						Eigen::MatrixXd alpha = eigenVectors.transpose() * submatrixDataView;
						//We need a view, to account for the fact that we're only using part of ac
						auto acView = ac.block(0, 0, alpha.rows(), controlIndicesThisVariable.size());
						//We need a view, to account for the fact that we're only using part of submatrixDataControls
						auto submatrixDataControlsView = submatrixDataControls.block(0, 0, nSubmatrixRows, controlIndicesThisVariable.size());
						//Further subset fullalpha, looking at just the columns corresponding to the control variables
						for(int controlCounter = 0; controlCounter < static_cast<int>(controlIndicesThisVariable.size()); controlCounter++)
						{
							acView.col(controlCounter) = alpha.col(controlIndicesThisVariable[controlCounter]);
							submatrixDataControlsView.col(controlCounter) = submatrixDataView.col(controlIndicesThisVariable[controlCounter]);
						}
						Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> currentPeptideW = submatrixDataControlsView * acView.transpose() * (acView * acView.transpose()).inverse();
						//Now store the results. We're putting the results back into the correct places in the results matrix.
						auto corrected = submatrixDataView.col(columnIndexWithinInput) - currentPeptideW * alpha.col(columnIndexWithinInput);
						//Now store the results. We're putting the results back into the correct places in the results matrix.
						int correctedCounter = 0;
						for(int row = 0; row < nRows; row++)
						{
							if(std::isnan(input(row, columnIndexWithinInput)))
							{
								resultsAsEigen(row, i) = std::numeric_limits<double>::quiet_NaN();
							}
							else 
							{
								resultsAsEigen(row, i) = corrected(correctedCounter);
								correctedCounter++;
							}
						}
						if(withW)
						{
							WValues[i].swap(currentPeptideW);
						}
					}
				}
				catch(...)
				{
					#pragma omp critical
					{
						std::cout << "Error in eigen decomposition for variable " << currentToCorrect << ", this variable was not normalised." << std::endl;
						for(int row = 0; row < nRows; row++)
						{
							resultsAsEigen(row, i) = std::numeric_limits<double>::quiet_NaN();
						}
					}
				}
			}
			if(withW)
			{
				residualDimensions[i] = currentResidualDimensions;
			}
		}
	}
	if(foundError) throw std::runtime_error("Error in RUVIIIC_Varying, check output");
	if(withW)
	{
		Rcpp::List returnValue;
		Rcpp::List WValues_R;
		for(int i = 0; i < nCorrections; i++)
		{
			if(WValues[i].rows() > 0)
			{
				WValues_R(toCorrectNames[i]) = WValues[i];
			}
		}
		returnValue["W"] = WValues_R;
		returnValue["newY"] = results;
		Rcpp::IntegerVector wrappedResidualDimensions = Rcpp::wrap(residualDimensions);
		wrappedResidualDimensions.names() = Rcpp::wrap(toCorrectNames);
		returnValue["residualDimensions"] = wrappedResidualDimensions;
		return returnValue;
	}
	else return results;
}
