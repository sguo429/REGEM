
#include "declars.h"


void printOutputHeader(int main_headers, int printStart, int printEnd, vector<string> covNames, string output, string outStyle);

int main(int argc, char* argv[]) {

	CommandLine cmd;
	cmd.processCommandLine(argc, argv);
	auto start_time = std::chrono::high_resolution_clock::now();
	string inFile  = cmd.inFile;
	string outFile = cmd.outFile;
	vector<string> exp  = cmd.exp;
	vector<string> icov = cmd.icov;
	int expSq  = cmd.numExpSelCol;
	int expSq1 = expSq+1;
	int intSq  = cmd.numIntSelCol;
	int intSq1 = intSq+1;
	int Sq1    = intSq1 + expSq;

	if (intSq > 0 ) {
		for (int i = 0; i < intSq; i++) {
			exp.push_back(icov[i]);
		}
	}
	for (size_t i = 0; i < exp.size(); i++) {
		exp[i] = "G-" + exp[i];
	}
	exp.insert(exp.begin(), "G");

    int printStart = 1; 
	int printEnd   = expSq1;
    bool printFull = false;
    if (cmd.outStyle.compare("meta") == 0) {
        printStart = 0; printEnd = Sq1;
    } else if (cmd.outStyle.compare("full") == 0) {
        printStart = 0; printEnd = Sq1; printFull = true;
    }

	std::ifstream res;
    res.open(inFile);
    if (!res.is_open()) {
        cerr << "\nERROR: Cannot open results file. \n\n" << endl;
        exit(1);
    }

	string line;
    getline(res, line);
    std::istringstream issHead(line);
    string header;
    int header_i = 0;
	std::unordered_map<std::string, int> colNames;
    while (getline(issHead, header, '\t')) {
        header.erase(std::remove(header.begin(), header.end(), '\r'), header.end());
		if (colNames.find(header) != colNames.end()) {
            cerr << "\nERROR: There are duplicate header names (" << header << ") in the results file.\n\n";
            exit(1);
        }
		colNames[header] = header_i;

        if (header.rfind("Beta_G-", 0) == 0) {
            string tmp_name = header;
            tmp_name.erase(0, 5);
            if (std::find(exp.begin(), exp.end(), tmp_name) == exp.end()) {
                exp.push_back(tmp_name);
            }
        }
        header_i++;
    }


	// Get each column header index
	int main_headers = 9;
	if (colNames.find("RSID") != colNames.end()) {
		main_headers = 10;
	}
	int pMarginalIndex = colNames["P_Value_Marginal"];

	int dim = exp.size();
	vector<int> Vindex(dim*dim);
	vector<int> Betaindex(dim);
	vector<int> Varindex(dim*dim);
	for (int i = 0; i < dim; i++) {
		Betaindex[i] = colNames["Beta_" + exp[i]];
		if (Betaindex[i] == 0) {
			cerr << "\nERROR: Column header Beta_" + exp[i] + " does not exist in input file.\n\n";
			exit(1);
		}
		for (int j = 0; j < dim; j++) {
			if (i == j) {
				Varindex[(i * dim) + j] = colNames["Var_Beta_" + exp[i]];
				if (Varindex[(i * dim) + j] == 0) {
					cerr << "\nERROR: Column header Var_Beta_" + exp[i] + " does not exist in the input file.\n\n";
					exit(1);
				}
			} else {
				Varindex[(i * dim) + j] = colNames["Cov_Beta_" + exp[i] + "_" + exp[j]];
				if (Varindex[(i * dim) + j] == 0) {
					cerr << "\nERROR: Column header Cov_Beta_" + exp[i] + "_" + exp[j] + " does not exist in the input file.\n\n";
					exit(1);
				}
			}
			Vindex[(i * dim) + j] = colNames["V_" + exp[i] + "_" + exp[j]];
			if (Vindex[(i * dim) + j] == 0) {
				cerr << "\nERROR: Column header V_" + exp[i] + "_" + exp[j] + " does not exist in the input file.\n\n";
				exit(1);
			}
		}	
	}


    int nvars = 0;
    while (getline(res, line)) nvars++;
    res.clear(); 
	res.seekg(0, res.beg);
    getline(res, line);	


	vector<double> Vvec(dim*dim);
	vector<double> Betavec(dim);
	vector<double> Varvec(dim*dim);
	double* Beta = &Betavec[0];
	double* Var  = &Varvec[0];
	double* V    = &Vvec[0];

	boost::math::chi_squared chisq_dist_M(1);
	boost::math::chi_squared chisq_dist_Int(expSq);
	boost::math::chi_squared chisq_dist_Joint(expSq1);

   
    printOutputHeader(main_headers, printStart, printEnd, exp, outFile, cmd.outStyle);
	std::ofstream results(outFile, std::ios_base::app);
	std::ostringstream oss;

	for (int p = 0; p < nvars; p++) {
		getline(res, line);
		std::istringstream iss(line);
		string value;
		vector <string> values;
		while (getline(iss, value, '\t')) values.push_back(value);
		for (int i = 0; i < dim; i++) {
			sscanf(values[Betaindex[i]].c_str(), "%lf", &Betavec[i]);
			for (int j = 0; j < dim; j++) {
				sscanf(values[Varindex[(i * dim) + j]].c_str(), "%lf", &Varvec[(i * dim) + j]);
				sscanf(values[Vindex[(i * dim) + j]].c_str(), "%lf", &Vvec[(i * dim) + j]);
			}
		}

		double* U = new double[dim];
		matvecprod(V, Beta, U, dim, dim);
		double* S = new double[dim*dim];
		matNmatNprod(V, Var, S, dim, dim, dim);
		double* O = new double[dim*dim];
		matNmatNprod(S, V, O, dim, dim, dim);
		double* O_dot = new double[Sq1*Sq1];
		subMatrix(O, O_dot, Sq1, Sq1, dim, Sq1, 0);


		double* V_i = new double[Sq1*Sq1];
		subMatrix(V, V_i, Sq1, Sq1, dim, Sq1, 0);
		matInv(V_i, Sq1);

		double* betaAll = new double[Sq1];
		matvecprod(V_i, U, betaAll, Sq1, Sq1);

		double* OV_i = new double[Sq1 * Sq1];
        matNmatNprod(O_dot, V_i, OV_i, Sq1, Sq1, Sq1);
		double* VarBetaAll = new double[Sq1*Sq1];
		matNmatNprod(V_i, OV_i, VarBetaAll, Sq1, Sq1, Sq1);

		double* invVarBetaE = new double[expSq * expSq];
		subMatrix(VarBetaAll, invVarBetaE, expSq, expSq, Sq1, expSq, Sq1 + 1);
		matInv(invVarBetaE, expSq);

		double* Stemp3 = new double[expSq];
		matvecSprod(invVarBetaE, betaAll, Stemp3, expSq, expSq, 1);
		double statInt = 0.0;
		for (int j = 1; j < expSq1; j++) {
			statInt += betaAll[j] * Stemp3[j-1];
		}
		double PvalInt = (isnan(statInt) || statInt <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_Int, statInt));

		double* invA = new double[expSq1 * expSq1];
		subMatrix(VarBetaAll, invA, expSq1, expSq1, Sq1, expSq1, 0);
		matInv(invA, expSq1);

		double* Stemp4 = new double[expSq1];
		matvecprod(invA, betaAll, Stemp4, expSq1, expSq1);
		double statJoint = 0.0;
		for (int k = 0; k < expSq1; k++) {
			statJoint += betaAll[k] * Stemp4[k];
		}
			
		double PvalJoint = (isnan(statJoint) || statJoint <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_Joint, statJoint));

		for (int i = 0; i < main_headers; i++) {
			oss << values[i] << "\t";
		}
		for (int ii = printStart; ii < printEnd; ii++) {
            oss << betaAll[ii] << "\t";
        }
		for (int ii = printStart; ii < printEnd; ii++) {
            oss << VarBetaAll[ii * Sq1 + ii] << "\t";
        }
		for (int ii = printStart; ii < printEnd; ii++) {
            for (int jj = printStart; jj < printEnd; jj++) {
                 if (ii != jj) {
                    oss << VarBetaAll[ii * Sq1 + jj] << "\t";
                }
            }
        }
		if (printFull) {
            for (int ii = printStart; ii < printEnd; ii++) {
                for (int jj = printStart; jj < printEnd; jj++) {
                    oss << V[ii*Sq1 + jj] << "\t";
                }
            }
        }
		oss << values[pMarginalIndex] << "\t" << PvalInt << "\t" << PvalJoint << "\n";

		delete[] U;
		delete[] S;
		delete[] O;
		delete[] O_dot;
		delete[] Stemp3;
		delete[] Stemp4;
		delete[] V_i;
		delete[] invA;
		
	}

	results << oss.str();
	oss.str(std::string());
	oss.clear();
	results.close();
	
	auto end_time = std::chrono::high_resolution_clock::now();
	cout << "Finished\n";
	return 0;
}




void printOutputHeader(int main_headers, int printStart, int printEnd, vector<string> covNames, string output, string outStyle) {

    std::ofstream results(output, std::ofstream::binary);
    if (main_headers == 10) {
        results << "SNPID" << "\t" << "RSID" << "\t" << "CHR" << "\t" << "POS" << "\t" << "Non_Effect_Allele" << "\t" << "Effect_Allele" << "\t" << "N_Samples" << "\t" << "AF" << "\t" << "Beta_Marginal" << "\t" << "Var_Beta_Marginal" << "\t";
    }
    else {
        results << "SNPID" << "\t" << "CHR" << "\t" << "POS" << "\t" << "Non_Effect_Allele" << "\t" << "Effect_Allele" << "\t" << "N_Samples" << "\t" << "AF" << "\t" << "Beta_Marginal" << "\t" << "Var_Beta_Marginal" << "\t";
    }

   

    for (int i = printStart; i < printEnd; i++) {
        results << "Beta_" << covNames[i] << "\t";
    }
    for (int i = printStart; i < printEnd; i++) {
        results << "Var_Beta_" << covNames[i] << "\t";  
    }
    for (int i = printStart; i < printEnd; i++) {
        for (int j = printStart; j < printEnd; j++) {
            if (i != j) {
                results << "Cov_Beta_" << covNames[i] << "_" << covNames[j] << "\t";  
            } 
        }
    }
    if (outStyle.compare("full") == 0) {
        for (int i = printStart; i < printEnd; i++) {
            for (int j = printStart; j < printEnd; j++) {
                results << "V_" << covNames[i] << "_" << covNames[j] << "\t";
            }
        }
    }
    results << "P_Value_Marginal" << "\t" << "P_Value_Interaction" << "\t" << "P_Value_Joint\n";


    results.close();
}