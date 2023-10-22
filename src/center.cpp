#include <vector>
#include <unordered_map>
#include <string>
#include <cmath>

void centerConversion(
    std::vector<double>& beta, 
    std::vector<double>& mb_v, 
    const std::unordered_map<std::string, double>& meanValues, 
    const std::string& interaction, 
    const size_t dim,
    const std::vector<int>& betaIntIndex, 
    const std::vector<int>& mbCovIndex) 
{
    // Create ordered mean values
    std::vector<double> orderedMeanValues;
    for (const char& c : interaction) {
        if (meanValues.find(std::string(1, c)) != meanValues.end()) {
            orderedMeanValues.push_back(meanValues.at(std::string(1, c)));
        } else {
            orderedMeanValues.push_back(0.0);  // Default to 0 if not provided
        }
    }

    // Create matrix C
    std::vector<std::vector<double>> C(dim, std::vector<double>(dim, 0.0));
    C[0][0] = 1.0;
    for (size_t i = 1; i < dim; i++) {
        C[i][i] = 1.0;
        C[0][i] = orderedMeanValues[i - 1];
    }

    // Create matrix M
    size_t mDim = dim*dim;
    std::vector<std::vector<double>> M(mDim, std::vector<double>(mDim, 0.0));

    for (size_t i = 0; i < mDim; i++) {
        M[i][i] = 1.0;
    }

    for (size_t i = 1; i < dim; i++) {
        M[0][i] = 2 * orderedMeanValues[i - 1];
        M[0][(i * dim) + i] = std::pow(orderedMeanValues[i - 1], 2);
        M[i][(i * dim) + i] = orderedMeanValues[i - 1];
        for (size_t j = i + 1; j < dim; j++) {
            M[0][(i * dim) + j] = 2 * orderedMeanValues[i - 1] * orderedMeanValues[j - 1];
            M[i][(i * dim) + j] = orderedMeanValues[j - 1];
        }
    }

    // Multiply beta with C
    std::vector<double> newBeta(dim, 0.0);
    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < dim; j++) {
            newBeta[i] += C[i][j] * beta[j];
        }
    }
    beta = newBeta;

    // Multiply mb_v with M
    std::vector<double> newMb_v(mDim, 0.0);
    for (size_t i = 0; i < mDim; i++) {
        for (size_t j = 0; j < mDim; j++) {
            newMb_v[i] += M[i][j] * mb_v[j];
        }
    }
    mb_v = newMb_v;
}



