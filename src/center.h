#include <vector>
#include <unordered_map>
#include <string>
#include <cmath>

void centerConversion(
    std::vector<double>& beta, 
    std::vector<double>& mb_v, 
    std::unordered_map<std::string, double> meanValues, 
    std::vector<std::string> interaction_names, 
    size_t dim,
    std::vector<int> betaIntIndex, 
    std::vector<int> mbCovIndex,
    int nExp,
    int center_in,
    int center_out);
