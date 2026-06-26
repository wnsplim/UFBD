#include "EmpiricalRateMatrix.hpp"
#include "Msg.hpp"

#include <fstream>
#include <string>
#include <vector>

EmpiricalRateMatrix::EmpiricalRateMatrix(const std::string& pamlFile, int numStates){
    std::ifstream f(pamlFile.c_str());
    if(f.is_open() == false)
        Msg::error("could not open empirical rate matrix file '" + pamlFile + "'");

    int n = numStates;
    int nLower = n * (n - 1) / 2;
    std::vector<double> lower;
    double v;
    while((int)lower.size() < nLower && (f >> v))
        lower.push_back(v);
    while((int)freq.size() < n && (f >> v))
        freq.push_back(v);
    if((int)lower.size() < nLower || (int)freq.size() < n)
        Msg::error("empirical rate matrix file '" + pamlFile + "' must contain " +
                   std::to_string(nLower) + " lower-triangular exchangeabilities followed by " +
                   std::to_string(n) + " equilibrium frequencies");

    exch.assign(nLower, 0.0);
    int k = 0;
    for(int i = 0; i < n; i++)
        for(int j = i + 1; j < n; j++)
            exch[k++] = lower[j * (j - 1) / 2 + i];
}
