#ifndef EmpiricalRateMatrix_hpp
#define EmpiricalRateMatrix_hpp

#include <string>
#include <vector>

class EmpiricalRateMatrix {

    public:
                                    EmpiricalRateMatrix(const std::string& pamlFile, int numStates);
        const std::vector<double>&  getExchangeabilities(void) const { return exch; }
        const std::vector<double>&  getFrequencies(void) const { return freq; }

    private:
        std::vector<double>         exch;
        std::vector<double>         freq;
};

#endif
