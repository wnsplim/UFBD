#ifndef Convergence_hpp
#define Convergence_hpp

#include <vector>

namespace Convergence {

    struct Diagnostic {
        double rhat;
        double bulkEss;
        double tailEss;
    };

    Diagnostic diagnose(const std::vector<std::vector<double>>& chains);

    bool meetsThresholds(const Diagnostic& d, double rhatMax, double essMin);
}

#endif
