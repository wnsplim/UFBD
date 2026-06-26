#ifndef ParameterDouble_hpp
#define ParameterDouble_hpp

#include "Parameter.hpp"
#include <deque>
#include <iosfwd>
#include <string>
#include <vector>

namespace Probability { enum class PriorFamily; }

class ParameterDouble : public Parameter {

    public:
                                    ParameterDouble(void) = delete;
                                    ParameterDouble(double prob, PhylogeneticModel* m, std::string n, double lb, double up); //for "unbounded RVs" pass in std::numeric_limits<double>::max() / 2 and std::numeric_limits<double>::min() / 2
        double                      getAcceptanceRatio(void) { return ((double)numAcceptances)/((double)(numAcceptances+numRejections));}
        bool                        getAdaptiveProposalActive(void) { return adaptiveProposalActive; }
        double                      getValue(void) { return value[0]; } // 0 is the one we update, 1 is the one we don't (last currently accepted value)
        void                        setValue(double v) { value[0] = v; value[1] = v; }
        void                        scaleProposed(double c) { value[0] = value[1] * c; }
        void                        commitProposed(void) { value[1] = value[0]; }
        void                        restoreProposed(void) { value[0] = value[1]; }
        void                        setPrior(Probability::PriorFamily f, double p1, double p2);
        double                      lnProbability(void);
        void                        print(void);
        double                      update(void);
        void                        updateForAcceptance(void);
        void                        updateForRejection(void);
        void                        writeState(std::ostream& os);
        void                        readState(std::istream& is);
    private:
        double                      updateBactrianScale(void);
        std::vector<double>         value;
        std::deque<bool>            recentAcceptRej;
        double                      lowerBound;
        double                      upperBound;
        Probability::PriorFamily    priorFamily;
        double                      priorP1;
        double                      priorP2;
        double                      targetAr;
        double                      windowSize;
        int                         numAcceptances;
        int                         numAdaptive;
        int                         numRejections;
};


#endif /* ParameterDouble_hpp */
