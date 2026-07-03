#ifndef Parameter_hpp
#define Parameter_hpp

#include <iosfwd>
#include <string>
class PhylogeneticModel;


class Parameter {

    public:
                            Parameter(void) = delete;
                            Parameter(double prob, PhylogeneticModel* m, std::string n);
        virtual double      getAcceptanceRatio(void) = 0;
        std::string         getName(void) { return parmName; }
        double              getProposalProbability(void) { return proposalProbability; }
        virtual double      lnProbability(void) = 0;
        virtual void        print(void) = 0;
        void                setProposalProbability(double x) { proposalProbability = x; }
        virtual double      update(void) = 0;
        virtual void        updateForAcceptance(void) = 0;
        virtual void        updateForRejection(void) = 0;
        virtual void        writeState(std::ostream& os) {}
        virtual void        readState(std::istream& is) {}

    protected:
        std::string         parmName;
        PhylogeneticModel*  model;
        double              proposalProbability;
};

#endif
