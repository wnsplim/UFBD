#ifndef PathModel_hpp
#define PathModel_hpp

#include <map>
#include <string>
#include <utility>
#include <vector>

#include "Eigen/Dense"
#include "FBDInput.hpp"
#include "PhylogeneticModel.hpp"

class FBDTreeModel;
class Node;
class ParameterOUField;
class RelaxedClockTreeModel;
struct SkylineRates;

enum class PseudoPriorCoord : char { LOG, IDENTITY, STANDARDIZED_DEVIATION };

struct PseudoPrior {
    std::vector<std::string>    names;
    Eigen::VectorXd             mean;
    Eigen::LLT<Eigen::MatrixXd> chol;
    std::vector<PseudoPriorCoord> coord;
    std::vector<int>            index;
    std::vector<std::pair<ParameterOUField*,int>> fields;
    Eigen::MatrixXd             K;
    Eigen::VectorXd             wMean;
    Eigen::MatrixXd             cov;
};

class PathModel : public PhylogeneticModel {

    public:
                                    PathModel(void) = delete;
                                    PathModel(PhylogeneticModel* a, PhylogeneticModel* b, std::vector<Fossil>& fossils, unsigned int seed);
        std::vector<std::string>    getParameterNames(void);
        std::vector<double>         getParameterString(void);
        std::vector<bool>           getParameterFixedMask(void);
        std::vector<std::string>    getLatentNames(void);
        std::vector<double>         getLatentString(void);
        std::vector<std::string>    getZoneLegend(void);
        bool                        treeIncludesFossils(void);
        int                         getNumFossils(void);
        double                      lnLikelihood(void);
        double                      lnPriorProbability(void);
        void                        invalidateLikelihoodCache(void);
        void                        invalidatePriorCache(void);
        int                         getLastMoveType(void);
        void                        setHeat(double h);
        void                        print(void);
        double                      update(void);
        void                        updateForAcceptance(void);
        void                        updateForRejection(void);
        void                        writeState(std::ostream& os);
        void                        readState(std::istream& is);
        void                        setStart(int run);
        static std::map<std::string,double> readLastRow(const std::string& path);
        static std::vector<std::vector<double>> readColumns(const std::string& path, const std::vector<std::string>& names, const std::string& key);
        static double               posteriorBurn(const std::string& base);
        static double               lnTypeLabels(const std::vector<std::string>& types, std::vector<Fossil>& fossils);
        static void                 fitRatePseudoPrior(FBDTreeModel* fm, PseudoPrior& d, const std::vector<bool>& spec, const std::string& path, double burn, const std::string& model);
        static double               lnPseudoPrior(FBDTreeModel* f, const PseudoPrior& d);
        static void                 writePseudoPriorState(std::ostream& os, const PseudoPrior& d);
        static int                  readPseudoPriorState(std::istream& is, PseudoPrior& d);

    private:
        double                      dataPart(int k);
        size_t                      rateEnd(void);
        std::vector<bool>           specificMask(FBDTreeModel* f);
        double                      lnRatePrior(FBDTreeModel* f);
        double                      lnPseudoPriorTerm(int k);
        void                        sync(int from);
        void                        setState(int k, const std::map<std::string,double>& mainRow, const std::map<std::string,double>& latRow, const std::string& mainLog, const std::string& latentLog);
        static void                 setRateString(SkylineRates& s, const std::vector<double>& v);
        std::vector<std::string>    clockNames(int k);
        std::vector<double>         clockValues(int k);
        void                        setClockValues(int k, const std::vector<double>& v);
        double                      specificPrior(int k);
        std::vector<Node*>          varyingAgeNodes(int k);
        Eigen::VectorXd             logAges(int k);
        double                      lnSubstPseudoPrior(int k);
        void                        readSubstPseudoPrior(int k, const std::string& path, double burn);
        double                      clockCoords(int k, const std::vector<double>& muV, const std::vector<double>& s2V, const std::vector<std::vector<double>>& r, Eigen::VectorXd& z);
        double                      lnClockPseudoPrior(int k);
        void                        readClockPseudoPrior(int k, const std::string& base, double burn);
        PhylogeneticModel*          m[2];
        FBDTreeModel*               f[2];
        RelaxedClockTreeModel*      rc[2];
        PseudoPrior                 pp[2];
        PseudoPrior                 ppc[2];
        PseudoPrior                 pps[2];
        double                      lab[2];
        bool                        shared[3];
        bool                        sharedClock, sharedSubst;
        double                      pathT;
        int                         mover = 0;
};

#endif
