#ifndef UserSettings_hpp
#define UserSettings_hpp

#include <string>

#include "Probability.hpp"

enum class Conditioning { CROWN, ORIGIN };
enum class ConditioningEvent { SURVIVAL, ANYSAMPLE, EXTINCT };
enum class Model { FBD, HEA14, UFBD };
enum class RateMode { IID, SMOOTH };

class UserSettings {

    public:
        static UserSettings&        userSettings(void)
                                        {
                                        static UserSettings us;
                                        return us;
                                        }
        void                        initializeSettings(int argc, const char* argv[], bool sbcMode = false);
        unsigned long               getChainLength(void) { checkSettings(); return chainLength; }
        int                         getNumChains(void) { return numChains; }
        int                         getNumThreads(void) { return numThreads; }
        std::string                 getTreeOutput(void) { checkSettings(); return treeOut; }
        std::string                 getParamOutput(void) { checkSettings(); return parametersOut; }
        std::string                 getTreeFile(void) { checkSettings(); return treeFile; }
        std::string                 getCladesFile(void) { checkSettings(); return cladesFile; }
        std::string                 getFossilFile(void) { checkSettings(); return fossilFile; }
        Conditioning                getConditioning(void) { checkSettings(); return conditioning; }
        ConditioningEvent           getConditioningEvent(void) { checkSettings(); return conditioningEvent; }
        Probability::PriorFamily         getConditionAgePrior(void) { checkSettings(); return conditionAgePrior; }
        double                      getConditionAgePriorP1(void) { checkSettings(); return conditionAgePriorP1; }
        double                      getConditionAgePriorP2(void) { checkSettings(); return conditionAgePriorP2; }
        bool                        getConditionAgePriorSet(void) { checkSettings(); return conditionAgePriorSet; }
        Probability::PriorSpec      getLambdaPrior(void) { checkSettings(); return lambdaPrior; }
        Probability::PriorSpec      getMuPrior(void) { checkSettings(); return muPrior; }
        Probability::PriorSpec      getPsiPrior(void) { checkSettings(); return psiPrior; }
        RateMode                    getLambdaMode(void) { checkSettings(); return lambdaMode; }
        RateMode                    getMuMode(void) { checkSettings(); return muMode; }
        RateMode                    getPsiMode(void) { checkSettings(); return psiMode; }
        double                      getHsmrfShifts(void) { checkSettings(); return hsmrfShifts; }
        double                      getHsmrfShiftSize(void) { checkSettings(); return hsmrfShiftSize; }
        bool                        getCpuTime(void) { checkSettings(); return cpuTime; }
        std::vector<double>         getSkylineTimes(void) { checkSettings(); return skylineTimes; }
        Model                       getModel(void) { checkSettings(); return model; }
        double                      getRho(void) { checkSettings(); return rho; }
        unsigned int                getSeed(void) { checkSettings(); return seed; }
        bool                        getSeedSet(void) { checkSettings(); return seedSet; }
        int                         getPrintFrequency(void) { checkSettings(); return printFrequency; }
        int                         getSampleFrequency(void) { checkSettings(); return sampleFrequency; }
        std::string                 getHessianFile(void) { checkSettings(); return hessianFile; }
        std::string                 getClockModelName(void) { checkSettings(); return clockModelName; }
        int                         getNStates(void) { checkSettings(); return nStates; }
        std::string                 getSequenceFile(void) { checkSettings(); return sequenceFile; }
        std::string                 getPartitionFile(void) { checkSettings(); return partitionFile; }
        int                         getNumCats(void) { checkSettings(); return numCats; }
        const double*               getRgeneGamma(void) { checkSettings(); return rgeneGamma; }
        const double*               getSigma2Gamma(void) { checkSettings(); return sigma2Gamma; }
        void                        print(void);
        void                        printHelp(void);

    private:
                                    UserSettings(void) { settingsInitialized = false; }
                                    UserSettings(const UserSettings& u);
        UserSettings&               operator=(const UserSettings&);
        void                        checkSettings(void);
        void                        parsePriorInto(const std::string& spec, Probability::PriorFamily& family, double& p1, double& p2);
        std::string                 executablePath;
        std::string                 treeOut;
        std::string                 parametersOut;
        std::string                 treeFile;
        std::string                 cladesFile;
        std::string                 fossilFile;
        Conditioning                conditioning;
        ConditioningEvent           conditioningEvent;
        bool                        conditioningSet;
        Probability::PriorFamily         conditionAgePrior;
        bool                        conditionAgePriorSet;
        double                      conditionAgePriorP1;
        double                      conditionAgePriorP2;
        Probability::PriorSpec      lambdaPrior;
        Probability::PriorSpec      muPrior;
        Probability::PriorSpec      psiPrior;
        RateMode                    lambdaMode;
        RateMode                    muMode;
        RateMode                    psiMode;
        double                      hsmrfShifts;
        double                      hsmrfShiftSize;
        bool                        cpuTime;
        std::vector<double>         skylineTimes;
        Model                       model;
        double                      rho;
        unsigned int                seed;
        bool                        seedSet;
        unsigned long               chainLength;
        int                         numChains;
        int                         numThreads;
        int                         printFrequency;
        int                         sampleFrequency;
        std::string                 hessianFile;
        std::string                 clockModelName;
        int                         nStates;
        std::string                 sequenceFile;
        std::string                 partitionFile;
        int                         numCats;
        double                      rgeneGamma[3];
        double                      sigma2Gamma[3];
        bool                        settingsInitialized;
};

#endif
