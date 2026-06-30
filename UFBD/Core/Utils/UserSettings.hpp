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
        unsigned long               getChainLength(void) { return chainLength; }
        int                         getNumCoupledChains(void) { return numCoupledChains; }
        int                         getNumRuns(void) { return numRuns; }
        bool                        getResume(void) { return resume; }
        bool                        getAutoChainLength(void) { return autoChainLength; }
        double                      getBurninFraction(void) { return burninFraction; }
        double                      getRhatThreshold(void) { return rhatThreshold; }
        double                      getEssThreshold(void) { return essThreshold; }
        int                         getCheckEverySamples(void) { return checkEverySamples; }
        unsigned long               getMaxGen(void) { return maxGen; }
        int                         getNumThreads(void) { return numThreads; }
        int                         getNumCores(void) { return numCores; }
        std::string                 getTreeOutput(void) { return treeOut; }
        std::string                 getParamOutput(void) { return parametersOut; }
        std::string                 getTreeFile(void) { return treeFile; }
        std::string                 getCladesFile(void) { return cladesFile; }
        std::string                 getFossilFile(void) { return fossilFile; }
        Conditioning                getConditioning(void) { return conditioning; }
        ConditioningEvent           getConditioningEvent(void) { return conditioningEvent; }
        Probability::PriorFamily         getConditionAgePrior(void) { return conditionAgePrior; }
        double                      getConditionAgePriorP1(void) { return conditionAgePriorP1; }
        double                      getConditionAgePriorP2(void) { return conditionAgePriorP2; }
        bool                        getConditionAgePriorSet(void) { return conditionAgePriorSet; }
        Probability::PriorSpec      getLambdaPrior(void) { return lambdaPrior; }
        Probability::PriorSpec      getMuPrior(void) { return muPrior; }
        Probability::PriorSpec      getPsiPrior(void) { return psiPrior; }
        RateMode                    getLambdaMode(void) { return lambdaMode; }
        RateMode                    getMuMode(void) { return muMode; }
        RateMode                    getPsiMode(void) { return psiMode; }
        double                      getHsmrfShifts(void) { return hsmrfShifts; }
        double                      getHsmrfShiftSize(void) { return hsmrfShiftSize; }
        bool                        getCpuTime(void) { return cpuTime; }
        std::vector<double>         getSkylineTimes(void);
        std::vector<double>         getLambdaSkylineTimes(void) { return lambdaSkylineTimes; }
        std::vector<double>         getMuSkylineTimes(void)     { return muSkylineTimes; }
        std::vector<double>         getPsiSkylineTimes(void)    { return psiSkylineTimes; }
        Model                       getModel(void) { return model; }
        double                      getRho(void) { return rho; }
        unsigned int                getSeed(void) { return seed; }
        bool                        getSeedSet(void) { return seedSet; }
        int                         getPrintFrequency(void) { return printFrequency; }
        int                         getSampleFrequency(void) { return sampleFrequency; }
        std::string                 getHessianFile(void) { return hessianFile; }
        std::string                 getClockModelName(void) { return clockModelName; }
        int                         getSigma2Pncp(void) { return sigma2Pncp; }
        int                         getNStates(void) { return nStates; }
        int                         getModelNStates(void) { return datatypeProvided ? (seqDataType == "aa" ? 20 : 4) : nStates; }
        std::string                 getSequenceFile(void) { return sequenceFile; }
        std::string                 getPartitionFile(void) { return partitionFile; }
        int                         getNumCats(void) { return numCats; }
        std::string                 getSeqDataType(void) { return seqDataType; }
        std::string                 getSubstModel(void) { return substModel; }
        std::string                 getFreqMode(void) { return freqMode; }
        bool                        getUseInvariant(void) { return useInvariant; }
        const double*               getRgeneGamma(void) { return rgeneGamma; }
        const double*               getSigma2Gamma(void) { return sigma2Gamma; }
        void                        print(void);
        void                        printHelp(void);

    private:
                                    UserSettings(void) {}
                                    UserSettings(const UserSettings& u);
        UserSettings&               operator=(const UserSettings&);
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
        std::vector<double>         lambdaSkylineTimes;
        std::vector<double>         muSkylineTimes;
        std::vector<double>         psiSkylineTimes;
        Model                       model;
        double                      rho;
        unsigned int                seed;
        bool                        seedSet;
        unsigned long               chainLength;
        int                         numCoupledChains;
        int                         numRuns;
        bool                        resume;
        bool                        autoChainLength;
        double                      burninFraction;
        double                      rhatThreshold;
        double                      essThreshold;
        int                         checkEverySamples;
        unsigned long               maxGen;
        int                         numThreads;
        int                         numCores;
        int                         printFrequency;
        int                         sampleFrequency;
        std::string                 hessianFile;
        std::string                 clockModelName;
        int                         sigma2Pncp;
        int                         nStates;
        std::string                 sequenceFile;
        std::string                 partitionFile;
        int                         numCats;
        std::string                 seqDataType;
        std::string                 substModel;
        std::string                 freqMode;
        bool                        useInvariant;
        bool                        nstatesProvided;
        bool                        datatypeProvided;
        bool                        coresProvided;
        double                      rgeneGamma[3];
        double                      sigma2Gamma[3];
};

#endif
