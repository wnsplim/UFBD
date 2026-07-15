#ifndef UserSettings_hpp
#define UserSettings_hpp

#include <map>
#include <string>
#include <vector>

#include "Probability.hpp"

enum class Conditioning { CROWN, ORIGIN };
enum class ConditioningEvent { SURVIVAL, ANYSAMPLE, EXTINCT };
enum class Model { RFBD, HEA14, UFBD };
enum class RateMode { IID, SMOOTH };
enum class Sigma2Param { PNCP, C, NC };

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
        bool                        getArLog(void) { return arLog; }
        bool                        getAutoChainLength(void) { return autoChainLength; }
        double                      getBurninFraction(void) { return burninFraction; }
        double                      getRhatThreshold(void) { return rhatThreshold; }
        double                      getEssThreshold(void) { return essThreshold; }
        unsigned long               getMaxGen(void) { return maxGen; }
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
        RateMode                    getLambdaMode(void) { return lambdaMode; }
        RateMode                    getMuMode(void) { return muMode; }
        double                      getHsmrfShifts(void) { return hsmrfShifts; }
        double                      getHsmrfShiftSize(void) { return hsmrfShiftSize; }
        bool                        getCpuTime(void) { return cpuTime; }
        std::vector<double>         getSkylineTimes(void);
        std::vector<double>         getLambdaSkylineTimes(void) { return lambdaSkylineTimes; }
        std::vector<double>         getMuSkylineTimes(void)     { return muSkylineTimes; }
        int                         getNumPsiTypes(void) { return psiTypeNames.empty() ? 1 : (int)psiTypeNames.size(); }
        const std::vector<std::string>& getPsiTypeNames(void) { return psiTypeNames; }
        Probability::PriorSpec      getPsiPrior(int t);
        std::vector<double>         getPsiSkylineTimes(int t);
        RateMode                    getPsiMode(int t);
        std::vector<int>            getLambdaGroups(void) { return lambdaGroups; }
        std::vector<int>            getMuGroups(void) { return muGroups; }
        std::vector<int>            getPsiGroups(int t);
        std::map<int,Probability::PriorSpec> getLambdaGroupPrior(void) { return lambdaGroupPrior; }
        std::map<int,Probability::PriorSpec> getMuGroupPrior(void) { return muGroupPrior; }
        std::map<int,Probability::PriorSpec> getPsiGroupPrior(int t);
        Model                       getModel(void) { return model; }
        double                      getRho(void) { return rho; }
        unsigned int                getSeed(void) { return seed; }
        int                         getThinning(void) { return thinning; }
        std::string                 getHessianFile(void) { return hessianFile; }
        std::string                 getConfigFile(void) { return configFilePath; }
        std::string                 getClockModelName(void) { return clockModelName; }
        const std::vector<int>&     getClockGroups(void) { return clockGroups; }
        int                         getModelNStates(void) {
                                        if(sequenceFile.empty() == false)
                                            return datatypeProvided ? (seqDataType == "aa" ? 20 : 4) : 4;
                                        if(nstatesProvided)
                                            return nStates;
                                        return datatypeProvided ? (seqDataType == "aa" ? 20 : 4) : nStates;
                                    }
        std::string                 getSequenceFile(void) { return sequenceFile; }
        std::string                 getPartitionFile(void) { return partitionFile; }
        int                         getNumCats(void) { return numCats; }
        std::string                 getSubstModel(void) { return substModel; }
        std::string                 getFreqMode(void) { return freqMode; }
        bool                        getUseInvariant(void) { return useInvariant; }
        const double*               getRgeneGamma(void) { return rgeneGamma; }
        const double*               getSigma2Gamma(void) { return sigma2Gamma; }
        Sigma2Param                 getSigma2Param(void) { return sigma2Param; }
        unsigned long               getPncpTuningGens(void) { return pncpTuningGens; }
        bool                        clockOrCtmcConfigured(void) { return clockModelName != "ucln" || substModel != "gtr" || datatypeProvided || nstatesProvided || partitionFile.empty() == false || useInvariant || numCats != 4 || freqMode != "model" || clockGroups.empty() == false || sigma2Param != Sigma2Param::PNCP; }
        void                        print(void);
        void                        printHelp(void);

    private:
                                    UserSettings(void) {}
                                    UserSettings(const UserSettings& u);
        UserSettings&               operator=(const UserSettings&);
        void                        parsePriorInto(const std::string& spec, Probability::PriorFamily& family, double& p1, double& p2);
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
        std::vector<std::string>    psiTypeNames;
        std::map<std::string, Probability::PriorSpec> psiPriorByName;
        std::map<std::string, std::vector<double>> psiTimesByName;
        std::map<std::string, RateMode> psiModeByName;
        std::vector<int>            lambdaGroups, muGroups, psiGroups;
        std::map<std::string, std::vector<int>> psiGroupsByName;
        std::map<int, Probability::PriorSpec> lambdaGroupPrior, muGroupPrior, psiGroupPrior;
        std::map<std::string, std::map<int, Probability::PriorSpec>> psiGroupPriorByName;
        Model                       model;
        double                      rho;
        unsigned int                seed;
        bool                        seedSet;
        unsigned long               chainLength;
        int                         numCoupledChains;
        int                         numRuns;
        bool                        resume;
        bool                        arLog;
        bool                        autoChainLength;
        double                      burninFraction;
        double                      rhatThreshold;
        double                      essThreshold;
        bool                        rhatThresholdSet = false;
        bool                        essThresholdSet = false;
        unsigned long               maxGen;
        int                         numCores;
        int                         thinning;
        std::string                 hessianFile;
        std::string                 configFilePath;
        std::string                 invocation;
        std::string                 clockModelName;
        Sigma2Param                 sigma2Param = Sigma2Param::PNCP;
        unsigned long               pncpTuningGens = 50000;
        std::vector<int>            clockGroups;
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
