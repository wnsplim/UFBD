#ifndef UserSettings_hpp
#define UserSettings_hpp

#include <string>


enum class Conditioning { CROWN, ORIGIN };
enum class Model { FBD, HEA14, UFBD };
enum class ConditionAgePriorFamily { EXP, GAMMA, LOGNORMAL, UNIFORM };

class UserSettings {

    public:
        static UserSettings&        userSettings(void)
                                        {
                                        static UserSettings us;
                                        return us;
                                        }
        void                        initializeSettings(int argc, const char* argv[]);
        unsigned long               getChainLength(void) { checkSettings(); return chainLength; }
        int                         getNumChains(void) { return numChains; }
        int                         getNumThreads(void) { return numThreads; }
        std::string                 getTreeOutput(void) { checkSettings(); return treeOut; }
        std::string                 getParamOutput(void) { checkSettings(); return parametersOut; }
        std::string                 getTreeFile(void) { checkSettings(); return treeFile; }
        std::string                 getCladesFile(void) { checkSettings(); return cladesFile; }
        std::string                 getFossilFile(void) { checkSettings(); return fossilFile; }
        Conditioning                getConditioning(void) { checkSettings(); return conditioning; }
        ConditionAgePriorFamily          getConditionAgePrior(void) { checkSettings(); return conditionAgePrior; }
        double                      getConditionAgePriorP1(void) { checkSettings(); return conditionAgePriorP1; }
        double                      getConditionAgePriorP2(void) { checkSettings(); return conditionAgePriorP2; }
        bool                        getConditionAgePriorSet(void) { checkSettings(); return conditionAgePriorSet; }
        Model                       getModel(void) { checkSettings(); return model; }
        double                      getRho(void) { checkSettings(); return rho; }
        unsigned int                getSeed(void) { checkSettings(); return seed; }
        bool                        getSeedSet(void) { checkSettings(); return seedSet; }
        int                         getPrintFrequency(void) { checkSettings(); return printFrequency; }
        int                         getSampleFrequency(void) { checkSettings(); return sampleFrequency; }
        void                        print(void);
        void                        printHelp(void);

    private:
                                    UserSettings(void) { settingsInitialized = false; }
                                    UserSettings(const UserSettings& u);
        UserSettings&               operator=(const UserSettings&);
        void                        checkSettings(void);
        void                        parseConditionAgePrior(const std::string& spec);
        std::string                 executablePath;
        std::string                 treeOut;
        std::string                 parametersOut;
        std::string                 treeFile;
        std::string                 cladesFile;
        std::string                 fossilFile;
        Conditioning                conditioning;
        bool                        conditioningSet;
        ConditionAgePriorFamily          conditionAgePrior;
        bool                        conditionAgePriorSet;
        double                      conditionAgePriorP1;
        double                      conditionAgePriorP2;
        Model                       model;
        double                      rho;
        unsigned int                seed;
        bool                        seedSet;
        unsigned long               chainLength;
        int                         numChains;
        int                         numThreads;
        int                         printFrequency;
        int                         sampleFrequency;
        bool                        settingsInitialized;
};

#endif
