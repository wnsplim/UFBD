#ifndef UserSettings_hpp
#define UserSettings_hpp

#include <string>


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
        int                         getPrintFrequency(void) { checkSettings(); return printFrequency; }
        int                         getSampleFrequency(void) { checkSettings(); return sampleFrequency; }
        void                        print(void);
        void                        printHelp(void);

    private:
                                    UserSettings(void) { settingsInitialized = false; }
                                    UserSettings(const UserSettings& u);
        UserSettings&               operator=(const UserSettings&);
        void                        checkSettings(void);
        std::string                 executablePath;
        std::string                 treeOut;
        std::string                 parametersOut;
        std::string                 treeFile;
        std::string                 cladesFile;
        std::string                 fossilFile;
        unsigned long               chainLength;
        int                         numChains;
        int                         numThreads;
        int                         printFrequency;
        int                         sampleFrequency;
        bool                        settingsInitialized;
};

#endif
