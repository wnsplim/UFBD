#ifndef ChainRunner_hpp
#define ChainRunner_hpp

#include <fstream>
#include <iostream>
#include <mutex>
#include <string>
#include <vector>

#include "WriteTSV.hpp"

class RandomVariable;
class Tree;

class ChainRunner {

    public:
        virtual                                 ~ChainRunner(void) {}
        virtual Tree*                           getTree(void) { return nullptr; }
        virtual void                            setVerbose(bool) {}
        virtual void                            setLabel(int) {}
        static void                             logLine(const std::string& s){ static std::mutex m; std::lock_guard<std::mutex> lk(m); std::cout << s; }
        virtual void                            init(void) = 0;
        virtual void                            advance(unsigned long nGens) = 0;
        virtual void                            finalize(void) = 0;
        void                                    setOutputPaths(const std::string& po, const std::string& to) { paramOut = po; treeOut = to; }
        const std::vector<std::vector<double>>& traceColumns(void) { return traceCols; }
        const std::vector<std::string>&         traceNames(void) { return traceNms; }
        const std::vector<bool>&                fixedTraceMask(void) { return fixedMask; }
        const std::vector<std::vector<double>>& latentColumns(void) { return latentCols; }
        const std::vector<std::string>&         latentNames(void) { return latentNms; }
        virtual bool                            treeHasFossils(void) { return false; }
        virtual void                            printMoveDiagnostics(int rep) {}
        virtual void                            writeCheckpoint(void) {}
        virtual bool                            loadCheckpoint(void) { return false; }
        void                                    resumeOutputs(void);
        unsigned long                           currentGen(void) { return gen; }

    protected:
        static std::ifstream                    openCheckpoint(const std::string& path);
        static void                             reconcileThinning(int storedSf, int& thinning);
        static void                             requireCheckpointIntact(std::istream& is, const std::string& path);
        virtual RandomVariable*                 resumeRng(void) = 0;
        virtual std::vector<std::string>        resumeParameterNames(void) = 0;
        virtual std::vector<std::string>        resumeLatentNames(void) = 0;
        void                                    resumeLatentOutput(void);

        WriteTSV                                params;
        WriteTSV                                trees;
        WriteTSV                                latent;
        std::string                             paramOut;
        std::string                             treeOut;
        std::string                             latentOut;
        bool                                    writeTrees;
        bool                                    writeLatent = false;
        unsigned long                           gen;
        std::vector<std::vector<double>>        traceCols;
        std::vector<std::string>                traceNms;
        std::vector<bool>                       fixedMask;
        std::vector<std::vector<double>>        latentCols;
        std::vector<std::string>                latentNms;
};

#endif
