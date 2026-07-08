#ifndef ChainRunner_hpp
#define ChainRunner_hpp

#include <iostream>
#include <mutex>
#include <string>
#include <vector>

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
        virtual void                            setOutputPaths(const std::string& po, const std::string& to) = 0;
        virtual const std::vector<std::vector<double>>& traceColumns(void) = 0;
        virtual const std::vector<std::string>& traceNames(void) = 0;
        virtual const std::vector<std::vector<double>>& latentColumns(void) { static std::vector<std::vector<double>> e; return e; }
        virtual const std::vector<std::string>& latentNames(void) { static std::vector<std::string> e; return e; }
        virtual bool                            treeHasFossils(void) { return false; }
        virtual void                            writeCheckpoint(void) {}
        virtual bool                            loadCheckpoint(void) { return false; }
        virtual void                            resumeOutputs(void) {}
        virtual unsigned long                   currentGen(void) { return 0; }
};

#endif
