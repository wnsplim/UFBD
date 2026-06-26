#ifndef ChainRunner_hpp
#define ChainRunner_hpp

#include <string>
#include <vector>

class ChainRunner {

    public:
        virtual                                 ~ChainRunner(void) {}
        virtual void                            init(void) = 0;
        virtual void                            advance(unsigned long nGens) = 0;
        virtual void                            finalize(void) = 0;
        virtual void                            setOutputPaths(const std::string& po, const std::string& to) = 0;
        virtual const std::vector<std::vector<double>>& traceColumns(void) = 0;
        virtual const std::vector<std::string>& traceNames(void) = 0;
};

#endif
