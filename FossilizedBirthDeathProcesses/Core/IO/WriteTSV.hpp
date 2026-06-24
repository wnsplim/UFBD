#ifndef WriteTSV_hpp
#define WriteTSV_hpp

#include <fstream>
#include <string>
#include <vector>

class WriteTSV{
    public:
                                    WriteTSV(void);
                                    WriteTSV(std::string filepath, bool overwrite);
                                   ~WriteTSV(void);
        void                        addColumnNamesTSV(std::vector<std::string> cn);
        void                        addFilepath(std::string fp, bool overwrite);
        void                        appendDataTSV(std::vector<double> data);
        void                        appendDataTSV(std::string data);
        void                        closeTSV(void);

    private:
        std::fstream                fout;
        int                         numCols;
        int                         numRows;
        std::string                 filepath;
};

#endif
