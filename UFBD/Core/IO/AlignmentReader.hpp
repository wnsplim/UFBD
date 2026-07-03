#ifndef AlignmentReader_hpp
#define AlignmentReader_hpp

#include <string>
#include <vector>

class AlignmentReader {

    public:
                                        AlignmentReader(const std::string& alignmentFile, const std::string& partitionFile, int numStates);
        int                             getNumPartitions(void) const { return (int)partitionPatterns.size(); }
        const std::vector<int>&         getPartitionGroups(void) const { return partitionGroup; }
        const std::vector<std::string>& getTaxa(void) const { return taxa; }
        const std::vector<std::vector<int>>& getPatterns(int p) const { return partitionPatterns[p]; }
        const std::vector<int>&         getWeights(int p) const { return partitionWeights[p]; }

    private:
        void                            readAlignment(const std::string& file);
        void                            readPhylip(std::ifstream& f);
        void                            readFasta(std::ifstream& f);
        void                            readNexus(std::ifstream& f);
        void                            readPartitions(const std::string& file);
        void                            parseClockPartition(const std::vector<std::string>& tok);
        void                            compress(void);
        int                             charBitmask(char c) const;

        int                             numStates;
        std::vector<std::string>                    taxa;
        std::vector<std::vector<int>>               matrix;
        std::vector<std::vector<int>>               partitionSites;
        std::vector<std::string>                    partitionNames;
        std::vector<int>                            partitionGroup;
        std::vector<std::vector<std::vector<int>>>  partitionPatterns;
        std::vector<std::vector<int>>               partitionWeights;
};

#endif
