#ifndef AdaptiveOperatorSampler_hpp
#define AdaptiveOperatorSampler_hpp

#include <vector>

class AdaptiveOperatorSampler {

    public:
                            AdaptiveOperatorSampler(void) = delete;
                            AdaptiveOperatorSampler(int nMoves);
        int                 pick(double u);
        void                reward(int move, double r);
        void                setWeight(int m, double w) { weight[m] = w; }
        void                setActive(int m, bool a) { active[m] = a ? 1 : 0; }
        void                setPinned(int m, double w) { pinned[m] = 1; weight[m] = w; }
        void                setAdapting(bool a) { adapting = a; }
        double              weightOf(int m) { return weight[m]; }

    private:
        void                adapt(void);
        int                 n;
        bool                adapting;
        std::vector<double> weight;
        std::vector<double> rewSum;
        std::vector<double> rewCnt;
        std::vector<char>   active;
        std::vector<char>   pinned;
        double              adaptN;
        long                picks;
};

#endif
