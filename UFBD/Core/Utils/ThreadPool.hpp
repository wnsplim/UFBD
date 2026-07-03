#ifndef ThreadPool_hpp
#define ThreadPool_hpp

#include <atomic>
#include <condition_variable>
#include <functional>
#include <iostream>
#include <mutex>
#include <queue>
#include <thread>
#include <vector>

enum { OP_GAMMA = 0, OP_CTMC, OP_FBD, OP_CLOCK, OP_NUM };

class ThreadPool {
    public:
                                        ThreadPool(int numThreads);
                                       ~ThreadPool();
    void                                enqueue(std::function<void()> task);
    int                                 size(void) { return (int)threads.size(); }
    void                                parallelFor(int opId, int n, const std::function<void(int, int)>& body, int maxThreads = 1000000);
    void                                setChainCap(int c) { chainCap = (c < 1) ? 1 : c; }
    static ThreadPool&                  shared(void);
    static ThreadPool&                  current(void);
    static void                         setCurrent(ThreadPool* p) { tlsCurrent = p; }
    private:
    static thread_local ThreadPool*     tlsCurrent;
    std::vector<std::thread>            threads;
    std::queue<std::function<void()>>   tasks;
    std::mutex                          queue_mutex;
    std::condition_variable             cv;
    bool                                stop;
    int                                 chainCap;
    std::atomic<double>                 opCost[OP_NUM];
    std::atomic<long>                   opCalls[OP_NUM];
};

#endif
