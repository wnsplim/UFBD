#include "ThreadPool.hpp"
#include "UserSettings.hpp"

#include <algorithm>
#include <chrono>
#include <condition_variable>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <mutex>
#include <queue>
#include <thread>
#include <utility>
#include <vector>

ThreadPool::ThreadPool(int numThreads) : stop(false), chainCap(numThreads){
    for(int i = 0; i < OP_NUM; i++){
        opCost[i] = -1.0;
        opCalls[i] = 0;
    }
    for (size_t i = 0; i < numThreads; ++i) {
        threads.emplace_back([this] {
            while (true) {
                std::function<void()> task;
                {
                    std::unique_lock<std::mutex> lock(queue_mutex);
                    cv.wait(lock, [this] {
                        return !tasks.empty() || stop;
                    });
                    if (stop && tasks.empty()) {
                        return;
                    }
                    task = std::move(tasks.front());
                    tasks.pop();
                }

                task();
            }
        });
    }
}

ThreadPool::~ThreadPool()
{
    {
        std::unique_lock<std::mutex> lock(queue_mutex);
        stop = true;
    }
    
    cv.notify_all();

    for (auto& thread : threads) {
        thread.join();
    }
}

void ThreadPool::enqueue(std::function<void()> task){
    {
        std::unique_lock<std::mutex> lock(queue_mutex);
        tasks.emplace(std::move(task));
    }
    cv.notify_one();
}

void ThreadPool::parallelFor(int opId, int n, const std::function<void(int, int)>& body, int maxThreads){
    if(n <= 0)
        return;
    int nt = (int)threads.size();
    if(nt > maxThreads) nt = maxThreads;
    if(nt > chainCap) nt = chainCap;
    long calls = ++opCalls[opId];
    double cost = opCost[opId];
    if(nt <= 1 || cost < 0.0 || (calls % 256) == 0){
        std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
        body(0, n);
        double per = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() / (double)n;
        opCost[opId] = (cost < 0.0) ? per : 0.8 * cost + 0.2 * per;
        return;
    }
    if((double)n * cost < 150e-6){
        body(0, n);
        return;
    }
    int grain = (int)(50e-6 / cost);
    if(grain < 1) grain = 1;
    int chunks = n / grain;
    if(chunks > nt) chunks = nt;
    if(chunks < 2){
        body(0, n);
        return;
    }
    int per = (n + chunks - 1) / chunks;
    std::vector<std::pair<int,int> > ranges;
    for(int s = 0; s < n; s += per)
        ranges.push_back(std::make_pair(s, std::min(s + per, n)));
    int remaining = (int)ranges.size();
    std::mutex m;
    std::condition_variable doneCv;
    for(size_t r = 0; r < ranges.size(); r++){
        int a = ranges[r].first, b = ranges[r].second;
        enqueue([&remaining, &m, &doneCv, &body, a, b](){
            body(a, b);
            std::lock_guard<std::mutex> lk(m);
            if(--remaining == 0)
                doneCv.notify_one();
        });
    }
    std::unique_lock<std::mutex> lk(m);
    doneCv.wait(lk, [&remaining]{ return remaining == 0; });
}

ThreadPool& ThreadPool::shared(void){
    static ThreadPool pool(UserSettings::userSettings().getNumCores());
    return pool;
}
