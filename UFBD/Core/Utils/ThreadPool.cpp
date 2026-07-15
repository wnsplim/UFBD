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
    dispatchCost.store(-1.0, std::memory_order_relaxed);
    for(int i = 0; i < OP_NUM; i++){
        opCost[i].store(-1.0, std::memory_order_relaxed);
        opCalls[i].store(0, std::memory_order_relaxed);
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

double ThreadPool::measureDispatch(int nt){
    if(nt <= 1)
        return 0.0;
    std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
    int remaining = nt;
    std::mutex m;
    std::condition_variable doneCv;
    for(int k = 0; k < nt; k++)
        enqueue([&remaining, &m, &doneCv](){
            std::lock_guard<std::mutex> lk(m);
            if(--remaining == 0)
                doneCv.notify_one();
        });
    std::unique_lock<std::mutex> lk(m);
    doneCv.wait(lk, [&remaining]{ return remaining == 0; });
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

void ThreadPool::parallelFor(int opId, int n, const std::function<void(int, int)>& body, int maxThreads){
    if(n <= 0)
        return;
    int nt = (int)threads.size();
    if(nt > maxThreads) nt = maxThreads;
    if(nt > chainCap) nt = chainCap;
    long calls = opCalls[opId].fetch_add(1, std::memory_order_relaxed) + 1;
    double cost = opCost[opId].load(std::memory_order_relaxed);
    if(nt <= 1 || cost < 0.0 || (calls % 256) == 0){
        std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
        body(0, n);
        double per = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() / (double)n;
        opCost[opId].store((cost < 0.0) ? per : 0.8 * cost + 0.2 * per, std::memory_order_relaxed);
        return;
    }
    double disp = dispatchCost.load(std::memory_order_relaxed);
    if(disp < 0.0 || (calls % 4099) == 0){
        double d = measureDispatch(nt);
        disp = (disp < 0.0) ? d : 0.7 * disp + 0.3 * d;
        dispatchCost.store(disp, std::memory_order_relaxed);
    }
    if((double)n * cost < 2.0 * disp){
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

void ThreadPool::parallelTasks(const std::vector<std::function<void()>>& tasks){
    int n = (int)tasks.size();
    if(n <= 0)
        return;
    if((int)threads.size() <= 1){
        for(int k = 0; k < n; k++)
            tasks[k]();
        return;
    }
    int remaining = n;
    std::mutex m;
    std::condition_variable doneCv;
    for(int k = 0; k < n; k++){
        enqueue([&tasks, k, &remaining, &m, &doneCv](){
            tasks[k]();
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

thread_local ThreadPool* ThreadPool::tlsCurrent = nullptr;

ThreadPool& ThreadPool::current(void){
    return tlsCurrent ? *tlsCurrent : shared();
}
