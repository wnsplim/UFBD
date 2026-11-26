#include "ChildJobs.hpp"
#include "Msg.hpp"

#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <iostream>
#include <map>
#include <spawn.h>
#include <sys/wait.h>
#include <unistd.h>

extern char** environ;

void runChildJobs(const std::vector<ChildJob>& jobs, int width, const std::string& what){
    setenv("UFBD_PARENT", std::to_string(getpid()).c_str(), 1);
    std::map<pid_t,size_t> running;
    size_t next = 0;
    std::string failed;
    while(next < jobs.size() || running.empty() == false){
        while(failed.empty() && next < jobs.size() && (int)running.size() < width){
            std::vector<char*> argv;
            for(const std::string& a : jobs[next].args)
                argv.push_back(const_cast<char*>(a.c_str()));
            argv.push_back(nullptr);
            posix_spawn_file_actions_t fa;
            posix_spawn_file_actions_init(&fa);
            if(jobs[next].stdoutToConsole)
                posix_spawn_file_actions_addopen(&fa, 1, jobs[next].console.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
            else
                posix_spawn_file_actions_addopen(&fa, 1, "/dev/null", O_WRONLY, 0);
            posix_spawn_file_actions_adddup2(&fa, 1, 2);
            pid_t pid;
            int rc = posix_spawnp(&pid, argv[0], &fa, nullptr, argv.data(), environ);
            posix_spawn_file_actions_destroy(&fa);
            if(rc != 0)
                Msg::error(what + ": could not start " + jobs[next].label + ": " + std::strerror(rc));
            std::cout << "started  " << jobs[next].label << "\n" << std::flush;
            running[pid] = next++;
        }
        if(running.empty())
            break;
        int status = 0;
        pid_t pid = waitpid(-1, &status, 0);
        size_t i = running[pid];
        running.erase(pid);
        if(WIFEXITED(status) && WEXITSTATUS(status) == 0)
            std::cout << "finished " << jobs[i].label << "\n" << std::flush;
        else if(failed.empty())
            failed = jobs[i].label + " failed (" + (WIFEXITED(status) ? "exit status " + std::to_string(WEXITSTATUS(status)) : "signal " + std::to_string(WTERMSIG(status))) + "); its console is " + jobs[i].console;
    }
    if(failed.empty() == false)
        Msg::error(what + ": " + failed + ".");
}
