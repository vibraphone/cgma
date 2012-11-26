
#include "CubitProcess.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <sstream>

#ifdef WIN32
#include <windows.h>
#include <io.h>
#include <direct.h>
#define access _access
#define getcwd _getcwd
#define X_OK 0
#define PATH_MAX _MAX_PATH
const char path_separator = '\\';
#else
#include <unistd.h>
#include <signal.h>
#include <sys/wait.h>
#include <limits.h>
const char path_separator = '/';
#endif

static std::string quote_string(const std::string& str)
{
  std::string ret = "\"";
  ret += str;
  ret += "\"";
  return ret;
}

static std::string join_args(const std::vector<std::string>& strings)
{
  std::string joined_strings;
  for(unsigned int i=0; i<strings.size(); i++)
  {
    if(i!=0)
      joined_strings += " ";
    
    bool need_quote = strings[i].find(' ') != std::string::npos;
    if(need_quote)
      joined_strings += "\"";
    joined_strings += strings[i];
    if(need_quote)
      joined_strings += "\"";
  }
  return joined_strings;
}

// given a relative or absolute path, make it absolute
static std::string make_path_absolute(const std::string& p)
{
  std::string ret;
  if(!p.empty())
  {
    if(p[0] == path_separator)
    {
      ret = p;
    }
    else
    {
      // if any '/' character is in it, its relative to current directory
      char tmp[PATH_MAX+1];
      ret = getcwd(tmp, PATH_MAX);
      ret += path_separator;
      ret += p;
    }
  }
  return ret;
}

#ifndef WIN32
static sigset_t oldsig;
#endif

PidType CubitProcess::start(const std::string& app, const std::vector<std::string>& args, bool hide)
{
#ifndef WIN32
  std::vector<const char*> c_args(args.size()+2);
  int idx = 0;
  c_args[idx++] = app.c_str();
  for(unsigned int i=0; i<args.size(); i++)
  {
    c_args[idx++] = args[i].c_str();
  }
  c_args[idx++] = NULL;
  
  std::string app_real = find_executable(app);

  // temporarily block delivery of child signals
  // TODO: does this overwrite currently set signals for the process?
  sigset_t newsig;
  sigemptyset(&newsig);
  sigaddset(&newsig, SIGCHLD);
  sigprocmask(SIG_BLOCK, &newsig, &oldsig);

  pid_t pid = fork();
  if(pid < 0)
    return 0;

  if(pid == 0)
  {
    execv(app_real.c_str(), const_cast<char**>(&c_args[0]));
    perror(app_real.c_str());
    exit(EXIT_FAILURE);
  }

  return pid;
#else

  STARTUPINFO si;
  PROCESS_INFORMATION pi;
  
  ZeroMemory( &si, sizeof(si) );
  si.cb = sizeof(si);
  ZeroMemory( &pi, sizeof(pi) );
  
  // hide child window
  if(hide)
  {
    si.dwFlags |= STARTF_USESHOWWINDOW;
    si.wShowWindow = SW_HIDE;
  }

  std::string real_app = find_executable(app);
  std::replace(real_app.begin(), real_app.end(), '/', '\\');
  
  std::string call = quote_string(app);
  std::string joined_args = join_args(args);
  call += " ";
  call += joined_args;


    // Start the child process. 
  if( CreateProcessA( real_app.c_str(),  // path to cubit executable 
                      (char*)call.c_str(), // Command line. 
                      NULL,      // Process handle not inheritable. 
                      NULL,      // Thread handle not inheritable. 
                      TRUE,     // Set handle inheritance to TRUE.
                      0,         // No creation flags. 
                      NULL,      // Use parent's environment block. 
                      NULL,      // Use parent's starting directory. 
                      &si,       // Pointer to STARTUPINFO structure.
                      &pi )      // Pointer to PROCESS_INFORMATION structure.
      )
  {
    return pi;
  }
  pi.hProcess = 0;
  return pi;
#endif
}

int CubitProcess::wait(PidType pid)
{
#ifndef WIN32
  int status;
  int result = -1;
  
  waitpid(pid, &status, 0);
  if(WIFEXITED(status))
    result = WEXITSTATUS(status);
  else if (WIFSIGNALED(status))
    result = 128 + WTERMSIG(status);
  
  sigprocmask(SIG_SETMASK, &oldsig, NULL);

  return result;
#else
  
  // Wait until child process exits.
  WaitForSingleObject( pid.hProcess, INFINITE );
    
  // Get its exit code
  DWORD exit_code = -1;
  GetExitCodeProcess(pid.hProcess, &exit_code);
    
  // Close process and thread handles. 
  CloseHandle( pid.hProcess );
  CloseHandle( pid.hThread );
  return exit_code;

#endif
}


int CubitProcess::execute(const std::string& app, const std::vector<std::string>& args, bool hide)
{
  PidType pid = CubitProcess::start(app, args, hide);
#if WIN32
  if(pid.hProcess == 0)
    return -1;
#else
  if(pid == 0)
    return -1;
#endif
  return CubitProcess::wait(pid);
}
  
std::string CubitProcess::find_executable(const std::string& app)
{
  std::string app_real = app;

#ifdef WIN32 
  if(app_real.rfind(".exe") == std::string::npos)
	  app_real += ".exe";
  const char path_delimiter = ';';
#else
  const char path_delimiter = ':';
#endif

#ifdef WIN32
  if(app_real[0] != path_separator && 
	  (app_real.size() > 1 && app_real[1] != ':'))
#else
  if(app_real[0] != path_separator)
#endif
  {
    if(app_real.find(path_separator) != std::string::npos)
    {
      // if any '/' character is in it, its relative to current directory
      app_real = make_path_absolute(app_real);
    }
    else
    {
      // search PATH env. for this executable (note PATH may have relative directories)
      std::vector<std::string> paths;
      std::stringstream ss(getenv("PATH"));
      std::string item;
      while(std::getline(ss, item, path_delimiter))
      {
        paths.push_back(item);
      }

      for(size_t i=0; i<paths.size(); i++)
      {
        std::string p = paths[i];
        p += "/";
        p += app;
        std::string abs_p = make_path_absolute(p);
        if(0 == access(abs_p.c_str(), X_OK))
        {
          app_real = abs_p;
          break;
        }
      }
    }
  }
  return app_real;
}
