#include <windows.h>

#include <cstdio>
#include <cstdlib>
#include <cwchar>
#include <string>

static std::wstring quote(const std::wstring& a) {   
  if (!a.empty() && a.find_first_of(L" \t\"") == std::wstring::npos) return a;
  std::wstring r = L"\"";
  for (size_t i = 0;; ++i) {
    size_t bs = 0;
    while (i < a.size() && a[i] == L'\\') ++bs, ++i;
    if (i == a.size()) {
      r.append(bs * 2, L'\\');
      break;
    }
    if (a[i] == L'"') {
      r.append(bs * 2 + 1, L'\\');
      r += L'"';
    } else {
      r.append(bs, L'\\');
      r += a[i];
    }
  }
  return r + L"\"";
}

int wmain(int argc, wchar_t** argv) {
  if (argc < 3) {
    std::fwprintf(stderr, L"usage: memcap <limit_MB> <program> [args...]\n");
    return 2;
  }
  const unsigned long long limit_mb = _wcstoui64(argv[1], nullptr, 10);
  if (limit_mb == 0) {
    std::fwprintf(stderr, L"memcap: bad limit %ls\n", argv[1]);
    return 2;
  }
  std::wstring cmd;
  for (int i = 2; i < argc; ++i) {
    if (i > 2) cmd += L' ';
    cmd += quote(argv[i]);
  }
  HANDLE job = CreateJobObjectW(nullptr, nullptr);
  HANDLE port = CreateIoCompletionPort(INVALID_HANDLE_VALUE, nullptr, 0, 1);
  if (!job || !port) {
    std::fwprintf(stderr, L"memcap: cannot create the job (error %lu)\n", GetLastError());
    return 127;
  }
  JOBOBJECT_ASSOCIATE_COMPLETION_PORT acp{};
  acp.CompletionKey = job;
  acp.CompletionPort = port;
  SetInformationJobObject(job, JobObjectAssociateCompletionPortInformation, &acp, sizeof acp);

  const char* lowpri_env = std::getenv("MEMCAP_LOWPRI");
  const bool lowpri = lowpri_env && lowpri_env[0] == '1';
  const char* aff_env = std::getenv("MEMCAP_AFFINITY");
  const unsigned long long affinity = aff_env ? std::strtoull(aff_env, nullptr, 16) : 0;
  const char* qos_env = std::getenv("MEMCAP_HIGHQOS");
  const bool highqos = qos_env && qos_env[0] == '1';
  JOBOBJECT_EXTENDED_LIMIT_INFORMATION li{};
  li.BasicLimitInformation.LimitFlags = JOB_OBJECT_LIMIT_PROCESS_MEMORY | JOB_OBJECT_LIMIT_JOB_MEMORY |
                                        JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE |
                                        (lowpri ? JOB_OBJECT_LIMIT_PRIORITY_CLASS : 0) |
                                        (affinity ? JOB_OBJECT_LIMIT_AFFINITY : 0);
  li.BasicLimitInformation.PriorityClass = BELOW_NORMAL_PRIORITY_CLASS;
  li.BasicLimitInformation.Affinity = static_cast<ULONG_PTR>(affinity);
  if (affinity || highqos) std::fprintf(stderr, "memcap: affinity=0x%llx highqos=%d\n", affinity, highqos ? 1 : 0);
  li.ProcessMemoryLimit = li.JobMemoryLimit = static_cast<SIZE_T>(limit_mb << 20);
  if (!SetInformationJobObject(job, JobObjectExtendedLimitInformation, &li, sizeof li)) {
    std::fwprintf(stderr, L"memcap: cannot set the limit (error %lu)\n", GetLastError());
    return 127;
  }
  STARTUPINFOW si{};
  si.cb = sizeof si;
  si.dwFlags = STARTF_USESTDHANDLES;
  si.hStdInput = GetStdHandle(STD_INPUT_HANDLE);
  si.hStdOutput = GetStdHandle(STD_OUTPUT_HANDLE);
  si.hStdError = GetStdHandle(STD_ERROR_HANDLE);
  for (HANDLE h : {si.hStdInput, si.hStdOutput, si.hStdError})
    if (h && h != INVALID_HANDLE_VALUE) SetHandleInformation(h, HANDLE_FLAG_INHERIT, HANDLE_FLAG_INHERIT);

  const bool has_path = std::wcspbrk(argv[2], L"/\\") != nullptr;
  PROCESS_INFORMATION pi{};
  if (!CreateProcessW(has_path ? argv[2] : nullptr, cmd.data(), nullptr, nullptr, TRUE,
                      CREATE_SUSPENDED | (lowpri ? BELOW_NORMAL_PRIORITY_CLASS : 0), nullptr, nullptr, &si, &pi)) {
    std::fwprintf(stderr, L"memcap: cannot start %ls (error %lu)\n", argv[2], GetLastError());
    return 127;
  }
  if (!AssignProcessToJobObject(job, pi.hProcess)) {
    std::fwprintf(stderr, L"memcap: cannot assign the job (error %lu)\n", GetLastError());
    TerminateProcess(pi.hProcess, 127);
    return 127;
  }
  if (highqos) {
    PROCESS_POWER_THROTTLING_STATE pt{};
    pt.Version = PROCESS_POWER_THROTTLING_CURRENT_VERSION;
    pt.ControlMask = PROCESS_POWER_THROTTLING_EXECUTION_SPEED;
    pt.StateMask = 0;   
    if (!SetProcessInformation(pi.hProcess, ProcessPowerThrottling, &pt, sizeof pt))
      std::fprintf(stderr, "memcap: cannot turn off power throttling (error %lu)\n", GetLastError());
  }
  ResumeThread(pi.hThread);
  CloseHandle(pi.hThread);
  WaitForSingleObject(pi.hProcess, INFINITE);
  DWORD code = 0;
  GetExitCodeProcess(pi.hProcess, &code);
  bool hit = false;
  DWORD msg;
  ULONG_PTR key;
  LPOVERLAPPED ov;
  while (GetQueuedCompletionStatus(port, &msg, &key, &ov, 0))
    if (msg == JOB_OBJECT_MSG_PROCESS_MEMORY_LIMIT || msg == JOB_OBJECT_MSG_JOB_MEMORY_LIMIT) hit = true;
  QueryInformationJobObject(job, JobObjectExtendedLimitInformation, &li, sizeof li, nullptr);
  std::fprintf(stderr, "memcap: limit_mb=%llu peak_mb=%llu exit=%lu\n", limit_mb,
               static_cast<unsigned long long>(li.PeakProcessMemoryUsed >> 20), code);
  if (hit) std::fprintf(stderr, "memcap: memory limit reached\n");
  return static_cast<int>(code);
}
