#pragma once
#include "export_snapshot.hpp"
#include <process.h>


static void checkpoint_barrier(const Solver& s, int batch, const Timer& clock,
                               double engine, int rounds, size_t dels, size_t ins) {
  const std::string path = std::string(std::getenv("TRACE_DIR")) + "/" + std::to_string(batch) + ".cert";
  Timer exporting;


  const auto selector = std::to_string(batch + 1);
  if (_putenv_s("KCDS_CERT_FILE", path.c_str()) || _putenv_s("KCDS_CERT_BATCH", selector.c_str()))
    throw std::runtime_error("export environment");
  export_certificate_if_requested(s, batch + 1);
  const double export_s = exporting.sec();
  const std::string python = std::getenv("TRACE_PYTHON");
  auto quote = [](const std::string& x) { return "\"" + x + "\""; };
  const auto program = quote(python), checker = quote(std::getenv("TRACE_CHECKER"));
  const auto native = quote(std::getenv("TRACE_NATIVE")), checkpoint = quote(path);
  const char* args[] = {program.c_str(), "-I", checker.c_str(), native.c_str(), checkpoint.c_str(), "1/100", nullptr};
  std::fflush(stdout);
  Timer verifying;
  const intptr_t status = _spawnv(_P_WAIT, python.c_str(), args);
  const double verify_s = verifying.sec();
  if (status != 0) throw std::runtime_error("rigorous verification rejected batch " + std::to_string(batch));
  const double elapsed = clock.sec();
  std::printf("accepted batch=%d total_s=%.9f engine_s=%.9f export_s=%.9f verify_s=%.9f rounds=%d dels=%zu ins=%zu\n",
              batch, elapsed, engine, export_s, verify_s, rounds, dels, ins);
  std::fflush(stdout);
}
