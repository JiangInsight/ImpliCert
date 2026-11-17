#define KCDS_SOLVE_MAX_IT 1
#define main covered_solver_main
#include "main.cpp"
#undef main

int main(int argc,char** argv) {
  if (argc>1) return covered_solver_main(argc,argv);
  guard_audit::boxes.assign(1,{});
  if (clampLam(301,0)!=300 || clampLam(-302,0)!=-300 || clampLam(300,0)!=300) return 1;
  auto b=guard_audit::boxes.front();
  if (b.writes!=3 || b.upper!=2 || b.lower!=1 || b.clipped!=2) return 1;
  if (guard_audit::choose_bound(5,3)!=10 || guard_audit::choose_bound(3,4)!=0 ||
      guard_audit::choose_bound(1000000,5)!=guard_audit::too_large) return 1;
  try { guard_audit::exact_count(guard_audit::exact_limit+1); return 1; }
  catch (const std::runtime_error&) {}
  std::puts("box counters, combination bound and integer-range negative control passed");
  return 0;
}
