#pragma once
#include <atomic>


inline bool claim_absorber(char& slot) {
  static_assert(std::atomic_ref<char>::required_alignment == alignof(char));
  std::atomic_ref<char> flag(slot);
  return !flag.load(std::memory_order_relaxed) &&
         !flag.exchange(1, std::memory_order_relaxed);
}
