#pragma once

#include <algorithm>
#include <chrono>
#include <functional>
#include <vector>

namespace rnla {

// Median wall-clock seconds over `reps` runs of f, after `warmup` untimed runs.
template <typename F>
double time_median(F&& f, int reps = 7, int warmup = 2) {
  for (int i = 0; i < warmup; ++i) f();

  std::vector<double> times;
  times.reserve(reps);
  for (int i = 0; i < reps; ++i) {
    const auto t0 = std::chrono::steady_clock::now();
    f();
    const auto t1 = std::chrono::steady_clock::now();
    times.push_back(std::chrono::duration<double>(t1 - t0).count());
  }
  std::sort(times.begin(), times.end());
  return times[times.size() / 2];
}

}  // namespace rnla