//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, OpenSees/Xara Developers
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
//===----------------------------------------------------------------------===//
//
#pragma once
#include <array>
#include <cassert>
#include <chrono>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <string_view>
#include <type_traits>

#include <utility/magic_enum.hpp>

class TimerBase {
  virtual void print(std::ostream&) const {};
};

template<typename Step>
requires std::is_enum_v<Step>
class Timer : public TimerBase {
private:
  using Clock    = std::chrono::steady_clock;
  using Duration = Clock::duration;

  inline static constexpr auto names_ =
      magic_enum::enum_names<Step>();

  static constexpr std::size_t step_count_ = names_.size();

  template<Step S>
  static consteval std::size_t index()
  {
    return magic_enum::enum_index<S>();
  }

  static consteval std::size_t name_width()
  {
    std::size_t width = std::string_view{"Step"}.size();

    for (const auto name : names_) {
      if (name.size() > width)
        width = name.size();
    }

    return width;
  }

public:
  static_assert(step_count_ > 0, "Step enum contains no reflected values.");

  template<Step S>
  void start() noexcept
  {
    constexpr std::size_t i = index<S>();

    assert(!active_[i] && "Timer step is already active.");

    start_times_[i] = Clock::now();
    active_[i] = true;
  }

  template<Step S>
  void stop() noexcept
  {
    constexpr std::size_t i = index<S>();

    assert(active_[i] && "Timer step is not active.");

    totals_[i] += Clock::now() - start_times_[i];
    ++counts_[i];
    active_[i] = false;
  }

  void reset() noexcept
  {
    totals_.fill(Duration::zero());
    counts_.fill(0);
    active_.fill(false);
  }

  void print(std::ostream& output) const override
  {

    constexpr int step_width =
        static_cast<int>(name_width());

    constexpr int time_width    = 14;
    constexpr int percent_width = 5;
    constexpr int total_width   =
        time_width + 2 + percent_width + 2;  // "<time> (<percent>%)"

    constexpr std::size_t line_width =
        step_width + 2 + time_width + 2 + total_width;

    const auto flags     = output.flags();
    const auto precision = output.precision();
    const auto fill      = output.fill();

    double grand_total = 0.0;

    for (const Duration duration : totals_) {
      grand_total +=
          std::chrono::duration<double, std::milli>(duration).count();
    }

    output << std::left
           << std::setw(step_width) << "Step"
           << "  "
           << std::right
           << std::setw(time_width) << "Avg. [ms]"
           << "  "
           << std::setw(total_width) << "Total [ms] (%)"
           << '\n';

    for (std::size_t i = 0; i < line_width; ++i)
      output.put('-');

    output.put('\n');

    for (std::size_t i = 0; i < step_count_; ++i) {
      const double total =
          std::chrono::duration<double, std::milli>(
              totals_[i]
          ).count();

      const double average =
          counts_[i] == 0
              ? 0.0
              : total / static_cast<double>(counts_[i]);

      const double percentage =
          grand_total == 0.0
              ? 0.0
              : 100.0 * total / grand_total;

      output << std::left
             << std::setw(step_width) << names_[i]
             << "  "
             << std::right
             << std::fixed
             << std::setprecision(0)
             << std::setw(time_width) << average
             << "  "
             << std::setw(time_width) << total
             << " ("
             << std::setprecision(1)
             << std::setw(percent_width) << percentage
             << "%)\n";
    }

    // Total row
    for (std::size_t i = 0; i < line_width; ++i)
      output.put('-');

    output.put('\n');

    output << std::left
           << std::setw(step_width) << "Total"
           << "  "
           << std::right
           << std::setw(time_width) << ""
           << "  "
           << std::fixed
           << std::setprecision(0)
           << std::setw(time_width) << grand_total
           << " ("
           << std::setprecision(1)
           << std::setw(percent_width) << (grand_total > 0.0 ? 100.0 : 0.0)
           << "%)\n";

    output.flags(flags);
    output.precision(precision);
    output.fill(fill);
  }

private:
  std::array<Clock::time_point, step_count_> start_times_{};
  std::array<Duration, step_count_> totals_{};
  std::array<std::size_t, step_count_> counts_{};
  std::array<bool, step_count_> active_{};
};
