#ifdef _MSC_VER
#pragma once
#endif

#ifndef PROGRESS_H
#define PROGRESS_H

#include <algorithm>
#include <format>
#include <string>
#include <chrono>
#include <utility>

#include "debug.h"

namespace tinymesh {

class ProgressBar {
    using Clock = std::chrono::steady_clock;

public:
    explicit ProgressBar(int total)
        : m_total(total)
        , start(Clock::now()) {
    }

    virtual ~ProgressBar() {
    }

    template <typename... Args>
    void set_description(std::format_string<Args...> fmt, Args &&...args) {
        m_description = STR_FMT(fmt, std::forward<Args>(args)...);
    }

    void step(int n = 1) {
        if (m_total <= 0) {
            return;
        }

        m_step = std::clamp(m_step + n, 0, m_total);

        const auto now = Clock::now();
        const double elapsed = std::chrono::duration<double, std::milli>(now - start).count();

        const double percent = 100.0 * m_step / m_total;
        const double time_msec_per_step = elapsed / std::max(1, m_step);
        const double rest_time = time_msec_per_step * (m_total - m_step);

        const int n_min = static_cast<int>(elapsed / 1000.0) / 60;
        const int n_sec = static_cast<int>(elapsed / 1000.0) % 60;
        const int r_min = static_cast<int>(rest_time / 1000.0) / 60;
        const int r_sec = static_cast<int>(rest_time / 1000.0) % 60;

        const double steps_per_sec = 1000.0 / time_msec_per_step;

        std::string it_text;
        if (steps_per_sec >= 1000.0) {
            it_text = "1000+it/s";
        } else if (steps_per_sec >= 1.0) {
            it_text = std::format("{}it/s", static_cast<int>(steps_per_sec));
        } else {
            it_text = std::format("{:.2f}s/it", time_msec_per_step / 1000.0);
        }

        const int tick = (int)(m_width * m_step / m_total);
        std::string pbar = std::string(tick, '=');
        if (tick != m_width) {
            pbar += ">";
            pbar += std::string(m_width - tick - 1, ' ');
        }

        if (m_step == m_total || m_step % std::max(1, (m_total / 1000)) == 0) {
            const std::string description = m_description.empty() ? "" : m_description + " ";
            Print("\r{}[{:3d}%]|{}| {}/{} [{:02d}:{:02d}<{:02d}:{:02d}, {}]", description, static_cast<int>(percent),
                  pbar, m_step, m_total, n_min, n_sec, r_min, r_sec, it_text);
        }

        if (m_step == m_total) {
            Print("\n");
        }
    }

    void finish() {
        if (m_step != m_total) {
            step(m_total - m_step);
        }
    }

private:
    const int m_width = 40;
    int m_step = 0;
    int m_total;
    std::string m_description = "";
    Clock::time_point start;
};

}  // namespace tinymesh

#endif  // PROGRESS_H
