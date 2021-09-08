#ifdef _MSC_VER
#pragma once
#endif 

#ifndef PROGRESS_H
#define PROGRESS_H

#include <iostream>
#include <string>
#include <chrono>
#include <sstream>

class ProgressBar {
public:
    ProgressBar() {}
    ProgressBar(int total) {
        m_step = 0;
        m_total = total;
        start = std::chrono::system_clock::now();
    }
    virtual ~ProgressBar() {}

    void step(int n = 1) {
        m_step += n;
        auto now = std::chrono::system_clock::now();

        const double percent = 100.0 * m_step / m_total;
        const int64_t elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count();
        const double time_msec_per_step = (double)elapsed / (double)m_step;
        const double rest_time = time_msec_per_step * (m_total - m_step);

        const int n_min = (int)(elapsed / 1000.0) / 60;
        const int n_sec = (int)(elapsed / 1000.0) % 60;
        const int r_min = (int)(rest_time / 1000.0) / 60;
        const int r_sec = (int)(rest_time / 1000.0) % 60;

        const int steps_per_sec = (int)(1000.0 / time_msec_per_step);
        std::ostringstream oss;
        oss << steps_per_sec;
        const std::string it_text = steps_per_sec < 1000 ? oss.str() : "1000+";

        const int tick = (int)(m_width * m_step / m_total);
        std::string pbar = std::string(tick, '=');
        if (tick != m_width) {
            pbar += ">";
            pbar += std::string(m_width - tick - 1, ' ');
        }

        printf("\r[%3d%%]|%s| %d/%d [%02d:%02d<%02d:%02d, %sit/s]",
               (int)percent, pbar.c_str(), m_step, m_total, n_min, n_sec, r_min, r_sec, it_text.c_str());

        if (m_step == m_total) {
            printf("\n");
        }
    }

    void finish() {
        printf("\n");
    }

private:
    const int m_width = 40;
    int m_step, m_total;
    std::chrono::system_clock::time_point start;
};

#endif  // PROGRESS_H