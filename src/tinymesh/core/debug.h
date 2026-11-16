#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_DEBUG_H
#define TINYMESH_DEBUG_H

#include <cstdlib>
#include <exception>
#include <format>
#include <iostream>
#include <sstream>
#include <utility>
#include <source_location>

#if defined(TINYMESH_PYTHON_MODULE)
#include <pybind11/pybind11.h>
namespace py = pybind11;
#endif  // TINYMESH_PYTHON_MODULE

// -----------------------------------------------------------------------------
// ANSI color
// -----------------------------------------------------------------------------
#define FG_RED "\033[31m"
#define FG_GRN "\033[32m"
#define FG_YEL "\033[33m"
#define FG_BLU "\033[34m"
#define FG_NC "\033[0m"

// -----------------------------------------------------------------------------
// Debug output stream
// -----------------------------------------------------------------------------
class PrintStream {
public:
    PrintStream() {
    }

    virtual ~PrintStream() {
        flush();
    }

    template <typename Type>
    PrintStream &operator<<(const Type &s) {
        ss << s;
        if (ss.str().size() >= bufSize) {
            flush();
        }
        return *this;
    }

    PrintStream &operator<<(std::ostream &(*f)(std::ostream &)) {
        ss << f;
        flush();
        return *this;
    }

    void flush() {
#if !defined(TINYMESH_PYTHON_MODULE)
        std::cout << ss.str();
        std::flush(std::cout);
#else
        py::print(ss.str(), py::arg("end") = "", py::arg("flush") = true);
#endif
        ss.str("");
    }

private:
    std::ostringstream ss;
    static const int bufSize = 128;
};

// -----------------------------------------------------------------------------
// Message handlers
// -----------------------------------------------------------------------------

template <typename... Args>
std::string STR_FMT(std::format_string<Args...> fmt, Args &&...args) {
    return std::format(fmt, std::forward<Args>(args)...);
}

template <typename... Args>
void Print(std::format_string<Args...> fmt, Args &&...args) {
    PrintStream() << STR_FMT(fmt, std::forward<Args>(args)...);
}

template <typename... Args>
void Info(std::format_string<Args...> fmt, Args &&...args) {
    PrintStream() << FG_BLU << "[INFO] " << FG_NC << STR_FMT(fmt, std::forward<Args>(args)...) << std::endl;
}

template <typename... Args>
void Debug(std::format_string<Args...> fmt, Args &&...args) {
    PrintStream() << FG_GRN << "[DEBUG] " << FG_NC << STR_FMT(fmt, std::forward<Args>(args)...) << std::endl;
}

template <typename... Args>
void Warn(std::format_string<Args...> fmt, Args &&...args) {
    PrintStream() << FG_YEL << "[WARNING] " << FG_NC << STR_FMT(fmt, std::forward<Args>(args)...) << std::endl;
}

template <typename... Args>
[[noreturn]] void Error(std::format_string<Args...> fmt, Args &&...args) {
#if defined(TINYMESH_PYTHON_MODULE)
    throw std::runtime_error(STR_FMT(fmt, std::forward<Args>(args)...));
#else
    PrintStream() << FG_RED << "[ERROR] " << FG_NC << STR_FMT(fmt, std::forward<Args>(args)...) << std::endl;
    std::terminate();
#endif
}

// -----------------------------------------------------------------------------
// Assertion with message
// -----------------------------------------------------------------------------

inline std::string ASRT_LABEL(const char *predicate, std::source_location location = std::source_location::current()) {
    return std::format("Assertion \"{}\" failed in {}:{} ({})", predicate, location.file_name(), location.line(),
                       location.function_name());
}

#if defined(TINYMESH_PYTHON_MODULE)
#define Assertion(PREDICATE, ...)                             \
    do {                                                      \
        if (!(PREDICATE)) {                                   \
            throw std::runtime_error(ASRT_LABEL(#PREDICATE)); \
        }                                                     \
    } while (false)
#elif !defined(TINYMESH_DISABLE_ASSERT)
#define Assertion(PREDICATE, ...)                             \
    do {                                                      \
        if (!(PREDICATE)) {                                   \
            std::cerr << ASRT_LABEL(#PREDICATE) << std::endl; \
            std::abort();                                     \
        }                                                     \
    } while (false)
#else  // TINYMESH_DISABLE_ASSERT
#define Assertion(PREDICATE, ...) \
    do {                          \
    } while (false)
#endif  // TINYMESH_DISABLE_ASSERT

#endif  // TINYMESH_DEBUG_H
