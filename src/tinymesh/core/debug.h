#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_DEBUG_H
#define TINYMESH_DEBUG_H

// -----------------------------------------------------------------------------
// Message handlers
// -----------------------------------------------------------------------------

#define Info(...)                     \
    do {                              \
        std::cout << "[INFO] ";       \
        fprintf(stdout, __VA_ARGS__); \
        std::cout << std::endl;       \
    } while (false);

#define Warn(...)                  \
    do {                              \
        std::cerr << "[WARNING] ";    \
        fprintf(stderr, __VA_ARGS__); \
        std::cerr << std::endl;       \
    } while (false);

#define FatalError(...)               \
    do {                              \
        std::cerr << "[ERROR] ";      \
        fprintf(stderr, __VA_ARGS__); \
        std::cerr << std::endl;       \
        std::abort();                 \
    } while (false);

#define Debug(...)                    \
    do {                              \
        std::cerr << "[DEBUG] ";      \
        fprintf(stdout, __VA_ARGS__); \
        std::cerr << std::endl;       \
} while (false);


#endif  // TINYMESH_DEBUG_H

