/*
 * Small macro for printing debug messages.
 * The level of detail is controlled by DLM_LOG_LEVEL, which can be specified at compile-time using:
 *   cmake -DDLM_LOG_LEVEL=40
 * By default, DLM_LOG_LEVEL=10, namely only errors and fatal errors are shown.
*/

#ifndef DLM_LOGGER_H
#define DLM_LOGGER_H

// Non solvable error that must break the execution of the program (throws std::runtime_error)
#define FATAL 0
// A problem that can be taken care of, in some way
#define ERROR 10
// Something that might lead to wrong or inconsistent output, deprecated code, experimental or unstable features. The
// user must be carefull and ensure that the output makes sense
#define WARN 20
// Status and progress of some calculation or process. Should not be too verbose
#define INFO 30
// Detailed output about technical stuff. Can be (very) verbose
#define DEBUG 40

#ifndef DLM_LOG_LEVEL
#define DLM_LOG_LEVEL ERROR  // By default, only show problems
#endif

#include <iostream>

#define DLM_LOG_LEVEL_STRING(level) \
    ((level) == FATAL   ? "FATAL"   \
     : (level) == ERROR ? "ERROR"   \
     : (level) == WARN  ? "WARNING" \
     : (level) == INFO  ? "INFO"    \
     : (level) == DEBUG ? "DEBUG"   \
                        : "UNKNOWN")

#define LOG(level, msg)                                                                                        \
    do {                                                                                                       \
        if (level <= DLM_LOG_LEVEL) {                                                                          \
            std::cerr << "[" << DLM_LOG_LEVEL_STRING(level) << "]"                                             \
                      << "[" << __FILE__ << " -- " << __func__ << ":" << __LINE__ << "] " << msg << std::endl; \
            if (level == FATAL) {                                                                              \
                throw std::runtime_error(std::string(msg));                                                    \
            }                                                                                                  \
        }                                                                                                      \
    } while (0)

#endif
