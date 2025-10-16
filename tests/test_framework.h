#ifndef TEST_FRAMEWORK_H
#define TEST_FRAMEWORK_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// Test statistics
static int tests_run = 0;
static int tests_passed = 0;
static int tests_failed = 0;

// Color codes for terminal output
#define COLOR_GREEN "\033[0;32m"
#define COLOR_RED "\033[0;31m"
#define COLOR_YELLOW "\033[0;33m"
#define COLOR_BLUE "\033[0;34m"
#define COLOR_RESET "\033[0m"

// Test assertion macros
#define TEST_ASSERT(condition, message) do { \
    tests_run++; \
    if (!(condition)) { \
        printf(COLOR_RED "✗ FAIL: %s\n" COLOR_RESET, message); \
        printf("  File: %s, Line: %d\n", __FILE__, __LINE__); \
        tests_failed++; \
        return 0; \
    } else { \
        tests_passed++; \
    } \
} while(0)

#define TEST_ASSERT_EQUAL(a, b, message) do { \
    tests_run++; \
    if ((a) != (b)) { \
        printf(COLOR_RED "✗ FAIL: %s\n" COLOR_RESET, message); \
        printf("  Expected: %d, Got: %d\n", (int)(b), (int)(a)); \
        printf("  File: %s, Line: %d\n", __FILE__, __LINE__); \
        tests_failed++; \
        return 0; \
    } else { \
        tests_passed++; \
    } \
} while(0)

#define TEST_ASSERT_DOUBLE_EQUAL(a, b, tolerance, message) do { \
    tests_run++; \
    if (fabs((a) - (b)) > (tolerance)) { \
        printf(COLOR_RED "✗ FAIL: %s\n" COLOR_RESET, message); \
        printf("  Expected: %.10e, Got: %.10e, Diff: %.10e\n", \
               (double)(b), (double)(a), fabs((a) - (b))); \
        printf("  File: %s, Line: %d\n", __FILE__, __LINE__); \
        tests_failed++; \
        return 0; \
    } else { \
        tests_passed++; \
    } \
} while(0)

#define TEST_ASSERT_NULL(ptr, message) do { \
    tests_run++; \
    if ((ptr) != NULL) { \
        printf(COLOR_RED "✗ FAIL: %s\n" COLOR_RESET, message); \
        printf("  File: %s, Line: %d\n", __FILE__, __LINE__); \
        tests_failed++; \
        return 0; \
    } else { \
        tests_passed++; \
    } \
} while(0)

#define TEST_ASSERT_NOT_NULL(ptr, message) do { \
    tests_run++; \
    if ((ptr) == NULL) { \
        printf(COLOR_RED "✗ FAIL: %s\n" COLOR_RESET, message); \
        printf("  File: %s, Line: %d\n", __FILE__, __LINE__); \
        tests_failed++; \
        return 0; \
    } else { \
        tests_passed++; \
    } \
} while(0)

// Test runner macros
#define RUN_TEST(test_func) do { \
    printf(COLOR_BLUE "Running: %s\n" COLOR_RESET, #test_func); \
    if (test_func()) { \
        printf(COLOR_GREEN "✓ PASS: %s\n\n" COLOR_RESET, #test_func); \
    } else { \
        printf(COLOR_RED "✗ FAIL: %s\n\n" COLOR_RESET, #test_func); \
    } \
} while(0)

#define TEST_SUITE_START(name) do { \
    printf("\n" COLOR_YELLOW "════════════════════════════════════════\n"); \
    printf("  TEST SUITE: %s\n", name); \
    printf("════════════════════════════════════════\n" COLOR_RESET); \
    tests_run = 0; \
    tests_passed = 0; \
    tests_failed = 0; \
} while(0)

#define TEST_SUITE_END() do { \
    printf(COLOR_YELLOW "════════════════════════════════════════\n"); \
    printf("  RESULTS\n"); \
    printf("════════════════════════════════════════\n" COLOR_RESET); \
    printf("Tests run:    %d\n", tests_run); \
    printf(COLOR_GREEN "Tests passed: %d\n" COLOR_RESET, tests_passed); \
    if (tests_failed > 0) { \
        printf(COLOR_RED "Tests failed: %d\n" COLOR_RESET, tests_failed); \
    } else { \
        printf("Tests failed: 0\n"); \
    } \
    double pass_rate = tests_run > 0 ? (100.0 * tests_passed / tests_run) : 0.0; \
    printf("Pass rate:    %.1f%%\n\n", pass_rate); \
    return (tests_failed == 0) ? 0 : 1; \
} while(0)

#endif // TEST_FRAMEWORK_H
