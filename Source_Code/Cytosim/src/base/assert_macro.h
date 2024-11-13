// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University
/**
 @file
 @brief Cytosim's assertions
 */

#ifndef ASSERT_MACRO_H
#define ASSERT_MACRO_H

#include <cstdio>
#include <cstdlib>
#include <cstring>

/**
 Assertions are used extensively to check the validity of data and function arguments.
 Defining NDEBUG disables:
    - the standard assert() macro and,
    - the custom assert_true() and assert_false() macros defined below.
 This makes the executable faster.
 Note that the value defined for NDEBUG is ignored: 0 or 1 have the same effects!
 */
//#define NDEBUG 1


/// strips the full pathname for a file name
#define SFILE strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__


/// print the current function name:
//#define SFUNC __func__
#define SFUNC __PRETTY_FUNCTION__

//------------------------------- ASSERTIONS -----------------------------------
/** 
  assert_true(X) stops the program if condition X is false.
  assert_false(X) stops if X is true, and prints the value of X.
*/

#ifdef NDEBUG

  #define assert_true(ignore)  ((void) 0)
  #define assert_false(ignore) ((void) 0)
  #define assert_small(ignore) ((void) 0)

#else

  #define assert_true(expression)\
        if (!(expression)) {\
            fprintf(stderr, "\n*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *\n");\
            fprintf(stderr, "Cytosim failed assert_true(%s)\n", #expression);\
            fprintf(stderr, "`%s' in %s:%d", SFUNC, SFILE, __LINE__);\
            fprintf(stderr, "\n*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *\n");\
            abort();\
        }

  #define assert_false(expression)\
        { int e = expression;\
        if (e) {\
            fprintf(stderr, "\n*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *\n");\
            fprintf(stderr, "Cytosim failed assert_false(%s) with value %i\n", #expression, e);\
            fprintf(stderr, "`%s' in %s:%d", SFUNC, SFILE, __LINE__);\
            fprintf(stderr, "\n*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *\n");\
            abort();\
        } }

  #define assert_small(expression)\
        { real e = expression;\
        if ( std::fabs(e) > 0.01 ) {\
            fprintf(stderr, "//  Cytosim failed assert_small(%s) with value %e\n", #expression, e);\
            fprintf(stderr, "\\\\  within `%s' in %s:%d\n", SFUNC, SFILE, __LINE__);\
        } }

#endif


//-------------------------- ERROR HANDLING MACROS -----------------------------

/// macro to abort after printing an error message
#define ABORT_NOW(message)\
    {\
    fprintf(stderr, "Cytosim ERROR `%s'\n", message);\
    fprintf(stderr, "@ `%s' in %s:%d\n", SFUNC, SFILE, __LINE__);\
    exit(EXIT_FAILURE);\
    }

#endif
