//
// MATLAB Compiler: 5.2 (R2014b)
// Date: Tue Jun  9 11:50:18 2015
// Arguments: "-B" "macro_default" "-g" "-G" "-W" "cpplib:libStraticounter"
// "-T" "link:lib" "-d" "MCR_Library/" "-I" "Subroutines/" "-I"
// "Subroutines/algorithm/" "-I" "Subroutines/batchresults/" "-I"
// "Subroutines/checkinput/" "-I" "Subroutines/layertemplates/" "-I"
// "Subroutines/matchmaker/" "-I" "Subroutines/preprocessing/" "-I"
// "Subroutines/showresults/" "-I" "Subroutines/syntheticdata/" "-I"
// "Subroutines/utilities/" "straticounter_scibox.m"
//

#ifndef __libStraticounter_h
#define __libStraticounter_h 1

#if defined(__cplusplus) && !defined(mclmcrrt_h) && defined(__linux__)
#  pragma implementation "mclmcrrt.h"
#endif
#include "mclmcrrt.h"
#include "mclcppclass.h"
#ifdef __cplusplus
extern "C" {
#endif

#if defined(__SUNPRO_CC)
/* Solaris shared libraries use __global, rather than mapfiles
 * to define the API exported from a shared library. __global is
 * only necessary when building the library -- files including
 * this header file to use the library do not need the __global
 * declaration; hence the EXPORTING_<library> logic.
 */

#ifdef EXPORTING_libStraticounter
#define PUBLIC_libStraticounter_C_API __global
#else
#define PUBLIC_libStraticounter_C_API /* No import statement needed. */
#endif

#define LIB_libStraticounter_C_API PUBLIC_libStraticounter_C_API

#elif defined(_HPUX_SOURCE)

#ifdef EXPORTING_libStraticounter
#define PUBLIC_libStraticounter_C_API __declspec(dllexport)
#else
#define PUBLIC_libStraticounter_C_API __declspec(dllimport)
#endif

#define LIB_libStraticounter_C_API PUBLIC_libStraticounter_C_API


#else

#define LIB_libStraticounter_C_API

#endif

/* This symbol is defined in shared libraries. Define it here
 * (to nothing) in case this isn't a shared library.
 */
#ifndef LIB_libStraticounter_C_API
#define LIB_libStraticounter_C_API /* No special import/export declaration */
#endif

extern LIB_libStraticounter_C_API
bool MW_CALL_CONV libStraticounterInitializeWithHandlers(
       mclOutputHandlerFcn error_handler,
       mclOutputHandlerFcn print_handler);

extern LIB_libStraticounter_C_API
bool MW_CALL_CONV libStraticounterInitialize(void);

extern LIB_libStraticounter_C_API
void MW_CALL_CONV libStraticounterTerminate(void);



extern LIB_libStraticounter_C_API
void MW_CALL_CONV libStraticounterPrintStackTrace(void);

extern LIB_libStraticounter_C_API
bool MW_CALL_CONV mlxStraticounter_scibox(int nlhs, mxArray *plhs[], int nrhs, mxArray
                                          *prhs[]);


#ifdef __cplusplus
}
#endif

#ifdef __cplusplus

/* On Windows, use __declspec to control the exported API */
#if defined(_MSC_VER) || defined(__BORLANDC__)

#ifdef EXPORTING_libStraticounter
#define PUBLIC_libStraticounter_CPP_API __declspec(dllexport)
#else
#define PUBLIC_libStraticounter_CPP_API __declspec(dllimport)
#endif

#define LIB_libStraticounter_CPP_API PUBLIC_libStraticounter_CPP_API

#else

#if !defined(LIB_libStraticounter_CPP_API)
#if defined(LIB_libStraticounter_C_API)
#define LIB_libStraticounter_CPP_API LIB_libStraticounter_C_API
#else
#define LIB_libStraticounter_CPP_API /* empty! */
#endif
#endif

#endif

extern LIB_libStraticounter_CPP_API void MW_CALL_CONV straticounter_scibox(const mwArray& settings_path, const mwArray& output_path);

#endif
#endif
