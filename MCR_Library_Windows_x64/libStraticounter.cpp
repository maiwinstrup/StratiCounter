//
// MATLAB Compiler: 5.2 (R2014b)
// Date: Fri Jun 19 13:44:11 2015
// Arguments: "-B" "macro_default" "-W" "cpplib:libStraticounter" "-T"
// "link:lib" "-d" "MCR_Library_Windows_x64/" "-I" "Subroutines/" "-I"
// "Subroutines/algorithm/" "-I" "Subroutines/batchresults/" "-I"
// "Subroutines/checkinput/" "-I" "Subroutines/layertemplates/" "-I"
// "Subroutines/matchmaker/" "-I" "Subroutines/preprocessing/" "-I"
// "Subroutines/showresults/" "-I" "Subroutines/syntheticdata/" "-I"
// "Subroutines/utilities/" "straticounter_scibox.m" 
//

#include <stdio.h>
#define EXPORTING_libStraticounter 1
#include "libStraticounter.h"

static HMCRINSTANCE _mcr_inst = NULL;


#if defined( _MSC_VER) || defined(__BORLANDC__) || defined(__WATCOMC__) || defined(__LCC__)
#ifdef __LCC__
#undef EXTERN_C
#endif
#include <windows.h>

static char path_to_dll[_MAX_PATH];

BOOL WINAPI DllMain(HINSTANCE hInstance, DWORD dwReason, void *pv)
{
    if (dwReason == DLL_PROCESS_ATTACH)
    {
        if (GetModuleFileName(hInstance, path_to_dll, _MAX_PATH) == 0)
            return FALSE;
    }
    else if (dwReason == DLL_PROCESS_DETACH)
    {
    }
    return TRUE;
}
#endif
#ifdef __cplusplus
extern "C" {
#endif

static int mclDefaultPrintHandler(const char *s)
{
  return mclWrite(1 /* stdout */, s, sizeof(char)*strlen(s));
}

#ifdef __cplusplus
} /* End extern "C" block */
#endif

#ifdef __cplusplus
extern "C" {
#endif

static int mclDefaultErrorHandler(const char *s)
{
  int written = 0;
  size_t len = 0;
  len = strlen(s);
  written = mclWrite(2 /* stderr */, s, sizeof(char)*len);
  if (len > 0 && s[ len-1 ] != '\n')
    written += mclWrite(2 /* stderr */, "\n", sizeof(char));
  return written;
}

#ifdef __cplusplus
} /* End extern "C" block */
#endif

/* This symbol is defined in shared libraries. Define it here
 * (to nothing) in case this isn't a shared library. 
 */
#ifndef LIB_libStraticounter_C_API
#define LIB_libStraticounter_C_API /* No special import/export declaration */
#endif

LIB_libStraticounter_C_API 
bool MW_CALL_CONV libStraticounterInitializeWithHandlers(
    mclOutputHandlerFcn error_handler,
    mclOutputHandlerFcn print_handler)
{
    int bResult = 0;
  if (_mcr_inst != NULL)
    return true;
  if (!mclmcrInitialize())
    return false;
  if (!GetModuleFileName(GetModuleHandle("libStraticounter"), path_to_dll, _MAX_PATH))
    return false;
    {
        mclCtfStream ctfStream = 
            mclGetEmbeddedCtfStream(path_to_dll);
        if (ctfStream) {
            bResult = mclInitializeComponentInstanceEmbedded(   &_mcr_inst,
                                                                error_handler, 
                                                                print_handler,
                                                                ctfStream);
            mclDestroyStream(ctfStream);
        } else {
            bResult = 0;
        }
    }  
    if (!bResult)
    return false;
  return true;
}

LIB_libStraticounter_C_API 
bool MW_CALL_CONV libStraticounterInitialize(void)
{
  return libStraticounterInitializeWithHandlers(mclDefaultErrorHandler, 
                                                mclDefaultPrintHandler);
}

LIB_libStraticounter_C_API 
void MW_CALL_CONV libStraticounterTerminate(void)
{
  if (_mcr_inst != NULL)
    mclTerminateInstance(&_mcr_inst);
}

LIB_libStraticounter_C_API 
void MW_CALL_CONV libStraticounterPrintStackTrace(void) 
{
  char** stackTrace;
  int stackDepth = mclGetStackTrace(&stackTrace);
  int i;
  for(i=0; i<stackDepth; i++)
  {
    mclWrite(2 /* stderr */, stackTrace[i], sizeof(char)*strlen(stackTrace[i]));
    mclWrite(2 /* stderr */, "\n", sizeof(char)*strlen("\n"));
  }
  mclFreeStackTrace(&stackTrace, stackDepth);
}


LIB_libStraticounter_C_API 
bool MW_CALL_CONV mlxStraticounter_scibox(int nlhs, mxArray *plhs[], int nrhs, mxArray 
                                          *prhs[])
{
  return mclFeval(_mcr_inst, "straticounter_scibox", nlhs, plhs, nrhs, prhs);
}

LIB_libStraticounter_CPP_API 
void MW_CALL_CONV straticounter_scibox(const mwArray& settings_path, const mwArray& 
                                       output_path)
{
  mclcppMlfFeval(_mcr_inst, "straticounter_scibox", 0, 0, 2, &settings_path, &output_path);
}

