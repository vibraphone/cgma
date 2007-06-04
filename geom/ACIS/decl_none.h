//*************************************************************
//  Added For ACIS 5.0
//
//*************************************************************

#ifndef DECL_NONE

#ifdef ACIS_DLL
# ifdef EXPORT_NONE
#  define DECL_NONE __declspec(dllexport)
# else
#  define DECL_NONE __declspec(dllimport)
# endif
#else
# define DECL_NONE
#endif

#ifndef EXPORT_NONE
# ifdef NT

#if CUBIT_ACIS_VERSION < 1100
#  pragma comment( lib, "kernel.lib") /* force link in VC++ */
#endif

# endif
#endif

#endif

