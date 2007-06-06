

#ifndef ATTRIB_GTC_EXPORT_H
#define ATTRIB_GTC_EXPORT_H

#if defined( _WIN32)
#if defined(gtcAttrib_EXPORTS)
#define DECL_GTCATTRIB __declspec(dllexport)
#else
#define DECL_GTCATTRIB __declspec(dllimport)
#endif
#else
#define DECL_GTCATTRIB
#endif


#endif
