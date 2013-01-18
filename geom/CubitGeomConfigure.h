
#ifndef CUBIT_GEOM_CONFIGURE_H
#define CUBIT_GEOM_CONFIGURE_H

/* #undef CUBIT_GEOM_BUILD_SHARED_LIBS */

#if defined(CUBIT_GEOM_BUILD_SHARED_LIBS)
#if defined(cubit_geom_EXPORTS)
# if defined(WIN32)
#  define CUBIT_GEOM_EXPORT __declspec(dllexport)
# elif defined(__GNUC__) && __GNUC__ >= 4
#  define CUBIT_GEOM_EXPORT __attribute__ ((visibility("default")))
# endif
#else
# if defined(WIN32)
#  define CUBIT_GEOM_EXPORT __declspec(dllimport)
# endif
#endif
#endif

#ifndef CUBIT_GEOM_EXPORT
#define CUBIT_GEOM_EXPORT
#endif

#endif

