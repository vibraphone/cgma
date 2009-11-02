#######################################################################################
# Macro to set the following variables:
# ACIS_VERSION
# ACIS_LIB_DIR
# ACIS_DEBUG_LIB_DIR
# ACIS_PLATFORM
# This macro expects ACIS_DIR to be set.
# If ACIS_SYSTEM is set, ACIS_LIB_DIR will
# be constructed using that value rather than
# via a trial-and-error search.
# If ACIS_VERSION is set, the corresponding test
# will be skipped.
#######################################################################################
AC_DEFUN([FATHOM_ACIS_ENV], [

AC_CHECK_FILE([$ACIS_DIR/include/acis.hxx], [], [AC_MSG_ERROR("Invalid ACIS path")])
 
if test "x" = "x$ACIS_SYSTEM"; then
  dir_list="$ACIS_DIR/bin/*"
else
  dir_list="$ACIS_DIR/bin/$ACIS_SYSTEM"
fi

AC_MSG_CHECKING([for ACIS library directory])
snl_acis_link_success=no
for dir in $dir_list; do
  if test "x$snl_acis_link_success" = "xyes"; then
    break
  fi
    # first try non-debug libraries
  case "$dir" in
    *_debug)
      ;;
    *)
      FATHOM_ACIS_LIBDIR( [$dir] )
      ;;
  esac
    # next try debug libraries
  case "$dir" in
    *_debug)
      FATHOM_ACIS_LIBDIR( [$dir] )
      ;;
  esac
done

if test "xyes" != "x$snl_acis_link_success"; then
  AC_MSG_ERROR([failed.])
else 
  ACIS_LIB_DIR="$dir"
  AC_MSG_RESULT([$ACIS_LIB_DIR])
  
  AC_MSG_CHECKING([ACIS version])
  AC_MSG_RESULT([$ACIS_VERSION])
  if test "x0" = "x$ACIS_VERSION"; then
    AC_MSG_ERROR([Failed to detect ACIS version.  Try --with-acis-version."])
  else 
    if test "0$ACIS_VERSION" -lt "0600"; then
      AC_MSG_ERROR([Invalid ACIS version.])
    fi
  fi
fi

# If either ACIS_LIB_DIR or ACIS_DEBUG_LIB_DIR is not defined,
# make both the same
if test "x" = "x$ACIS_LIB_DIR"; then
  ACIS_LIB_DIR="$ACIS_DEBUG_LIB_DIR"
elif test "x" = "x$ACIS_DEBUG_LIB_DIR"; then
  ACIS_DEBUG_LIB_DIR="$ACIS_LIB_DIR"
fi

# Determine the ACIS platform name from the lib dir
AC_MSG_CHECKING([ACIS platform name])
case "$ACIS_LIB_DIR" in
  $ACIS_DIR/bin/aix*)  
    ACIS_PLATFORM=aix
    ;;
  $ACIS_DIR/bin/hp700*)
    ACIS_PLATFORM=hp700
    ;;
  $ACIS_DIR/bin/linux*)
    ACIS_PLATFORM=linux
    ;;
  $ACIS_DIR/bin/mac*)
    ACIS_PLATFORM=mac
    ;;
  $ACIS_DIR/bin/osf1*)
    ACIS_PLATFORM=osf1
    ;;
  $ACIS_DIR/bin/sgi*)
    ACIS_PLATFORM=sgi
    ;;
  $ACIS_DIR/bin/solaris*)
    ACIS_PLATFORM=solaris
    ;;
  *)
    AC_MSG_ERROR([Cannot determine ACIS platform name.])
    ;;
esac
AC_MSG_RESULT([$ACIS_PLATFORM])
])



#######################################################################################
# Macro to check for ACIS translators.
#######################################################################################
AC_DEFUN([FATHOM_ACIS_TRANSLATOR], [
AC_REQUIRE([FATHOM_ACIS_ENV])
case "$ACIS_VERSION" in
  11??)
    ACIS_XLATE_LIBS='-lxacis2k -lxcore2k -lxutil'
    ;;
  1[23]??)
    ACIS_XLATE_LIBS='-lSPAXAssemblyRep -lSPAXInterop -lSPAXBase -lxacis2k -lxcore2k -lxutil'
    ;;
  14??)
    ACIS_XLATE_LIBS='-lSPAXAssemblyRep -lSPAXInterop -lSPAXAcisBase -lSPAXDefaultGeometryRep -lSPAXGeometryRepresentation -lSPAXPropertiesBRepImporter -lSPAXPropertiesBase -lSPAXBase -lxacis2k -lxcore2k'
    ;;
  *)
    ACIS_XLATE_LIBS='-lSPAXEBOMBase -lSPAXEBOMNewAssemblyExporter -lSPAXEBOMNewAssemblyImporter -lSPAXXMLTk -lpthread -lSPAXPropertiesNewAssemblyImporter -lSPAXAssemblyRep -lSPAXInterop -lSPAXAcisBase -lSPAXDefaultGeometryRep -lSPAXGeometryRepresentation -lSPAXPropertiesBRepImporter -lSPAXPropertiesBase -lSPAXBase -lxacis2k -lxcore2k'
    ;;
esac
old_LIBS="$LIBS"
LIBS="$LIBS -L$ACIS_LIB_DIR $ACIS_XLATE_LIBS $ACIS_BASE_LIBS"
AC_CHECK_LIB([xstep],[main],[ACIS_STEP_TRANSLATOR=-DACIS_STEP_TRANSLATOR])
AC_CHECK_LIB([xiges],[main],[ACIS_IGES_TRANSLATOR=-DACIS_IGES_TRANSLATOR])
LIBS="$old_LIBS"
])


#######################################################################################
# *******************************************************************************
# **************************** INTERNAL STUFF ***********************************
# *******************************************************************************
#######################################################################################

#######################################################################################
# Macro to get ACIS_VERSION
# Expected arguments: ACIS library path
# Sets variables:
#  snl_acis_link_success=yes/no
#  ACIS_VERSION=version/0
#  If ACIS_VERSION is alreay set, it will NOT be modified by
#  this function.
#######################################################################################
AC_DEFUN([FATHOM_ACIS_LIBDIR], [
AC_REQUIRE([AC_PROG_LIBTOOL])
AC_LANG_PUSH([C++])
snl_acis_libdir=$1
snl_acis_link_success=no
old_LDFLAGS="$LDFLAGS"
old_LIBS="$LIBS"
LDFLAGS="-L$1"
LIBS="-lSpaBase"
AC_LINK_IFELSE(
 [AC_LANG_PROGRAM([[
#include <stdio.h>
int get_major_version();
int get_minor_version();
int get_point_version();
]],[[
printf("%d\n", 
100*get_major_version() +
 10*get_minor_version() +
    get_point_version());  
]])],
[snl_acis_link_success=yes
 if test "x" = "x$ACIS_VERSION"; then
   old_LD_LIBRARY_PATH="$LD_LIBRARY_PATH"
   old_SHLIB_PATH="$SH_LIBPATH"
   old_PATH="$PATH"  # windows!
   LD_LIBRARY_PATH="${snl_acis_libdir}:$LD_LIBRARY_PATH"
   SHLIB_PATH="${snl_acis_libdir}:$SHLIB_PATH"
   PATH="${snl_acis_libdir}:$PATH"
   export LD_LIBRARY_PATH
   export SHLIB_PATH
   export PATH
   ACIS_VERSION=`./conftest` || ACIS_VERSION=0
   LD_LIBRARY_PATH="$old_LD_LIBRARY_PATH"
   SHLIB_PATH="$old_SHLIB_PATH"
   PATH="$old_PATH"
   export LD_LIBRARY_PATH
   export SHLIB_PATH
   export PATH
 fi
], [snl_acis_link_success=no])
AC_LANG_POP
]) # FATHOM_ACIS_LIBDIR
