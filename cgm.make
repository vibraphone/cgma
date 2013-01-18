#==================================================================
#
# To build CGM applications, do the following:
#
# 0. include this file in your makefile (using 'include cgm.make')
# 1. Insert '${CGM_INCLUDES}' (w/o quotes) in your compile command, e.g.
#      ${CXX} ${CGM_INCLUDES} -c mysource.cpp
# 2. Insert '${CGM_LIBS_LINK}' (w/o quotes) on your link line
#
# That's it! No need to look at the code below, unless you're curious.
#
#==================================================================

# The following are initialized here for a non-installed CGM.
# These values will be overridden when this file is installed.
CGM_INT_LDFLAGS = -L/mnt/disk2b/jhu/merge-cubit14.0/.libs
CGM_INT_LTFLAGS = -R/mnt/disk2b/jhu/merge-cubit14.0/.libs
CGM_INT_INCLUDE = -I/mnt/disk2b/jhu/merge-cubit14.0 \
                  -I/mnt/disk2b/jhu/merge-cubit14.0/util \
                  -I/mnt/disk2b/jhu/merge-cubit14.0/util \
                  -I/mnt/disk2b/jhu/merge-cubit14.0/init \
                  -I/mnt/disk2b/jhu/merge-cubit14.0/init \
                  -I/mnt/disk2b/jhu/merge-cubit14.0/geom \
                  -I/mnt/disk2b/jhu/merge-cubit14.0/geom \
                  -I/mnt/disk2b/jhu/merge-cubit14.0/amendment \
                  -I/mnt/disk2b/jhu/merge-cubit14.0/amendment \
                  -I/mnt/disk2b/jhu/merge-cubit14.0/geom/ACIS \
                  -I/mnt/disk2b/jhu/merge-cubit14.0/geom/virtual \
                  -I/mnt/disk2b/jhu/merge-cubit14.0/geom/facet \
                  -I/mnt/disk2b/jhu/merge-cubit14.0/geom/facetbool \
                  -I/mnt/disk2b/jhu/merge-cubit14.0/geom/Cholla 

# Pre-processor flags
CGM_DEFINES =   -DTEMPLATE_DEFS_INCLUDED  -DHAVE_OCC -DHAVE_OCC_IGES -DHAVE_OCC_STEP -DHAVE_OCC_STL
CGM_INCLUDES =  -I/mnt/disk2b/jhu/build653/inc -D_OCC64 -DHAVE_IOSTREAM -DHAVE_IOMANIP -DHAVE_FSTREAM -DHAVE_LIMITS_H $(CGM_INT_INCLUDE)
CGM_CPPFLAGS = $(CGM_DEFINES) $(CGM_INCLUDES)

# Link flags
CGM_LIBS = -lcgm  -lTKSTL -lTKSTEP -lTKSTEP209 -lTKSTEPAttr -lTKSTEPBase -lTKXSBase -lTKIGES -lTKXSBase -lTKBinL -lTKLCAF -lTKCDF -lTKCAF -lTKHLR -lTKOffset -lTKShHealing -lTKFillet -lTKFeat -lTKBool -lTKBO -lTKPrim -lTKMesh -lTKTopAlgo -lTKGeomAlgo -lTKBRep -lTKGeomBase -lTKG3d -lTKG2d -lTKMath -lTKernel
CGM_LDFLAGS = $(CGM_INT_LDFLAGS)  -L/mnt/disk2b/jhu/build653/lib
CGM_LTFLAGS = $(CGM_INT_LTFLAGS) 
CGM_LIBS_LINK = $(CGM_LDFLAGS) $(CGM_LIBS)

# Build-generated values appended after this line

