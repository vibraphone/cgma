
iGeom_SIDL_includes = \
  -I/mnt/disk2b/jhu/merge-cubit14.0/itaps/SIDL \
  -I/mnt/disk2b/jhu/merge-cubit14.0/itaps/SIDL/gserver \
  -I/mnt/disk2b/jhu/merge-cubit14.0/itaps/SIDL/Cclient \
  -I/mnt/disk2b/jhu/merge-cubit14.0/itaps/SIDL/Fclient \
  -I/mnt/disk2b/jhu/merge-cubit14.0/itaps/SIDL/.. 
  
iGeom_SIDL_INCLUDES = \
  $(iGeom_SIDL_includes) \
  -I/include
  
iGeom_SIDL_LDFLAGS = \
  -L/mnt/disk2b/jhu/merge-cubit14.0/itaps/SIDL/../.libs \
  -L/mnt/disk2b/jhu/merge-cubit14.0/itaps/SIDL/../../../.libs \
  -L/mnt/disk2b/jhu/merge-cubit14.0/itaps/SIDL/gserver/.libs \
  -L/mnt/disk2b/jhu/merge-cubit14.0/itaps/SIDL/Cclient/.libs \
  -L/mnt/disk2b/jhu/merge-cubit14.0/itaps/SIDL/Fclient/.libs
  
iGeom_SIDL_LIBS = -L/lib $(iGeom_SIDL_LDFLAGS) $(iGeom_SIDL_libs) \
                  -liGeomserver -lsidl -liGeom -lcgm  -lTKSTL -lTKSTEP -lTKSTEP209 -lTKSTEPAttr -lTKSTEPBase -lTKXSBase -lTKIGES -lTKXSBase -lTKBinL -lTKLCAF -lTKCDF -lTKCAF -lTKHLR -lTKOffset -lTKShHealing -lTKFillet -lTKFeat -lTKBool -lTKBO -lTKPrim -lTKMesh -lTKTopAlgo -lTKGeomAlgo -lTKBRep -lTKGeomBase -lTKG3d -lTKG2d -lTKMath -lTKernel -lmpichcxx -lmpich -lopa -lmpl -lrt -lpthread -lstdc++ -lm -lgcc_s -lc -lgcc_s

iGeom_SIDL_FILES = /mnt/disk2b/jhu/merge-cubit14.0/itaps/SIDL/libiGeomserver.la
