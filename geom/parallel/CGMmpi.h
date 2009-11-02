#ifndef CGM_MPI_H
#define CGM_MPI_H
#include "CGMmpi_config.h"

#ifndef __cplusplus
#  include <mpi.h>
#elif !defined(CGM_MPI_CXX_CONFLICT)
#  ifndef MPICH_IGNORE_CXX_SEEK
#    define MPICH_IGNORE_CXX_SEEK
#  endif
#  include <mpi.h>
#else
#  include <stdio.h>
#  ifdef SEEK_SET
#    undef SEEK_SET
#    ifdef CGM_SEEK_SET
#      define CGM_RESTORE_SEEK_SET
#    endif
#  endif
#  ifdef SEEK_CUR
#    undef SEEK_CUR
#    ifdef CGM_SEEK_CUR
#      define CGM_RESTORE_SEEK_CUR
#    endif
#  endif
#  ifdef SEEK_END
#    undef SEEK_END
#    ifdef CGM_SEEK_END
#      define CGM_RESTORE_SEEK_END
#    endif
#  endif
#  include <mpi.h>
#  ifdef CGM_RESTORE_SEEK_SET
#    undef CGM_RESTORE_SEEK_SET
#    define SEEK_SET CGM_SEEK_SET
#  endif
#  ifdef CGM_RESTORE_SEEK_CUR
#    undef CGM_RESTORE_SEEK_CUR
#    define SEEK_CUR CGM_SEEK_CUR
#  endif
#  ifdef CGM_RESTORE_SEEK_END
#    undef CGM_RESTORE_SEEK_END
#    define SEEK_END CGM_SEEK_END
#  endif
#endif


#endif
