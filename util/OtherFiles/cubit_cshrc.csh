#!/bin/csh
#This is meant to be included in your .cshrc file 
if( ! $?MACHINE_TYPE ) then
	echo "Trying to set your MACHINE_TYPE"
        set UNAME=`which uname`
        if( $status || ! -x $UNAME ) then
                echo "Cannot determine platform"
                echo "You must manually set MACHINE_TYPE"
                exit
        endif

        set platform=`uname`
        if( "$platform" == "SunOS" ) then
                setenv MACHINE_TYPE ss
        else if ( "$platform" == "IRIX" ) then
                setenv MACHINE_TYPE sg
        else if ( "$platform" == "HP-UX" ) then
                setenv MACHINE_TYPE hp
        else if ( "$platform" == "AIX" ) then
                setenv MACHINE_TYPE ibm
        else if ( "$platform" == "Linux" ) then
                setenv MACHINE_TYPE lin
        else if ( "$platform" == "OSF1" ) then
                setenv MACHINE_TYPE da
        else 
                echo "Unknown platform $platform"
                echo "You must manually set MACHINE_TYPE"
                exit
        endif
endif

if( ! $?MACHINE_TYPE ) then
	echo "Your MACHINE_TYPE is not set"
	echo "You can't compile Cubit"
	exit
endif
if( ! $?LD_LIBRARY_PATH ) setenv LD_LIBRARY_PATH ""
if( ! $?MANPATH ) setenv MANPATH ""

#DEC/Compaq Alpha
if ( $MACHINE_TYPE == "da" ) then
	setenv ACIS_VERSION_PATH 7.0.5
	setenv LD_LIBRARY_PATH /usr/local/eng_sci/cubit/acis/acis$ACIS_VERSION_PATH/lib/osf1_64_so:/usr/local/eng_sci/cubit/hoops/hoops-6.20/lib:$LD_LIBRARY_PATH
	setenv CUBIT_STEP_PATH /usr/local/eng_sci/cubit/acis/acis$ACIS_VERSION_PATH/step/tools/osf
	setenv CUBIT_IGES_PATH /usr/local/eng_sci/cubit/acis/acis$ACIS_VERSION_PATH/iges/igac/iges_04.dbt
	setenv ARCH osf1
	setenv arch $ARCH

#Linux machines (lot's of differences in this part--customize to fit)
else if ( $MACHINE_TYPE == "lin" ) then
	setenv ACIS_VERSION_PATH 7.0.5
	setenv LD_LIBRARY_PATH /usr/local/eng_sci/cubit/hoops/hoops-6.20/lib:/usr/local/eng_sci/cubit/acis/acis$ACIS_VERSION_PATH/lib/linux_so:$LD_LIBRARY_PATH 
	setenv CUBIT_STEP_PATH /usr/local/eng_sci/cubit/acis/acis$ACIS_VERSION_PATH/step/tools/linux
	setenv CUBIT_IGES_PATH /usr/local/eng_sci/cubit/acis/acis$ACIS_VERSION_PATH/iges/igac/iges_04.dbt
	setenv ARCH linux
	setenv arch Linux

#Solaris 8 with workshop 6 update 2
else if ( $MACHINE_TYPE == "ss" ) then
	setenv ACIS_VERSION_PATH 7.0.5
	set path = (/net/sasn232-atm/opt/SUNWspro6.0_update2/SUNWspro/bin $path )
	setenv LD_LIBRARY_PATH /usr/local/eng_sci/cubit/acis/acis$ACIS_VERSION_PATH/lib/solaris_so:/net/sasn232-atm/opt/SUNWspro6.0_update2/SUNWspro/WS6U2/lib:$LD_LIBRARY_PATH
	setenv MANPATH /net/sasn232-atm/opt/SUNWspro6.0_update2/SUNWspro/man:$MANPATH
	setenv CUBIT_STEP_PATH /usr/local/eng_sci/cubit/acis/acis$ACIS_VERSION_PATH/step/tools/solaris
	setenv CUBIT_IGES_PATH /usr/local/eng_sci/cubit/acis/acis$ACIS_VERSION_PATH/iges/igac/iges_04.dbt
	setenv ARCH solaris
	setenv arch $ARCH

#SGI IRIX 6.5
else if ( $MACHINE_TYPE == "sg" ) then
#  For 32-bit, use the following line:
#	setenv DEFAULT_BINARY_VERSION ""
#  For 64-bit, use the following line:
	setenv DEFAULT_BINARY_VERSION 64_
	setenv ACIS_VERSION_PATH 7.0.5
	setenv LD_LIBRARY_PATH /usr/local/eng_sci/cubit/acis/acis$ACIS_VERSION_PATH/lib/sgi_$DEFAULT_BINARY_VERSION\so:$LD_LIBRARY_PATH
	set path = (/usr/bsd $path)
	setenv MANPATH $MANPATH\:/opt/modules/2.2.2.4/man:/usr/ToolTalk:/usr/share/man:/usr/local/man
	setenv CUBIT_STEP_PATH /usr/local/eng_sci/cubit/acis/acis$ACIS_VERSION_PATH/step/tools/sgi
	setenv CUBIT_IGES_PATH /usr/local/eng_sci/cubit/acis/acis$ACIS_VERSION_PATH/iges/igac/iges_04.dbt
	setenv ARCH sgi
	setenv arch $ARCH

#HP-UX 10.20
else if ( $MACHINE_TYPE == "hp" ) then
	setenv ACIS_VERSION_PATH 7.0.5
	setenv MANPATH /usr/share/man/%L:/usr/share/man:/usr/contrib/man/%L:/usr/contrib/man:/usr/local/man/%L:/usr/local/man:/opt/pd/share/man/%L:/opt/pd/share/man:/opt/ignite/share/man/%L:/opt/ignite/share/man:/opt/resmon/share/man:/opt/hparray/share/man/%L:/opt/hparray/share/man:/opt/graphics/starbase/share/man:/opt/audio/share/man:/opt/blinklink/share/man:/opt/ansic/share/man/%L:/opt/ansic/share/man:/opt/langtools/share/man/%L:/opt/langtools/share/man:/opt/video/share/man:/opt/videoout/share/man:/opt/fortran/share/man/%L:/opt/fortran/share/man:/opt/imake/man:/opt/graphics/PEX5/share/man:/opt/aCC/share/man/%L:/opt/CC/share/man:/opt/softbench/share/man/%L:/opt/softbench/share/man:/opt/aCC/share/man:$MANPATH
	set path = ( /usr/local/bin /opt/ansic/bin $path)
	if !($?SHLIB_PATH) setenv SHLIB_PATH ""
	setenv SHLIB_PATH $SHLIB_PATH\:/usr/local/eng_sci/cubit/hoops/hoops-6.1/lib:/usr/local/eng_sci/cubit/acis/acis$ACIS_VERSION_PATH/lib/hp700_so
	setenv LD_LIBRARY_PATH $SHLIB_PATH:$LD_LIBRARY_PATH
	setenv CUBIT_STEP_PATH /usr/local/eng_sci/cubit/acis/acis$ACIS_VERSION_PATH/step/tools/hp700
	setenv CUBIT_IGES_PATH /usr/local/eng_sci/cubit/acis/acis$ACIS_VERSION_PATH/iges/igac/iges_04.dbt
	setenv ARCH hp700
	setenv arch $ARCH

endif
