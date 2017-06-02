
 #
 # Open-MPI ("stable" install) environment settings
 #
 # TJN - Sun Oct 21 15:15:04 EDT 2012
 #

basedir="$HOME/projects/ompi/openmpi-stable/local"
ompi_lib_path="${basedir}/lib"
ompi_bin_path="${basedir}/bin"

#### END CONFIG ####

if ! `echo $PATH | grep -q ${ompi_bin_path}` ; then
    export PATH="${ompi_bin_path}:$PATH"
fi

if [ "x$LD_LIBRARY_PATH" == "x" ] ; then
    export LD_LIBRARY_PATH="${ompi_lib_path}"
else
    if ! `echo $LD_LIBRARY_PATH | grep -q ${ompi_lib_path}` ; then
        export LD_LIBRARY_PATH="${ompi_lib_path}:$LD_LIBRARY_PATH"
    fi
fi

