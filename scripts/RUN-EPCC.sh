#!/bin/bash
# Thomas Naughton <naughtont@ornl.gov>

SCRIPT_VERSION=1.1

# Path to EPCC hybrid benchmark executable
#EPCC_BIN=/home/tjn/projects/ompi-ecp/apps/epcc/openmpmpibench_C/mixedModeBenchmark
EPCC_BIN=./mixedModeBenchmark

# Path to MPIRUN executable
MPIRUN=mpirun

# (Optional) Any default/additional options to pass to MPIRUN
#MPIRUN_OPTIONS=" --map-by node --hostfile /home/3t4/projects/ompi-ecp/t3ompix/apps/epcc/source/openmpmpibench_C/hosts4 "
MPIRUN_OPTIONS=

usage () {
    echo "v$SCRIPT_VERSION"
    echo "Usage: $0  [-h,-d] -t OMP_NUM_THREADS -r MPI_NUM_RANKS -f INPUT_FILE"
    echo "    -t       Number of OpenMP threads"
    echo "    -r       Number of MPI ranks"
    echo "    -f       EPCC input filename"
    echo ""
    echo "Optional flags:"
    echo "    -d       Enable debug output"
    echo "    -h       Print this help info"
}

###
# MAIN
###

if [ $# -eq 0 ] ; then
    usage
    exit 1
fi

#
# Process ARGV/cmd-line
#
DEBUG=0
OPTIND=1
while getopts hdt:r:f: opt ; do
    case "$opt" in
        t)  num_threads="$OPTARG";;           # '-t' num OMP threads
        r)  num_ranks="$OPTARG";;             # '-r' num MPI ranks
        f)  input_file="$OPTARG";;            # '-f' input filename
        d)  DEBUG=1;;                         # '-d' enable debug output
        h)  print_help_info=1;;               # '-h' print help/usage info
    esac
done

shift $(($OPTIND - 1))
if [ "$1" = '--' ]; then
    shift
fi

[ $DEBUG -ne 0 ] && echo "DBG: OMP num_threads='$num_threads'"
[ $DEBUG -ne 0 ] && echo "DBG: MPI   num_ranks='$num_ranks'"
[ $DEBUG -ne 0 ] && echo "DBG: EPCC input_file='$input_file'"

if [ "$print_help_info" ]; then 
    usage 
    exit 0
fi

if [ "x$num_threads" == "x" ] ; then
    echo "Error: missing '-t' option with num threads"
    usage
    exit 1
fi

if [ "x$num_ranks" == "x" ] ; then
    echo "Error: missing '-r' option with num ranks"
    usage
    exit 1
fi

if [ "x$input_file" == "x" ] ; then
    echo "Error: missing '-f' option with EPCC input filename"
    usage
    exit 1
fi

if [ ! -f "$input_file" ] ; then
    echo "Error: missing file '$input_file'"
    usage
    exit 1
fi

#######################################################
rc=0


if [ $DEBUG -ne 0 ] ; then
    echo ""
    echo "*** DEBUG MODE (NO EXECUTE, DISPLAY ONLY) ***"
    echo "# INFO: OMP_NUM_THREADS = $num_threads"
    echo "# INFO: MPI_NUM_RANKS   = $num_ranks"
    echo "# INFO:    INPUT_FILE   = $input_file"
    echo ""
    echo "   CMD: export OMP_NUM_THREADS=$num_threads"
    echo "   CMD: $MPIRUN $MPIRUN_OPTIONS -np $num_ranks $EPCC_BIN $input_file"
    echo ""
else
    export OMP_NUM_THREADS=$num_threads
    $MPIRUN $MPIRUN_OPTIONS -np $num_ranks $EPCC_BIN $input_file
    rc=$?
fi

exit $rc


