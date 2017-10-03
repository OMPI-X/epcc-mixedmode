#------------------------------------------------#
#                                                #
# Makefile for mixed mode Benchmark code.        #
#                                                #
# Set the CC and CFLAGS for you platform and     #
# then type 'make'.                              #
#------------------------------------------------#

CC=mpicc
CFLAGS=-fopenmp -O3 -g

#--------------------------------------#
# Build rules                          #
#--------------------------------------#

EXE=	mixedModeBenchmark

SRC=	parallelEnvironment.c \
	benchmarkSetup.c \
	output.c \
	pt_to_pt_pingpong.c \
	pt_to_pt_pingping.c \
	pt_to_pt_multiPingpong.c \
	pt_to_pt_multiPingping.c \
	pt_to_pt_haloexchange.c \
	collective_barrier.c \
	collective_broadcast.c \
	collective_scatterGather.c \
	collective_reduction.c \
	collective_alltoall.c \
	mixedModeBenchmarkDriver.c

INC=	parallelEnvironment.h \
	benchmarkSetup.h \
	output.h \
	pt_to_pt_pingpong.h \
	pt_to_pt_pingping.h \
	pt_to_pt_multiPingpong.h \
	pt_to_pt_multiPingping.h \
	pt_to_pt_haloexchange.h \
	collective_barrier.h \
	collective_broadcast.h \
	collective_scatterGather.h \
	collective_reduction.h \
	collective_alltoall.h

#
# No need to edit below this line
#

.SUFFIXES:
.SUFFIXES: .c .o

OBJ=	$(SRC:.c=.o)

.c.o:
	$(CC) $(CFLAGS) -c $<


all:	$(EXE)

$(EXE):	$(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(LFLAGS) -lpmix

$(OBJ):	$(INC) $(MF)

clean:
	rm -f $(OBJ) $(EXE) core 
