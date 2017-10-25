F90      = gfortran
CXX      = g++
FFLAG    = -Iother/qdlib/include/qd
LFLAG    = -Lother/qdlib/lib -lqd_f_main -lqdmod -lqd \
           -lgfortran -lquadmath 
#OPTIONS  = -fcheck=all -Wall -Wextra -Warray-temporaries -Wrealloc-lhs-all -pedantic -std=f2008
OPTIONS = -O3

NAME = drake
OBJ  = file_INOUT.o precision.o precision_qd.o misc.o memory.o time.o \
       eproblem.o eproblem_qd.o lproblem.o lproblem_qd.o \
       commontypes.o inputread.o systemdef.o \
       int12core.o optimizer.o \
       SCFint.o SCFdriver.o CCint2.o CCint3.o CCdriver.o calcdriver.o \
       main.o

$(NAME) : $(OBJ)
	$(CXX) -o $@ $^ $(LFLAG)

%.o : %.mod
%.o : %.f90
	$(F90) -c $(OPTIONS) $(FFLAG) $<

precision_qd.o   : precision.o
memory.o         : file_INOUT.o precision.o
time.o           : file_INOUT.o precision.o
eproblem.o       : file_INOUT.o precision.o memory.o
eproblem_qd.o    : file_INOUT.o precision.o precision_qd.o memory.o
lproblem.o       : file_INOUT.o precision.o memory.o
lproblem_qd.o    : file_INOUT.o precision.o precision_qd.o memory.o
commontypes.o    : file_INOUT.o precision.o memory.o
inputread.o      : file_INOUT.o precision.o misc.o memory.o commontypes.o
systemdef.o      : file_INOUT.o precision.o misc.o memory.o commontypes.o
int12core.o      : file_INOUT.o precision.o precision_qd.o
optimizer.o      : file_INOUT.o precision.o memory.o commontypes.o
SCFint.o         : file_INOUT.o precision.o misc.o memory.o commontypes.o \
                   int12core.o
SCFdriver.o      : file_INOUT.o precision.o memory.o eproblem.o eproblem_qd.o \
                   commontypes.o \
                   SCFint.o
CCint2.o         : file_INOUT.o precision.o misc.o memory.o commontypes.o \
                   int12core.o
CCint3.o         : file_INOUT.o precision.o memory.o commontypes.o
CCdriver.o       : file_INOUT.o precision.o memory.o lproblem.o \
                   commontypes.o CCint2.o CCint3.o
calcdriver.o     : file_INOUT.o precision.o memory.o \
                   commontypes.o optimizer.o SCFdriver.o CCdriver.o
main.o           : file_INOUT.o precision.o memory.o time.o \
                   commontypes.o inputread.o systemdef.o calcdriver.o

.PHONY : clean
clean :
	rm *.o *.mod $(NAME)
