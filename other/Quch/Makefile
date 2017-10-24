COMPILER = /home/mitek/programy/gcc-7.1.0/install
F90      = $(COMPILER)/bin/gfortran
CXX      = $(COMPILER)/bin/g++
FFLAG    = -I../qdlib/include/qd
LFLAG    = -L../qdlib/lib -lqdmod -lqd \
           -L$(COMPILER)/lib64 -lgfortran -lquadmath \
           -Wl,-rpath -Wl,$(COMPILER)/lib64
#OPTIONS = -fcheck=all -Wall -Wextra -Warray-temporaries -Wrealloc-lhs-all -Wno-compare-reals -pedantic -std=f2008
OPTIONS = -O3

NAME = threel
OBJ  = indices.o korobov.o pachucki.o f_main.o c_main.o

$(NAME) : $(OBJ)
	$(CXX) -o $@ $^ $(LFLAG)

%.o : %.mod
%.o : %.f90
	$(F90) -c $(OPTIONS) $(FFLAG) $<
%.o : %.cpp
	$(CXX) -c $<

pachucki.o : indices.o korobov.o
f_main.o   : indices.o pachucki.o
c_main.o   : f_main.o

.PHONY : clean
clean :
	rm *.o *.mod $(NAME)
