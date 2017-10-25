F90      = gfortran
CXX      = g++
FFLAG    = -I../qdlib/include/qd
LFLAG    = -L../qdlib/lib -lqdmod -lqd \
           -lgfortran -lquadmath 
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
