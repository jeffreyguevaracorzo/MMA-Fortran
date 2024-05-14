# Makefile compilation parameters
CC = gfortran
CFLAGS = -O1 -fbounds-check -fbacktrace -fcheck=all -g -Wall -Wextra -Wrealloc-lhs-all
TARGET = Execute
EXT = -llapack -lblas
FILES = MMA_Routines.f90 \
		MMA_Variables.f90 \
		MMA_Main.f90

# object files
OBJ1 = ${FILES:.f90=.o}

# compilation and cleanup command
%.o : %.f90
	${CC} ${CFLAGS} -o $@ -c $<
	
${TARGET} : ${OBJ1}
	${CC} ${CFLAGS} -o $@ ${OBJ1} $(EXT)

.PHONY : clean
clean :
	@rm -f *.o *.mod ${TARGET} ${OBJ1}