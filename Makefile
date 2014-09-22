CC     = gcc
CFLAGS = -c -O2
LFLAGS = -L. -lm -lRNA -o
OBJ    = fold_vars.o read_epars.o energy_par.o params.o utils.o foldMod.o get_barrier.o
EXE    = get_barrier

all: get_barrier

get_barrier: get_barrier.o fold_vars.o foldMod.o utils.o energy_par.o read_epars.o params.o 
	$(CC) get_barrier.o fold_vars.o foldMod.o utils.o energy_par.o \
	read_epars.o params.o $(LFLAGS) get_barrier

fold_vars.o: fold_vars.c fold_vars.h
	$(CC) $(CFLAGS) fold_vars.c

get_barrier.o: get_barrier.c get_barrier.h utils.h
	$(CC) $(CFLAGS) get_barrier.c

foldMod.o: foldMod.c utils.h energy_par.h fold_vars.h pair_mat.h params.h
	$(CC) $(CFLAGS) foldMod.c

utils.o: utils.c config.h
	$(CC) $(CFLAGS) utils.c

energy_par.o: energy_par.c energy_const.h intloops.h
	$(CC) $(CFLAGS) energy_par.c

clean:	
	rm -f *~ $(OBJ)

small:
	make clean
	rm $(EXE)
