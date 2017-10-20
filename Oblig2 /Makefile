CC = mpicc
CFLAGS = -O2 -fopenmp -lm
COMFLAGS = "small_matrix_a.bin" "small_matrix_b.bin" "out_matrix.bin"

PROJ = oblig2
OBJS = oblig2.o

all : oblig2 $(PROJ)

oblig2: $(OBJS)
	$(CC) $(CFLAGS) $< -o $@

run: oblig2
	mpirun -n 5 ./oblig2 $(COMFLAGS)

obj-clean:
	$(RM) *.o

exec-clean:
	$(RM) $(PROJ)

autosave-clean:
	$(RM) *~

clean:
	$(MAKE) obj-clean
	$(MAKE) exec-clean
	$(MAKE) autosave-clean
