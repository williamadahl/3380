CC = mpicc
CFLAGS = -O2 -Wall -Wextra -g
LDFLAGS = -L../simple-jpeg/ -lsimplejpeg
VFLAGS =  --leak-check=full --track-origins=yes -v --show-leak-kinds=all
COMFLAGS = 0.2 100 "mona_lisa_noisy.jpeg" "out.jpeg"

PROJ = parallel_main
OBJS = parallel_main.o

all : simple-jpeg $(PROJ)

parallel_main : $(OBJS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

simple-jpeg :
	$(MAKE) -C ../simple-jpeg

run: parallel_main
	mpirun -n 18 ./parallel_main $(COMFLAGS)

test : parallel_main
	valgrind $(VFLAGS) ./parallel_main $(COMFLAGS)

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
