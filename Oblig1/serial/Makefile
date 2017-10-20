CC = gcc
CFLAGS = -O2
LDFLAGS = -L../simple-jpeg/ -lsimplejpeg
VFLAGS =  --leak-check=full --track-origins=yes -v
COMFLAGS = 0.2 50 "mona_lisa_noisy.jpeg" "out.jpeg"

PROJ = serial_main
OBJS = serial_main.o

all : simple-jpeg $(PROJ)

serial_main : $(OBJS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

simple-jpeg :
	$(MAKE) -C ../simple-jpeg

run: serial_main
	./serial_main $(COMFLAGS)

test : serial_main
		valgrind $(VFLAGS) ./serial_main $(COMFLAGS)

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
