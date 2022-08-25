CC=gcc
CFLAGS=-I. -lm -lfftw3 -g -O2 -I/usr/lib/jvm/java-8-openjdk-amd64/include/ -I/usr/lib/jvm/java-8-openjdk-amd64/include/linux/
DEPS = mosse-eqs-fftC.h quasse_eqs-fftC.h rfftw.h
OBJ = main.o
LFLAGS = -shared -fPIC

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

libtest.so: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LFLAGS)


.PHONY: clean

clean:
	rm -f *.o main