CC=/usr/local/bin/gcc-6
CFLAGS=-Wall -Werror -O0
LIBDIRS=-I /usr/local/opt/openblas/include -L /usr/local/opt/openblas/lib
LIBS=-l blas

all: neopt opt blas

neopt: solver_neopt.c main.c utils.h
	$(CC) $(CFLAGS) $^ -o $@

opt: solver_opt.c main.c utils.h
	$(CC) $(CFLAGS) $^ -o $@

blas: solver_blas.c main.c utils.h
	$(CC) $(CFLAGS) $^ $(LIBDIRS) $(LIBS) -o $@

clean:
	rm -rf neopt opt blas out0 out1 out2 out3 out4
