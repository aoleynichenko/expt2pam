CC=clang
CFLAGS=-g
OBJECTS=main.o basis.o expt_parse.o lex.yy.o error.o bfn.o basis_lib.o \
	elements.o molecule.o ecp.o ecp_lib.o

all: expt2pam.x

clean:
	rm -rf *.o *.x

expt2pam.x: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o expt2pam.x

main.o: main.c
	$(CC) $(CFLAGS) main.c -c

basis.o: basis.c
	$(CC) $(CFLAGS) basis.c -c

expt_parse.o: expt_parse.c
	$(CC) $(CFLAGS) expt_parse.c -c

lex.yy.o: lex.yy.c
	$(CC) $(CFLAGS) lex.yy.c -c

error.o: error.c
	$(CC) $(CFLAGS) error.c -c

bfn.o: bfn.c
	$(CC) $(CFLAGS) bfn.c -c

basis_lib.o: basis_lib.c
	$(CC) $(CFLAGS) basis_lib.c -c

elements.o: elements.c
	$(CC) $(CFLAGS) elements.c -c

molecule.o: molecule.c
	$(CC) $(CFLAGS) molecule.c -c

ecp.o: ecp.c
	$(CC) $(CFLAGS) ecp.c -c

ecp_lib.o: ecp_lib.c
	$(CC) $(CFLAGS) ecp_lib.c -c
