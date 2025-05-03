CC=			gcc
CFLAGS=		-g -Wall -Wc++-compat -O3
CPPFLAGS=
INCLUDES=
OBJS=		reader.o misc.o
PROG=		minisplice
LIBS=		-lz -lm

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address -ldl
endif

.SUFFIXES:.c .o
.PHONY:all clean depend

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

minisplice:$(OBJS) main.o
		$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

clean:
		rm -fr *.o a.out $(PROG) *~ *.a *.dSYM

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(CPPFLAGS) -- *.c)

# DO NOT DELETE

main.o: ketopt.h
misc.o: minisp.h
reader.o: msppriv.h minisp.h kseq.h
