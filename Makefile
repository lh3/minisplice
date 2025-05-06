CC=			gcc
CFLAGS=		-g -std=c99 -Wall -Wno-unused-function -Wc++-compat -O3
CPPFLAGS=	-DHAVE_PTHREAD
INCLUDES=
OBJS=		kautodiff.o kann.o reader.o misc.o strmap.o bed.o train.o apply.o
PROG=		minisplice
LIBS=		-lpthread -lz -lm

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

apply.o: kann.h kautodiff.h msppriv.h minisp.h ketopt.h
bed.o: msppriv.h minisp.h ksort.h
kann.o: kann.h kautodiff.h
kautodiff.o: kautodiff.h
main.o: msppriv.h minisp.h ketopt.h
misc.o: msppriv.h minisp.h
reader.o: msppriv.h minisp.h kseq.h
strmap.o: msppriv.h minisp.h khashl.h
train.o: kann.h kautodiff.h msppriv.h minisp.h ketopt.h
