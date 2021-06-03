# Makefile for ...su/main


##Compiler on my PC: gcc
CC=gcc

#Libaries from my PC
INCPAR=-I/home/zhouw/Codes/GitHub/SeisUnix/include/
LIBPAR=-L/home/zhouw/Codes/GitHub/SeisUnix/lib -lsu -lpar -lcwp -lm



LINEHDRFLAG =
XDRFLAG = # -DSUXDR
ENDIANFLAG = -DCWP_LITTLE_ENDIAN
LARGE_FILE_FLAG = -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE


CWP_FLAGS = -std=c99 -openmp  $(LARGE_FILE_FLAG) $(ENDIANFLAG) $(XDRFLAG) $(LINEHDRFLAG)


FUN = main


.PHONY: all

all : $(patsubst %.c,%,$(FUN))

%: %.c
	$(CC) $< -O2 $(CWP_FLAGS) $(INCPAR) $(LIBPAR) -o ./$@

clean:
	rm -f main a.out *.mod junk* JUNK* core
	rm -f in out
