CC = gcc
OPT= -O2 -lm -std=c99 -g -fopenmp
#OPT+= -DLONGIDS
OPTIONS = $(OPT)
GSL_LIB = $(GSL_LIBDIR)
GSL_INCL =$(GSL_INCDIR)

#INCL = lmmin.c lm_eval.c libgad.c

CFLAGS = $(OPTIONS)
DEPS= libgad.h
OBJS = libgad.o lmmin.o lm_eval.o KdTree.o ompfuncs.o OctTree.o

SOURCES=coarsen.c cut_subfind.c cutgadget.c friends2idlist.c get_group_catalogue.c getboundvol.c id2ascii.c pos2ascii.c id2pos.c join_gadfiles.c mk_id_list.c nfw.c phead.c trace.c trhalo.c cutsphere.c
EXECS=$(SOURCES:.c=)


all: makefile $(EXECS)

bin:
	mkdir bin

%.o: %.c $(DEPS)
	$(CC) $(OPTIONS) -c -o $@ $<

$(EXECS): $(OBJS) bin 
	$(CC) $(OPTIONS) -I$(GSL_INCDIR) -L$(GSL_INCDIR) $(OBJS) $@.c -lgsl -lgslcblas -lm -o bin/$@

clean:
	rm -f $(OBJS) $(EXEC) bin/*
