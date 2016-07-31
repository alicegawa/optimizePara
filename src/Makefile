#Makefile for cmaes_mpi

TARGET = \
	estimate_main \

#CC = mpifccpx
CC = mpicc
#CFLAGS = -Kfast
CFLAGS = -O3
LDFLAGS = -lm
INCLUDES = -I.
LIBS = -lm
OBJS = cmaes.o my_boundary_transformation.o boundary_transformation.o
TOBJS = $(addsuffix .o, $(TARGET))

all: $(TARGET)

$(TARGET): $(TOBJS) $(OBJS)
	$(CC) $(CFLAGS) -o ./$@ $@.o $(OBJS) $(LDFLAGS)

.c.o: 
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) -c $<

cmaes_interface.h: cmaes.h
$(TOBGS): cmaes_interface.h my_boundary_transformation.h
my_boundary_transformation.h: boundary_transformation.h
my_boundary_transformation.o: my_boundary_transformation.h
boundary_transformation.o: boundary_transformation.h
cmaes.o: cmaes.h

clean:
	rm *.o
