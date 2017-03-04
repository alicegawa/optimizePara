TARGET = estimate_main
CC = mpicc
CFLAGS = -O3 -Wall
#CC = mpifccpx
#CFLAGS = -Kfast
LDFLAGS = -lm
INCLUDES = -I.
LIBS = -lm
OBJS = my_boundary_transformation.o boundary_transformation.o mplmcma.o
TOBJS = $(addsuffix .o, $(TARGET))

all: $(TARGET)

$(TARGET): $(TOBJS) $(OBJS)
	$(CC) $(CFLAGS) -o ./$@ $@.o $(OBJS) $(LDFLAGS)
.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) -c $< $(LDFLAGS)
$(TOBGS): my_boundary_transformation.h mplmcma.h
my_boundary_transformation.h: boundary_transformation.h
my_boundary_transformation.o: my_boundary_transformation.h
boundary_transformation.o: boundary_transformation.h
mplmcma.o: mplmcma.h

clean:
	rm *.o
