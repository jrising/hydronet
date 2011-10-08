SYS = Linux

CC = mpic++
DEBUG = -g
OWD = ~/projects/openworld
CFLAGS = -Wall -c $(DEBUG) -I $(OWD)
LFLAGS = -Wall $(DEBUG)

OBJS_MAIN = main.o $(OWD)/dims/Dims.o $(OWD)/dims/GlobalDimensions.o

all: $(OBJS_MAIN)
	$(CC) $(LFLAGS) $(OBJS_MAIN) -o model

#Inference rule - states a general rule for compiling .o files
%.o : %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	\rm $(OBJS_MAIN)
	\rm *.o
