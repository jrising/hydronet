SYS = Linux

CC = mpic++
DEBUG = -g
OWD = ~/projects/openworld
CFLAGS = -Wall -c $(DEBUG) -I $(OWD) -O3 
LFLAGS = -Wall $(DEBUG) -O3 

OBJS_MAIN = basinmask.o $(OWD)/geotiff/tiffIO.o $(OWD)/geotiff/commonLib.o $(OWD)/datastr/GeographicMap.o

all: $(OBJS_MAIN)
	$(CC) $(LFLAGS) $(OBJS_MAIN) -o basinmask

#Inference rule - states a general rule for compiling .o files
%.o : %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	\rm $(OBJS_MAIN)
	\rm *.o basinmask
