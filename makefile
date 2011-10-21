SYS = Linux

CC = mpic++
DEBUG = -g
OWD = ~/projects/openworld
CFLAGS = -Wall -c $(DEBUG) -I $(OWD)
LFLAGS = -Wall $(DEBUG)

OBJS_MAIN = $(OWD)/geotiff/tiffIO.o $(OWD)/geotiff/commonLib.o $(OWD)/datastr/GeographicMap.o $(OWD)/memory/Transients.o $(OWD)/dims/Dims.o $(OWD)/indicator/Inds.o $(OWD)/dims/GlobalDimensions.o

all: basinmask.o basinscale.o demscale.o $(OBJS_MAIN)
	$(CC) $(LFLAGS) basinmask.o $(OBJS_MAIN) -o basinmask
	$(CC) $(LFLAGS) basinscale.o $(OBJS_MAIN) -o basinscale
	$(CC) $(LFLAGS) demscale.o $(OBJS_MAIN) -o demscale

#Inference rule - states a general rule for compiling .o files
%.o : %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	\rm $(OBJS_MAIN)
	\rm *.o basinmask
