SYS = Linux

CC = mpic++
DEBUG = -g
OWD = ~/projects/openworld
CFLAGS = -Wall -c $(DEBUG) -I $(OWD)
LFLAGS = -Wall $(DEBUG)

OBJS_MAIN = $(OWD)/geotiff/tiffIO.o $(OWD)/geotiff/commonLib.o $(OWD)/datastr/GeographicMap.o $(OWD)/memory/Transients.o $(OWD)/dims/Dims.o $(OWD)/dims/Dimensions.o $(OWD)/dims/Dimensionless.o $(OWD)/measure/Units.o $(OWD)/measure/Inds.o $(OWD)/dims/GlobalDimensions.o

all: $(OBJS_MAIN) basinmask.o
	$(CC) $(LFLAGS) basinmask.o $(OBJS_MAIN) -o basinmask

#Inference rule - states a general rule for compiling .o files
%.o : %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	\rm $(OBJS_MAIN)
	\rm *.o basinmask
