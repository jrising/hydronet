SYS = Linux

CC = mpic++
DEBUG = -g
OWD = ~/projects/openworld
ODS = ~/projects/opendims
CFLAGS = -std=c++11 -Wall -c $(DEBUG) -I $(OWD) -I $(ODS)
LFLAGS = -Wall $(DEBUG)

OBJS_MAIN = $(OWD)/geotiff/tiffIO.o $(OWD)/geotiff/commonLib.o $(OWD)/datastr/GeographicMap.o $(OWD)/memory/Transients.o $(ODS)/dims/Dims.o $(ODS)/dims/Dimensions.o $(ODS)/dims/Dimensionless.o $(ODS)/measure/Units.o $(ODS)/measure/Inds.o $(ODS)/dims/GlobalDimensions.o

all: $(OBJS_MAIN) basinmask.o
	$(CC) $(LFLAGS) basinmask.o $(OBJS_MAIN) -o basinmask

#Inference rule - states a general rule for compiling .o files
%.o : %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	\rm $(OBJS_MAIN)
	\rm *.o basinmask
