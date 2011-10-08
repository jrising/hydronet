SYS = Linux

CC = mpic++
DEBUG = -g
OWD = ~/projects/openworld
CFLAGS = -Wall -c $(DEBUG) -I $(OWD)
LFLAGS = -Wall $(DEBUG)

OBJS_MAIN = main.o HydroModel.o $(OWD)/geotiff/tiffIO.o $(OWD)/geotiff/commonLib.o $(OWD)/datastr/GeographicMap.o $(OWD)/memory/Transients.o $(OWD)/dims/Dims.o $(OWD)/indicator/Inds.o $(OWD)/dims/GlobalDimensions.o

all: $(OBJS_MAIN) $(OWD)/datastr/Matrix.h $(OWD)/datastr/GeographicMap.h $(OWD)/datastr/PartialConfidenceTemporalGeographicMap.h
	$(CC) $(LFLAGS) $(OBJS_MAIN) -o model

#Inference rule - states a general rule for compiling .o files
%.o : %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	\rm $(OBJS_MAIN)
	\rm *.o model
