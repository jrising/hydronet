SYS = Linux

CC = mpic++
DEBUG = -g
OWD = ~/projects/openworld
CFLAGS = -Wall -c $(DEBUG) -I $(OWD) -O3 -I $(OWD)/ext/gzstream
LFLAGS = -Wall $(DEBUG) -O3 -L. -lgzstream -lz

OBJS_MAIN = ../SJHydroNetModel.o $(OWD)/geotiff/tiffIO.o $(OWD)/geotiff/commonLib.o $(OWD)/datastr/GeographicMap.o $(OWD)/memory/Transients.o $(OWD)/dims/Dims.o $(OWD)/dims/Dimensions.o $(OWD)/dims/Dimensionless.o $(OWD)/measure/Units.o $(OWD)/measure/Inds.o $(OWD)/dims/GlobalDimensions.o $(OWD)/utils/Timer.o $(OWD)/ext/gzstream/gzstream.o

all: main.o $(OBJS_MAIN) $(OWD)/datastr/Matrix.h $(OWD)/datastr/GeographicMap.h $(OWD)/datastr/PartialConfidenceTemporalGeographicMap.h
	$(CC) $(LFLAGS) main.o $(OBJS_MAIN) -o model

cal: calibrate.o BhakraTest.h $(OBJS_MAIN) $(OWD)/datastr/Matrix.h $(OWD)/datastr/GeographicMap.h $(OWD)/datastr/PartialConfidenceTemporalGeographicMap.h
	$(CC) calibrate.o $(OBJS_MAIN) $(LFLAGS) -o model

load: testload.o BhakraTest.h $(OBJS_MAIN) $(OWD)/datastr/Matrix.h $(OWD)/datastr/GeographicMap.h $(OWD)/datastr/PartialConfidenceTemporalGeographicMap.h
	$(CC) testload.o $(OBJS_MAIN) $(LFLAGS) -o testload

future: future.o $(OBJS_MAIN) $(OWD)/datastr/Matrix.h $(OWD)/datastr/GeographicMap.h $(OWD)/datastr/PartialConfidenceTemporalGeographicMap.h
	$(CC) $(LFLAGS) future.o $(OBJS_MAIN) -o model

#Inference rule - states a general rule for compiling .o files
%.o : %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	\rm $(OBJS_MAIN)
	\rm *.o model
