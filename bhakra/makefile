SYS = Linux

UNAME := $(shell uname)

ifeq ($(UNAME), Darwin)
	OSARGS = -DSKIP_GETLINE
else
	OSARGS = ""
endif

CC = mpic++
DEBUG = -g
OWD = ~/projects/openworld
ODM = ~/projects/opendims
CFLAGS = -Wall -c $(DEBUG) -I $(OWD) -I $(ODM) -O3 -I $(OWD)/ext/gzstream $(OSARGS)
LFLAGS = -Wall $(DEBUG) -O3 -L. -L $(OWD)/ext/gzstream -lgzstream -lz

OBJS_MAIN = ../SJHydroNetModel.o ../HydroNet.o $(OWD)/geotiff/tiffIO.o $(OWD)/geotiff/commonLib.o $(OWD)/datastr/GeographicMap.o $(OWD)/memory/Transients.o $(ODM)/dims/Dims.o $(ODM)/dims/Dimensions.o $(ODM)/dims/Dimensionless.o $(ODM)/measure/Units.o $(ODM)/measure/Inds.o $(ODM)/dims/GlobalDimensions.o $(OWD)/utils/Timer.o $(OWD)/ext/gzstream/gzstream.o

all: main.o $(OBJS_MAIN) $(OWD)/datastr/Matrix.h $(OWD)/datastr/GeographicMap.h $(OWD)/datastr/PartialConfidenceTemporalGeographicMap.h
	$(CC) main.o $(OBJS_MAIN) -o model $(LFLAGS)

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
