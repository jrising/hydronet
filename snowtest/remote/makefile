SYS = Linux

DEFS = -DSKIP_TIFF -DSKIP_GETLINE -DSKIP_GZSTREAM
CC = g++
DEBUG = -g
OWD = ~/projects/openworld
CFLAGS = -Wall -c $(DEBUG) -I $(OWD) -O3 $(DEFS) -Wgnu-static-float-init
# others  -I $(OWD)/ext/gzstream
LFLAGS = -Wall $(DEBUG) -O3 -L. -lz
# others: -lgzstream

OBJS_BASE = $(OWD)/datastr/GeographicMap.o $(OWD)/memory/Transients.o $(OWD)/dims/Dims.o $(OWD)/dims/Dimensions.o $(OWD)/dims/Dimensionless.o $(OWD)/measure/Units.o $(OWD)/measure/Inds.o $(OWD)/dims/GlobalDimensions.o
OBJS_MAIN = ../SJHydroNetModel.o $(OWD)/geotiff/tiffIO.o $(OWD)/geotiff/commonLib.o $(OBJS_BASE) $(OWD)/utils/Timer.o $(OWD)/ext/gzstream/gzstream.o

all: main.o $(OBJS_MAIN) $(OWD)/datastr/Matrix.h $(OWD)/datastr/GeographicMap.h $(OWD)/datastr/PartialConfidenceTemporalGeographicMap.h
	$(CC) $(LFLAGS) main.o $(OBJS_MAIN) -o model

extractbin: extractbin.o $(OBJS_BASE) $(OWD)/datastr/GeographicMap.h
	$(CC) $(LFLAGS) extractbin.o $(OBJS_BASE) -o extractbin

extractpoint: extractpoint.o $(OBJS_BASE) $(OWD)/datastr/GeographicMap.h
	$(CC) $(LFLAGS) extractpoint.o $(OBJS_BASE) -o extractpoint

#Inference rule - states a general rule for compiling .o files
%.o : %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	\rm $(OBJS_MAIN)
	\rm *.o model
