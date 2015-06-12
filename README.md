hydronet
========

Hydrological Network Model

Compiling
---------

The master makefile is currently the one within the bhakra subdirectory.  So you can create executables using:

```
cd bhakra
make
```

basinmask
---------

`basinmask` is a simple tool for constructing a basin mask for a given
drainage location from a TauDEM D-8 or D-infinity GeoTIFF.

The makefile for basinmask is currently in the root directory:

```
make basinmask
```

Then call basinmask as follows:

```
./basinmask <ang.tiff> <mask.tsv> <OUTLAT> <OUTLON> <LAT0> <LAT1> <DLAT> <LON0> <LON1> <DLON>
```

Where `<ang.tiff>` is a D-8 or D-infinity angle file from TauDEM; `<mask.tsv>` is the output file (a TSV file); `<OUTLAT>, <OUTLON>` is the outlet flow location; `<LAT0>, <LON0>, <LAT1>, <LON1>` describe the bounding box of the DEM file, and `<DLAT>` and `<DLON>` is the resolution of the DEM in degrees of latitude and longitude, respectively.
