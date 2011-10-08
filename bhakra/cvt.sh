# Note: I had to compile my own libgeotiff to make sure that head -c 2 final.tiff read "II"

~/bin/geotiff-bin-114-macos10.1/bin/listgeo data.tiff > info.txt
convert data.tiff +compress -endian lsb uncomp.tiff
geotifcp -g info.txt uncomp.tiff final.tiff
../../taudem/pitremove final.tiff