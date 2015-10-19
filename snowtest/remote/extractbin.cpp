#include <math.h>
#include <stdio.h>
#include <dirent.h>
#include <fstream>

#define byte unsigned char

#include <datastr/GeographicMap.h>
#include <datastr/Matrix.h>
#include <measure/Inds.h>

using namespace openworld;

typedef void (*callback)(const char* infile, Matrix<byte>& matrix, const char* outfile);

int extractBin(const char* infile, const char* outfile, callback call);
int extractDirectory(const char* directory, const char* outfile, callback call);
void callbackSave(const char* infile, Matrix<byte>& matrix, const char* outfile);
void callbackSingle(const char* infile, Matrix<byte>& matrix, const char* outfile);

float lat0;
float lon0;

int main(int argc, const char* argv[])
{
  const char* infile = argv[1];
  const char* outfile = argv[2];

  if (argc == 3) {
    // Call as extractbin directory outfile
    cout << "Extracting entire map." << endl;
    return extractBin(infile, outfile, callbackSave);
  }

  lat0 = atof(argv[3]);
  lon0 = atof(argv[4]);

  if (argc == 5) {
    // Call as extractbin directory outfile lat0 lon0
    cout << "Extracting single location." << endl;
    return extractDirectory(infile, outfile, callbackSingle);
  }

  // Call as extractbin directory outfile lat0 lon0 lat1 lon1
/*
  float lat1 = atof(argv[5]);
  float lon1 = atof(argv[6]);
*/
}

int extractBin(const char* infile, const char* outfile, callback call) {
  DividedRange latitudes = DividedRange::rounded(-90, 90, .25, Inds::lat);
  DividedRange longitudes = DividedRange::rounded(0, 180, .25, Inds::lon);
  MatrixGeographicMap<byte>* map = MatrixGeographicMap<byte>::loadBinary(latitudes, longitudes, infile, 1, NULL);
  Matrix<byte>& matrix = map->getValues();
  call(infile, matrix, outfile);
  delete map;

  return 0;
}

int extractDirectory(const char* directory, const char* outfile, callback call) {
  // Iterate through files in directory
  DIR *dir;
  struct dirent *ent;
  if ((dir = opendir(directory)) != NULL) {
    while ((ent = readdir(dir)) != NULL) {
      cout << "Checking " << ent->d_name << endl;
      if (strlen(ent->d_name) <= 2)
        continue;

      string fullpath = directory + (string) ent->d_name;
      try {
        extractBin(fullpath.c_str(), outfile, call);
      } catch (exception& e) {
        ofstream outfp;
        outfp.open(outfile, ios_base::app);
        outfp << "NA" << endl;
        outfp.close();
      }
    }

    closedir(dir);
  } else {
    // could not open directory
    perror("Could not open directory.");
    return EXIT_FAILURE;
  }

  return 0;
}

void callbackSave(const char* infile, Matrix<byte>& matrix, const char* outfile) {
  matrix.saveDelimited(outfile, FileFormatter<byte>::formatSimple, ',');
}

void callbackSingle(const char* infile, Matrix<byte>& matrix, const char* outfile) {
  unsigned row = (unsigned) (360.0 + 360.0 * ((90.0 - lat0) / 90.0) * cos(M_PI * lon0 / 180.0));
  unsigned col = (unsigned) (360.0 + 360.0 * ((90.0 - lat0) / 90.0) * sin(M_PI * lon0 / 180.0));
  cout << row << ", " << col << endl;

  ofstream outfp;
  outfp.open(outfile, ios_base::app);
  outfp << infile << "," << matrix.getDouble(row, col) << endl;
  outfp.close();
}
