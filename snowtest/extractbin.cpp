#include <stdio.h>
#include <dirent.h>
#include "../datastr/Matrix.h"
#include "../datastr/FileFormats.h"

int main(int argc, const char* argv[])
{
  // Call as extractbin directory lat0 lon0 lat1 lon1
  char* directory = argv[1];
  float lat0 = atof(argv[2]);
  float lon0 = atof(argv[3]);
  float lat1 = atof(argv[4]);
  float lon1 = atof(argv[5]);

  // Iterate through files in directory
  DIR *dir;
  struct dirent *ent;
  if ((dir = opendir(directory)) != NULL) {
    while ((ent = readdir(dir)) != NULL) {
      extractBin(ent->d_name);
    }
    closedir(dir);
  } else {
    /* could not open directory */
    perror("Could not open directory.");
    return EXIT_FAILURE;
  }
}
