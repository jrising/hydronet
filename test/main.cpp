#include <iostream>
#include <dims/GlobalDimensions.h>
//#include <datastr/GeographicMap.h>

using namespace std;
using namespace openworld;

int main(int argc, const char* argv[])
{
  cout << "A" << endl;
  /*GeographicMap<double>* xxx = GeographicMap<double>::loadDelimited(DividedRange(29.74583, 33.9125, 0.008333334, Inds::lat),
                                                                    DividedRange(74.77084, 85.17084, 0.008333334, Inds::lon),
                                                                    "../regress/elevation.tsv", NULL, '\t');*/
  //cout << Dims::s << endl;
  cout << GlobalDimensions::get("hi") << endl;

  cout << "Z" << endl;
}
