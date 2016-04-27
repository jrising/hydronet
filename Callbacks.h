#include "SJHydroNetModel.h"
#include <exception>

class SJHydroNetModelSaveAllCells : public SJHydroNetModelStepCallback {
 protected:
  string allcells;
  
 public:
  SJHydroNetModelSaveAllCells(string allcells) {
    this->allcells = allcells;
  }

  virtual void setup(SJHydroNetModel& model) {
    ofstream cellfile;
    cellfile.open(allcells.c_str(), ios::out);
    for (unsigned rr = 1; rr < model.getLatitudes().count() - 1; rr++)
      for (unsigned cc = 1; cc < model.getLongitudes().count() - 1; cc++)
	cellfile << model.getLatitudes().getCellCenter(rr).getValue() << "\t" << model.getLongitudes().getCellCenter(cc).getValue() << "\t";
    cellfile << endl;
    cellfile.close();
  }
  
  virtual void post(SJHydroNetModel& model) {
    ofstream cellfile;
    cellfile.open(allcells.c_str(), ios::out | ios::app);
    for (unsigned rr = 1; rr < model.getLatitudes().count() - 1; rr++)
      for (unsigned cc = 1; cc < model.getLongitudes().count() - 1; cc++) {
        double precipVolume, meltVolume;
        model.getHydroNet().sumNodeMapVolumes(rr, cc, precipVolume, meltVolume);
      
        cellfile << precipVolume << "\t" << meltVolume << "\t";
      }
    cellfile << endl;
    cellfile.close();
  }
};

class SJHydroNetModelSaveSomeCells : public SJHydroNetModelStepCallback {
 protected:
  string filename;
  vector<double> latitudes;
  vector<double> longitudes;
  
 public:
  SJHydroNetModelSaveSomeCells(string filename) {
    this->filename = filename;
  }

  void addLocation(double lat, double lon) {
    latitudes.push_back(lat);
    longitudes.push_back(lon);
  }
  
  virtual void setup(SJHydroNetModel& model) {
    ofstream fp;
    fp.open(filename.c_str(), ios::out);
    fp << "time,location,variable,value" << endl;
    for (unsigned ii = 0; ii < latitudes.size(); ii++) {
      fp << "NA," << ii + 1 << ",latitude," << latitudes[ii] << endl;
      fp << "NA," << ii + 1 << ",longitude," << longitudes[ii] << endl;
    }
    fp.close();
  }
  
  virtual void post(SJHydroNetModel& model) {
    ofstream fp;
    fp.open(filename.c_str(), ios::out | ios::app);
    for (unsigned ii = 0; ii < latitudes.size(); ii++) {
      int rr = model.getLatitudes().inRange(latitudes[ii]);
      int cc = model.getLongitudes().inRange(longitudes[ii]);
      double precipVolume, meltVolume;
      model.getHydroNet().sumNodeMapVolumes(rr, cc, precipVolume, meltVolume);

      fp << model.getTime() << "," << ii + 1 << ",precip," << precipVolume << endl;
      fp << model.getTime() << "," << ii + 1 << ",melt," << meltVolume << endl;
      try {
	fp << model.getTime() << "," << ii + 1 << ",volume," << model.getSnowModel().getVolumes().getCell(rr, cc) << endl;
      } catch (exception& ex) {}
    }
    fp.close();
  }
};

class SJHydroNetModelStoreSingleVolume : public SJHydroNetModelStepCallback {
 protected:
  double latitude;
  double longitude;
  vector<time_t> times;
  vector<double> volumes;
  
 public:
  SJHydroNetModelStoreSingleVolume(double latitude, double longitude) {
    this->latitude = latitude;
    this->longitude = longitude;
  }

  vector<time_t> getTimes() {
    return times;
  }

  vector<double> getVolumes() {
    return volumes;
  }
  
  virtual void setup(SJHydroNetModel& model) {
    times.clear();
    volumes.clear();
  }
  
  virtual void post(SJHydroNetModel& model) {
    int rr = model.getLatitudes().inRange(latitude);
    int cc = model.getLongitudes().inRange(longitude);
    double volume = 0;
    try {
      volume = model.getSnowModel().getVolumes().getCell(rr, cc);
    } catch (exception& ex) {}

    volumes.push_back(volume);
    times.push_back(model.getTime());
  }
};
