#ifndef HYDROPOINTMODEL_H
#define HYDROPOINTMODEL_H

using namespace openworld;

class SJHydroPointModel : public SJHydroModel<double, double, bool> {
 protected:
  // Inputs
  PartialConfidenceTimeSeries<double>* precipitation;
  PartialConfidenceTimeSeries<double>* surfaceTemp;
  SnowPointModel* snowModel;
  double elevation;
  double mmdayToVolume;

  list<double> outFlowRain;
  list<double> outFlowMelt;
  Measure now;

 public:
  SJHydroPointModel(Indicator timeind, double meltDegreeDayFactor, double meltDegreeDaySlope, double rainRunoffCoefficient, double meltRunoffCoefficient, double groundCoefficient, double groundToBaseflowDay, double rainOnSnowCoefficient, double surfaceEvaporationFactor, double riverEvaporationFactor, doubl mmdayToVolume);
  SJHydroPointModel(SJHydroPointModel& copy);

  ~SJHydroPointModel();

  void setPrecipitation(PartialConfidenceTimeSeries<double>* precipitation);
  void setTemperature(PartialConfidenceTimeSeries<double>* surfaceTemp);
  void setSnowModel(SnowPointModel* snowModel);
  void setElevation(double elevation);

  void stepDay();
};

#endif
