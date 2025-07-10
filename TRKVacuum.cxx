void TRKServices::createVacuumCompositeShape()
{
  Double_t pipeRIn = 1.8f;
  Double_t A3IPLength = 1000.f;
  Double_t vacuumVesselRIn = 5.6f;
  Double_t vacuumVesselThickness = 0.08f;
  Double_t vacuumVesselLength = 76.f;
  Double_t windowLength = 250;

  // Vacuum for A and C Side
  Double_t vacuumASideLength = A3IPLength / 2. - vacuumVesselLength / 2. - windowLength;
  Double_t vacuumCSideLength = A3IPLength / 2. + vacuumVesselLength / 2.;

  // Vacuum tubes
  TGeoTube* vacuumASide = new TGeoTube("VACUUM_Ash", 0., pipeRIn, vacuumASideLength / 2.);
  TGeoTube* vacuumCSide = new TGeoTube("VACUUM_Csh", 0., vacuumVesselRIn, vacuumCSideLength / 2.);
  TGeoCone* vacuumWindow = new TGeoCone("VACUUM_WINDOW", windowLength/2., 0., vacuumVesselRIn, 0., pipeRIn);

  // Vacuum positions
  TGeoTranslation* posVacuumASide = new TGeoTranslation("VACUUM_ASIDE_POSITION", 0, 0, vacuumVesselLength / 2. + vacuumASideLength / 2. + windowLength);
  posVacuumASide->RegisterYourself();
  TGeoTranslation* posVacuumCSide = new TGeoTranslation("VACUUM_CSIDE_POSITION", 0, 0, vacuumVesselLength / 2. - vacuumCSideLength / 2.);
  posVacuumCSide->RegisterYourself();
  TGeoTranslation* posVacuumWindow = new TGeoTranslation("VACUUM_WINDOW_POSITION", 0, 0, vacuumVesselLength / 2. + windowLength/2.);
  posVacuumWindow->RegisterYourself();
  
  mVacuumCompositeFormula =
   "VACUUM_Ash:VACUUM_ASIDE_POSITION"
   "+VACUUM_Csh:VACUUM_CSIDE_POSITION"
   "+VACUUM_WINDOW:VACUUM_WINDOW_POSITION";

}
