
#include <string>
#include <unordered_set>

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include "FCTSimulation/Detector.h"
#include "FCTBase/FCTBaseParam.h"
#include "ITSMFTSimulation/Hit.h"
#include "SimulationDataFormat/BaseHits.h"
#include "Steer/MCKinematicsReader.h"
#endif

#include "TCanvas.h"
#include "TEveGeoNode.h"
#include "TEveGeoShape.h"
#include "TEveGeoShapeExtract.h"
#include "TEveManager.h"
#include "TEvePointSet.h"
#include "TEveTrack.h"
#include "TFile.h"
#include "TGeoManager.h"
#include "TGeoNode.h"
#include "TGeoVolume.h"
#include "TGLOrthoCamera.h"
#include "TGLViewer.h"
#include "TPad.h"
#include "TStyle.h"

#include <fmt/format.h>

TEveTrackPropagator* g_prop = 0;

struct Hit {
  int32_t event;
  int32_t track;
  int32_t pdg;
  float x;
  float y;
  float z;
  float px;
  float py;
  float pz;
};

void display_with_tracks()
{
  //__________________________________________
  // event manager
  fmt::print("Creating TEveManager\n");
  TEveManager::Create();

  //__________________________________________
  // geometry
  fmt::print("Importing geometry\n");
  gStyle->SetCanvasPreferGL(true);
  gGeoManager->Import("o2sim_geometry.root");

  auto* topv = gGeoManager->GetTopVolume();
  fmt::print("Top volume:{}\n", topv->GetName());

  gGeoManager->SetVisLevel(4);

  int numberOfLayers = 9;

  auto layerColor = kBlack;

  for (int32_t i = 0; i < numberOfLayers; ++i) {
    std::string layerName = std::string("FCTSensor") + std::string("_") + std::to_string(i);
    TGeoVolume* volLay = gGeoManager->GetVolume(layerName.c_str());
    volLay->SetTransparency(70);
    volLay->SetFillColor(layerColor);
    volLay->SetLineColor(layerColor);
  }

  TGeoNode* tnode = gGeoManager->GetVolume("barrel")->FindNode("FCTV_2");
  fmt::print("Found FCTV node:{}\n", tnode->GetName());

  TEveGeoTopNode* evenode = new TEveGeoTopNode(gGeoManager, tnode);
  evenode->SetVisLevel(3);
  gEve->AddGlobalElement(evenode);
  fmt::print("Added FCT to viewer\n");

  //__________________________________________
  // reading and drawing hits
  fmt::print("Reading FCT hits\n");

  std::vector<Hit> points;

  std::ifstream rfile;
  rfile.open("hits.txt");
  std::string line;
  while (std::getline(rfile, line)) {
    std::istringstream iss(line);
    int32_t event, track, pdg;
    float x, y, z, px, py, pz;
    Hit hit{};
    iss >> hit.event >> hit.track >> hit.pdg >>
      hit.x >> hit.y >> hit.z >>
      hit.px >> hit.py >> hit.pz;
    points.emplace_back(hit);
  }

  auto markerStyle = kFullCircle;
  auto markerColor = kRed;
  float markerSize = 1.f;

  int32_t prEvent = points[0].event;
  int32_t prTrack = points[0].track;
  std::string s = fmt::format("event_{}_hits_{}", prEvent, prTrack);
  TEvePointSet* evePoints = new TEvePointSet(s.data());
  for (const auto& p : points) {
    int32_t event = p.event;
    int32_t track = p.track;
    if (event != prEvent || track != prTrack) {
      evePoints->SetMarkerColor(markerColor);
      evePoints->SetMarkerStyle(markerStyle);
      evePoints->SetMarkerSize(markerSize);
      gEve->AddElement(evePoints, 0);
      prEvent = event;
      prTrack = track;
      std::string s = fmt::format("event_{}_hits_{}", event, track);
      evePoints = new TEvePointSet(s.data());
    }
    evePoints->SetNextPoint(p.x, p.y, p.z);
  }
  evePoints->SetMarkerColor(markerColor);
  evePoints->SetMarkerStyle(markerStyle);
  evePoints->SetMarkerSize(markerSize);
  gEve->AddElement(evePoints, 0);

  //__________________________________________
  // setup view

  gEve->GetDefaultGLViewer()->UseLightColorSet();
  gEve->GetDefaultGLViewer()->UpdateScene();

  double c[3] = {0, 0, -466};
  gEve->GetDefaultGLViewer()->CurrentCamera().Configure(9., 800., c, -20. / 180. * 3.14, 150. / 180. * 3.14);

  gEve->GetDefaultGLViewer()->SetGuideState(TGLUtil::kAxesEdge, true, false, nullptr);
  gEve->GetDefaultGLViewer()->DrawGuides();

  gEve->Redraw3D(false);

  //__________________________________________
  // "tracks"

  std::set<std::pair<int32_t, int32_t>> drawnTracks;

  auto list = new TEveTrackList();
  auto prop = g_prop = list->GetPropagator();
  prop->SetFitDaughters(kFALSE);
  prop->SetStepper(TEveTrackPropagator::kRungeKutta);
  prop->SetMaxZ(500);
  prop->SetMaxOrbs(0.1);
  prop->SetMaxStep(0.1);
  prop->SetDelta(0.01);
  list->SetName("RK Propagator");
  list->SetLineColor(kMagenta);

  prop->SetMagFieldObj(new TEveMagFieldConst(-0.05, 0., 0.));
  list->SetElementName(Form("%s, constB", list->GetElementName()));

  int32_t count = 0;

  for (const auto& p : points) {
    int32_t event = p.event;
    int32_t track = p.track;
    if (drawnTracks.find(std::make_pair(event, track)) != drawnTracks.end())
      continue;
    auto rc = new TEveRecTrackD();
    rc->fV.Set(p.x, p.y, p.z);
    rc->fP.Set(p.px, p.py, p.pz);
    rc->fSign = -p.pdg; // electrons -> negative
    auto* evetrack = new TEveTrack(rc, prop);
    evetrack->SetName(fmt::format("event_{}_track_{}", event, track).c_str());
    evetrack->SetLineColor(list->GetLineColor());
    evetrack->SetLineWidth(2);
    drawnTracks.insert(std::make_pair(event, track));
    evetrack->MakeTrack();
    list->AddElement(evetrack);
    //
    ++count;
    if (count == 40)
      break;
  }

  gEve->AddElement(list);

  //__________________________________________

  // gEve->GetDefaultGLViewer()->SavePicture("fct_conversions.png");
}
