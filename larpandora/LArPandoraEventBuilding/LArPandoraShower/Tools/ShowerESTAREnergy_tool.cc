//############################################################################
//### Name:        ShowerESTAREnergy                                       ###
//### Author:      Tom Ham			                                           ###
//### Date:        06/09/2021                                              ###
//### Description: Tool for finding the energy of a shower by using        ###
//###              the ESTAR database.                                     ###
//###              Derived from the method used by ArgoNeuT.               ###
//############################################################################

//Framework Includes
#include "art/Utilities/ToolMacros.h"

//LArSoft Includes
#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Tools/IShowerTool.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"

// ROOT includes
#include "TFile.h"
#include "TGraph2D.h"

namespace ShowerRecoTools {

  class ShowerESTAREnergy:IShowerTool {

  public:

  ShowerESTAREnergy(const fhicl::ParameterSet& pset);
    
  ~ShowerESTAREnergy(); 
        
  //Physics Function. Calculate the shower Energy.
  int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
  art::Event& Event,
  reco::shower::ShowerElementHolder& ShowerElementHolder
  ) override;

  private:
    
  double CalculateEnergy(const detinfo::DetectorClocksData& clockData,
           const detinfo::DetectorPropertiesData& detProp, 
           const std::vector<art::Ptr<recob::Hit> >& hits,
           const geo::PlaneID::PlaneID_t plane);

  art::InputTag fPFParticleLabel;
  int fVerbose;

  //Services
  art::ServiceHandle<geo::Geometry> fGeom;
  calo::CalorimetryAlg              fCalorimetryAlg;    

  //fcl params 
  std::string fShowerEnergyOutputLabel;
  std::string fShowerBestPlaneOutputLabel; 
  bool fSCECorrectEField;
  std::string fESTARLookupCurve;
  std::string fESTARTGraph;

  // Read in the ESTAR lookup curve
  // Need to change this once the lookup curve has been committed to sbnd_data
  TFile *estar = new TFile(fESTARLookupCurve.c_str());
  TGraph2D *t_estar = (TGraph2D*)estar->Get(fESTARTGraph.c_str());

  double nominal_Efield = 0;         // Nominal E-field in the detector
  double localEfield = 0;            // Local E-field at the geometric centre of the shower 
  double localEfield_cweighted = 0;  // Local E-field at the charge weighted centre of the shower

  std::vector<std::tuple<int, double>> sp_hits_tuple; // tuple holding the key and local E-field for all the sp hits 

  }; // class

  ShowerESTAREnergy::ShowerESTAREnergy(const fhicl::ParameterSet& pset) :
    IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools")),
    fPFParticleLabel(pset.get<art::InputTag>("PFParticleLabel")),
    fVerbose(pset.get<int>("Verbose")),
    fCalorimetryAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg")),
    fShowerEnergyOutputLabel(pset.get<std::string>("ShowerEnergyOutputLabel")),
    fShowerBestPlaneOutputLabel(pset.get<std::string>("ShowerBestPlaneOutputLabel")),
    fSCECorrectEField(pset.get<bool>("SCECorrectEField")),
    fESTARLookupCurve(pset.get<std::string>("ESTARLookupCurve")),
    fESTARTGraph(pset.get<std::string>("ESTARTGraph"))
  {
  }

  ShowerESTAREnergy::~ShowerESTAREnergy()
  {
  }

  int ShowerESTAREnergy::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
        art::Event& Event,
        reco::shower::ShowerElementHolder& ShowerEleHolder
  ){

    if(fVerbose){
      std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Shower ESTAR Energy Tool ~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
      if(fSCECorrectEField){std::cout << "Applying SCE correction" << std::endl;}
    }
      
    // Get the associated pfParicle vertex PFParticles
    auto const pfpHandle = Event.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleLabel);

    // Get the clusters
    auto const clusHandle = Event.getValidHandle<std::vector<recob::Cluster>>(fPFParticleLabel);

    const art::FindManyP<recob::Cluster>& fmc = ShowerEleHolder.GetFindManyP<recob::Cluster>(pfpHandle, Event, fPFParticleLabel);
    std::vector<art::Ptr<recob::Cluster> > clusters = fmc.at(pfparticle.key());

    // Get the hits association from the clusters
    const art::FindManyP<recob::Hit>& fmhc = ShowerEleHolder.GetFindManyP<recob::Hit>(clusHandle, Event, fPFParticleLabel);
    
    std::map<geo::PlaneID::PlaneID_t, std::vector<art::Ptr<recob::Hit>>> planeHits;

    // Get the spacepoints
    auto const spHandle = Event.getValidHandle<std::vector<recob::SpacePoint>>(fPFParticleLabel);

    const art::FindManyP<recob::SpacePoint>& fmsp = ShowerEleHolder.GetFindManyP<recob::SpacePoint>(pfpHandle, Event, fPFParticleLabel);
    std::vector<art::Ptr<recob::SpacePoint> > spacepoints = fmsp.at(pfparticle.key());

    // Get the hits association from the spacepoints
    const art::FindManyP<recob::Hit>& fmhsp = ShowerEleHolder.GetFindManyP<recob::Hit>(spHandle, Event, fPFParticleLabel);

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(Event);
    auto const detProp   = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(Event, clockData);
    
    // Get E-field -  Electric Field in the drift region in KV/cm
    nominal_Efield  = detProp.Efield(); // Nominal E-field 
   
    // Not all hits from the clusters will match to Space points. If that's the case let's use some sort of shower centre
    // Get shower centre 
    TVector3 showercentre = IShowerTool::GetLArPandoraShowerAlg().ShowerCentre(spacepoints);
    
    // Get the charge weighted shower centre - sum(charge*position)/sum(charge)
    // Charge is obtained from the lifetime corrected hits and hits are obtained from the spacepoints 
    TVector3 chargecentre = IShowerTool::GetLArPandoraShowerAlg().ShowerCentre(clockData, detProp, spacepoints, fmhsp);
    
    // Get space charge corrected E-field
    if(fSCECorrectEField){
      // Check the shower centre exists and is sensible
      // Think the charge weighted centre can take nan values which causes the E-field calculation to fail
      if(showercentre.Mag() >= 0 && chargecentre.Mag() >=0){ 
        localEfield = IShowerTool::GetLArPandoraShowerAlg().SCECorrectEField(nominal_Efield, showercentre);
        localEfield_cweighted = IShowerTool::GetLArPandoraShowerAlg().SCECorrectEField(nominal_Efield, chargecentre);
      }
      else{
        mf::LogWarning("ShowerESTAREnergy") << "Shower centre calculation doesn't look to be sensible. Reconstruction is probably dodgy." << std::endl;
      }
    }

    // Loop over the clusters in the plane and get the hits
    for(auto const& cluster: clusters){
          
      // Get the hits
      std::vector<art::Ptr<recob::Hit> > hits = fmhc.at(cluster.key());
      // Get the plane.
      const geo::PlaneID::PlaneID_t plane(cluster->Plane().Plane);
       
      planeHits[plane].insert(planeHits[plane].end(),hits.begin(),hits.end());
    }

    // Loop over the spacepoints and get the hits
    sp_hits_tuple.clear();
    for(auto const& spacepoint : spacepoints){
      //Get the hits
      std::vector<art::Ptr<recob::Hit> > hits = fmhsp.at(spacepoint.key());

      // get local Efield at SP_pos
      TVector3 sp_pos = IShowerTool::GetLArPandoraShowerAlg().SpacePointPosition(spacepoint); 
      double localEfield_SP = IShowerTool::GetLArPandoraShowerAlg().SCECorrectEField(nominal_Efield, sp_pos);

      // Store the localEfield of all the hits with associated SP's so we can compare with the hits from the clusters later
      for(auto const& h : hits){
        sp_hits_tuple.push_back(std::make_tuple(h.key(), localEfield_SP));
      }           
    }

    // Calculate the energy for each plane && best plane
    geo::PlaneID::PlaneID_t bestPlane  = std::numeric_limits<geo::PlaneID::PlaneID_t>::max();
    unsigned int bestPlaneNumHits = 0;

    //Holder for the final product
    std::vector<double> ShowerESTAREnergy(fGeom->Nplanes(), -999);
    std::vector<double> ShowerESTAREnergyError(fGeom->Nplanes(), -999);

    for(auto const& [plane, hits] : planeHits){

      unsigned int planeNumHits = hits.size();
        
      //Calculate the Energy
      double Energy = CalculateEnergy(clockData, detProp, hits, plane);

      if(Energy > 0){    
        ShowerESTAREnergy.at(plane) = Energy;
      }

      if(planeNumHits > bestPlaneNumHits){
        bestPlane        = plane;
        bestPlaneNumHits = planeNumHits;
      }
      if(fVerbose){
        std::cout << "plane: " << plane
                  << "  hits: " << hits.size()
                  << "  energy: " << Energy << std::endl;  
      }

    }

    ShowerEleHolder.SetElement(ShowerESTAREnergy,ShowerESTAREnergyError,fShowerEnergyOutputLabel);

    // Only set the best plane if it has some hits in it
    if (bestPlane < fGeom->Nplanes()){
      // Need to cast as an int for legacy default of -999
      // have to define a new variable as we pass-by-reference when filling
      int bestPlaneVal(bestPlane);
      ShowerEleHolder.SetElement(bestPlaneVal, fShowerBestPlaneOutputLabel);
      if(fVerbose){
        std::cout << "best plane: " << bestPlane << ",  best plane number of hit: " << bestPlaneNumHits << std::endl;
      }
    }
     
    return 0;
  }


  // function to calculate the reco energy
  double ShowerESTAREnergy::CalculateEnergy(const detinfo::DetectorClocksData& clockData,
         const detinfo::DetectorPropertiesData& detProp,
         const std::vector<art::Ptr<recob::Hit> >& hits, const geo::PlaneID::PlaneID_t plane){
   
    double totalEnergy = 0;
    double efield = 0; 
     
    // Loop over the hits
    for(auto const& h : hits){
      // If we're correcting for SCE, check if the hit is one that was we found to have an associated SP
      if(fSCECorrectEField){
        auto it = std::find_if(sp_hits_tuple.begin(), sp_hits_tuple.end(), [&](const std::tuple<unsigned int, double>& e){
        return std::get<0>(e) == h.key();
      });
      
      // do a per hit SCE correction since we have associated space points
      if(it != sp_hits_tuple.end()){
        efield = std::get<1>(*it);
      }

      // no associated space points :( so correct these hits using the shower centre position 
      else{
        efield = localEfield_cweighted;
      }
    }

    else{
      efield = nominal_Efield; 
    }

    double totalCharge = (h)->Integral() * fCalorimetryAlg.LifetimeCorrection(clockData, detProp, (h)->PeakTime()); // obtain charge and correct for lifetime
    double nElectrons = fCalorimetryAlg.ElectronsFromADCArea(totalCharge, plane);
    if(nElectrons !=0 && t_estar->Interpolate(nElectrons, efield) <= 0){
      mf::LogWarning("ShowerESTAREnergy") << "Interpolation is returning an energy of 0 even though we do have some collected charge. " 
                                          << "Probably due to an out of bounds entry for the charge or the efield." << std::endl; 
    }
    totalEnergy += t_estar->Interpolate(nElectrons, efield);

    }
    return totalEnergy; 
  }
} // namespace

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerESTAREnergy)










