//############################################################################
//### Name:        ShowerLinearEnergy                                      ###
//### Author:      Dominic Barker                                          ###
//### Date:        13.05.19                                                ###
//### Description: Tool for finding the Energy of the shower. Derived      ###
//###              from the linear energy algorithm, written for           ###
//###              the EMShower_module.cc                                  ###
//############################################################################

//Framework Includes
#include "art/Utilities/ToolMacros.h"

//LArSoft Includes
#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Tools/IShowerTool.h"
#include "lardataobj/RecoBase/Cluster.h"

// ROOT includes
#include "TFile.h"
#include "TTree.h" 
#include "TH1D.h" 

//C++ Includes
#include <iostream>
#include <vector>
#include <tuple>

namespace ShowerRecoTools {

  class ShowerLinearEnergy:IShowerTool {

    public:

      ShowerLinearEnergy(const fhicl::ParameterSet& pset);

      ~ShowerLinearEnergy(); 

      //Physics Function. Calculate the shower Energy.
      int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
          art::Event& Event,
          reco::shower::ShowerElementHolder& ShowerElementHolder
          ) override;
    private:

      double CalculateEnergy(const detinfo::DetectorClocksData& clockData,
          const detinfo::DetectorPropertiesData& detProp,
          const std::vector<art::Ptr<recob::Hit> >& hits,
          const geo::PlaneID::PlaneID_t plane) const;

      void FindLargestShower(std::vector<std::tuple<unsigned int, unsigned int, unsigned int, int, double>> recoenergy_largest_shower);
      //fcl parameters
      unsigned int        fNumPlanes;
      std::vector<double> fGradients;   //Gradient of the linear fit of total charge to total energy
      std::vector<double> fIntercepts;  //Intercept of the linear fit of total charge to total energy

      art::InputTag fPFParticleLabel;
      int           fVerbose;

      std::string fShowerEnergyOutputLabel;
      std::string fShowerBestPlaneOutputLabel;

      //Services
      art::ServiceHandle<geo::Geometry> fGeom;

      double Energy = 0;
      TTree *fOutputTree;
      std::string Tree_name = "Reco";
      unsigned int showernum = 0;
      art::SubRunNumber_t subRunN;
      art::EventNumber_t EventN;
      int hitsize;
      // vec to store subrun #, event #, shower #, # of hits and energy
      std::vector<std::tuple<unsigned int, unsigned int, unsigned int, int,  double>> n_hit_energy; // more useful when making plots

  };

  ShowerLinearEnergy::ShowerLinearEnergy(const fhicl::ParameterSet& pset) :
    IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools")),
    fGradients(pset.get<std::vector<double> >("Gradients")),
    fIntercepts(pset.get<std::vector<double> >("Intercepts")),
    fPFParticleLabel(pset.get<art::InputTag>("PFParticleLabel")),
    fVerbose(pset.get<int>("Verbose")),
    fShowerEnergyOutputLabel(pset.get<std::string>("ShowerEnergyOutputLabel")),
    fShowerBestPlaneOutputLabel(pset.get<std::string>("ShowerBestPlaneOutputLabel"))
  {
    fNumPlanes = fGeom->Nplanes();
    if (fNumPlanes!=fGradients.size() || fNumPlanes!=fIntercepts.size()){
      throw cet::exception("ShowerLinearEnergy")
        << "The number of planes does not match the size of the fcl parametes passed: Num Planes: "
        << fNumPlanes << ", Gradients size: " << fGradients.size() << ", Intercpts size: "
        << fIntercepts.size();
    }

    art::ServiceHandle<art::TFileService> tfs;
    fOutputTree = tfs->make<TTree>(Tree_name.c_str(),Tree_name.c_str());
           
    //add branches                                                                          
    fOutputTree->Branch("subrun", &subRunN, "subrun/i");
    fOutputTree->Branch("event", &EventN, "event/i");
    fOutputTree->Branch("showernum", &showernum, "showernum/i");
    fOutputTree->Branch("numhits", &hitsize, "numhits/I");
    fOutputTree->Branch("energy", &Energy, "energy/d");

  }

  ShowerLinearEnergy::~ShowerLinearEnergy()
  {
    FindLargestShower(n_hit_energy);
  }

  int ShowerLinearEnergy::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
      art::Event& Event, reco::shower::ShowerElementHolder& ShowerEleHolder
      ){

    // get shower number per event
    showernum = ShowerEleHolder.GetShowerNumber();
    // get subrun number
    subRunN = Event.subRun();
    // get event number 
    EventN = Event.id().event();

    // Get the assocated pfParicle vertex PFParticles
    auto const pfpHandle = Event.getValidHandle<std::vector<recob::PFParticle> >(fPFParticleLabel);

    //Get the clusters
    auto const clusHandle = Event.getValidHandle<std::vector<recob::Cluster> >(fPFParticleLabel);

    const art::FindManyP<recob::Cluster>& fmc = ShowerEleHolder.GetFindManyP<recob::Cluster>(
        pfpHandle, Event, fPFParticleLabel);
    // art::FindManyP<recob::Cluster> fmc(pfpHandle, Event, fPFParticleLabel);
    std::vector<art::Ptr<recob::Cluster> > clusters = fmc.at(pfparticle.key());

    //Get the hit association
    const art::FindManyP<recob::Hit>& fmhc = ShowerEleHolder.GetFindManyP<recob::Hit>(
        clusHandle, Event, fPFParticleLabel);
    // art::FindManyP<recob::Hit> fmhc(clusHandle, Event, fPFParticleLabel);

    std::map<geo::PlaneID::PlaneID_t, std::vector<art::Ptr<recob::Hit> > > planeHits;

    //Loop over the clusters in the plane and get the hits
    for(auto const& cluster: clusters){

      //Get the hits
      std::vector<art::Ptr<recob::Hit> > hits = fmhc.at(cluster.key());

      //Get the plane.
      const geo::PlaneID::PlaneID_t plane(cluster->Plane().Plane);

      planeHits[plane].insert(planeHits[plane].end(),hits.begin(),hits.end());
    }

    // Calculate the energy for each plane && best plane
    geo::PlaneID::PlaneID_t bestPlane  = std::numeric_limits<geo::PlaneID::PlaneID_t>::max();
    unsigned int bestPlaneNumHits = 0;

    //Holder for the final product
    std::vector<double> energyVec(fNumPlanes, -999.);
    std::vector<double> energyError(fNumPlanes, -999.);

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(Event);
    auto const detProp   = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(Event, clockData);

    for(auto const& [plane, hits]: planeHits){

      unsigned int planeNumHits = hits.size();

      //Calculate the Energy for
      Energy = CalculateEnergy(clockData, detProp, hits,plane);
      // If the energy is negative, leave it at -999
      if (Energy>0)
        energyVec.at(plane) = Energy;

      if (planeNumHits > bestPlaneNumHits) {
        bestPlane        = plane;
        bestPlaneNumHits = planeNumHits;
      }
      if(fVerbose){
        std::cout << "plane: " << plane
                  << "  hits: " << hits.size()
                  << "  energy: " << Energy << std::endl;  
      }

      if(plane == 2){   
           
        fOutputTree->Fill();
        n_hit_energy.push_back(std::make_tuple(subRunN, EventN, showernum, hits.size(), Energy)); //save info for collection plane
      }

    }

    ShowerEleHolder.SetElement(energyVec, energyError, fShowerEnergyOutputLabel);
    // Only set the best plane if it has some hits in it
    if (bestPlane < fGeom->Nplanes()){
      // Need to cast as an int for legacy default of -999
      // have to define a new variable as we pass-by-reference when filling
      int bestPlaneVal(bestPlane);
      ShowerEleHolder.SetElement(bestPlaneVal, fShowerBestPlaneOutputLabel);
    }

    return 0;
  }

  //Function to calculate the energy of a shower in a plane. Using a linear map between charge and Energy.
  //Exactly the same method as the ShowerEnergyAlg.cxx. Thanks Mike.
  double ShowerLinearEnergy::CalculateEnergy(const detinfo::DetectorClocksData& clockData,
          const detinfo::DetectorPropertiesData& detProp,
          const std::vector<art::Ptr<recob::Hit> >& hits,
          const geo::PlaneID::PlaneID_t plane) const {

    double totalCharge = 0, totalEnergy = 0;

    for (auto const& hit: hits){
      totalCharge += (hit->Integral() * std::exp((sampling_rate(clockData)* hit->PeakTime()) / (detProp.ElectronLifetime()*1e3) ) );
    }

    totalEnergy = (totalCharge * fGradients.at(plane)) + fIntercepts.at(plane);

    return totalEnergy;

  }








// Function to find only the largest shower
void ShowerLinearEnergy::FindLargestShower(std::vector<std::tuple<unsigned int, unsigned int, unsigned int, int, double>> recoenergy_largest_shower){
    // Cut showers so we are only left with the biggest (most hits) from each event.
    // Feel like there should be a more ROOT-y way to do this...
    unsigned int i = 0;                                                                 
    while(i < recoenergy_largest_shower.size()){  

        // If the shower at (i-1)'th position from the same event has fewer hits, delete it
        if(std::get<2>(recoenergy_largest_shower[i]) != 0 && std::get<3>(recoenergy_largest_shower[i]) > std::get<3>(recoenergy_largest_shower[(i-1)])){   

            recoenergy_largest_shower.erase(recoenergy_largest_shower.begin() + (i-1));   
        }
        // Delete any remaining i'th non primary shower (i.e. the non primary shower which has fewer hits than the one at (i-1))
        else if(std::get<2>(recoenergy_largest_shower[i]) != 0){ 
    
            recoenergy_largest_shower.erase(recoenergy_largest_shower.begin() + (i));
        }

        else{
            i++;
        }
    }


    // Add new tree to root file which has only the largest shower from each event
    TTree *recoenergy_LS = new TTree((Tree_name + "LS").c_str(), (Tree_name + "LS").c_str());
    art::ServiceHandle<art::TFileService> tfs;
    recoenergy_LS    = tfs->make<TTree>("recoenergy_LS","Reco Energy Largest Shower");

    recoenergy_LS->Branch("subrun", &subRunN, "subrun/i");                                                                              
    recoenergy_LS->Branch("event", &EventN, "event/i");    
    recoenergy_LS->Branch("showernum", &showernum, "showernum/i");
    recoenergy_LS->Branch("numhits", &hitsize, "numhits/I");                                                     
    recoenergy_LS->Branch("energy", &Energy, "energy/d"); 

    for(unsigned int i = 0; i < recoenergy_largest_shower.size(); i++){
        subRunN   = std::get<0>(recoenergy_largest_shower[i]);
        EventN    = std::get<1>(recoenergy_largest_shower[i]);
        showernum = std::get<2>(recoenergy_largest_shower[i]);
        hitsize   = std::get<3>(recoenergy_largest_shower[i]);
        Energy    = std::get<4>(recoenergy_largest_shower[i]);

        recoenergy_LS->Fill();
    }

}





}


DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerLinearEnergy)


