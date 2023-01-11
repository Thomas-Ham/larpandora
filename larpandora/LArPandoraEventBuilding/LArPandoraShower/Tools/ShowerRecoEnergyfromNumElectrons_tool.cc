//############################################################################
//### Name:        ShowerRecoEnergyfromNumElectrons                                ###
//### Author:      Tom Ham                                                 ###
//### Date:        01/04/2020                                              ###
//### Description: Tool for finding the Energy of the shower by going      ###
//###              from number of hits -> number of electrons -> energy.   ###
//###              Derived from the linear energy algorithm, written for   ###
//###              the EMShower_module.cc                                  ###
//############################################################################

//Framework Includes
#include "art/Utilities/ToolMacros.h"

//LArSoft Includes
#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Tools/IShowerTool.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"

//C++ Includes
#include <tuple>

#include "TFile.h"

namespace ShowerRecoTools {

  class ShowerRecoEnergyfromNumElectrons:IShowerTool {

    public:

      ShowerRecoEnergyfromNumElectrons(const fhicl::ParameterSet& pset);
      ~ShowerRecoEnergyfromNumElectrons();

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

     void FindLargestShower(std::vector<std::tuple<unsigned int, unsigned int, unsigned int, int, std::vector<double>, std::vector<int>, double>> recoenergy_largest_shower);
     double ModBoxRecombination(const detinfo::DetectorPropertiesData& detProp, double EField);

      art::InputTag fPFParticleLabel;
      int fVerbose;

      std::string fShowerEnergyOutputLabel;
      std::string fShowerBestPlaneOutputLabel;

      //Services
      art::ServiceHandle<geo::Geometry> fGeom;
      calo::CalorimetryAlg              fCalorimetryAlg;

      // Declare stuff
      double fRecombinationFactor;

      // Declare variables etc.
      double energy                = 0;
      double nominal_Efield        = 0;
      double localEfield           = 0;
      double localEfield_cweighted = 0;
      
      unsigned int showernum = 0;
      art::SubRunNumber_t subrun;
      art::EventNumber_t event;
      int hitsize;
      std::vector<double> hit_energy;
      int SP_hits;
      int SP_hits_temp;
      std::vector<int> hits_key;

      std::vector<std::tuple<unsigned int, unsigned int, unsigned int, int, std::vector<double>, std::vector<int>, double>> n_hit_energy; // more useful when making plots
      TTree *fOutputTree;
  };

  ShowerRecoEnergyfromNumElectrons::ShowerRecoEnergyfromNumElectrons(const fhicl::ParameterSet& pset) :
    IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools")),
    fPFParticleLabel(pset.get<art::InputTag>("PFParticleLabel")),
    fVerbose(pset.get<int>("Verbose")),
    fShowerEnergyOutputLabel(pset.get<std::string>("ShowerEnergyOutputLabel")),
    fShowerBestPlaneOutputLabel(pset.get<std::string>("ShowerBestPlaneOutputLabel")),
    fCalorimetryAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg")),
    fRecombinationFactor(pset.get<double>("RecombinationFactor"))
  {
    // save info 
    art::ServiceHandle<art::TFileService> tfs; 
    fOutputTree = tfs->make<TTree>("Reco","Reco");  
    fOutputTree->Branch("subrun", &subrun);    
    fOutputTree->Branch("event", &event);    
    fOutputTree->Branch("showernum", &showernum);    
    fOutputTree->Branch("numhits", &hitsize);    
    fOutputTree->Branch("hit_energy", &hit_energy);    
    fOutputTree->Branch("energy", &energy);    
    /*
    fOutputTree->Branch("sphits", &sphits);    
    fOutputTree->Branch("nosphits", &nosphits);    
    fOutputTree->Branch("hit_goodness_of_fit", &hit_goodness_of_fit);    
    fOutputTree->Branch("hit_peak_time", &hit_peak_time);    
    fOutputTree->Branch("hit_rms", &hit_rms);    
    fOutputTree->Branch("hit_peak_amplitude", &hit_peak_amplitude);    
    fOutputTree->Branch("hit_summed_adc", &hit_summed_adc);    
    fOutputTree->Branch("hit_integral", &hit_integral);    
    fOutputTree->Branch("hit_multiplicity", &hit_multiplicity);    
    fOutputTree->Branch("hit_charge", &hit_charge);    
    */

  }

  
  ShowerRecoEnergyfromNumElectrons::~ShowerRecoEnergyfromNumElectrons()
  {
    bool find_largest_shower = true;                                                                                         
    if(find_largest_shower){                                                                            
      FindLargestShower(n_hit_energy);
    }  
  }

  int ShowerRecoEnergyfromNumElectrons::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
      art::Event& Event,
      reco::shower::ShowerElementHolder& ShowerEleHolder
      ){

    if (fVerbose)
      std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Shower Reco Energy Tool ~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    
    showernum = ShowerEleHolder.GetShowerNumber();
    subrun = Event.subRun();
    event = Event.id().event();

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
    std::vector<double> energyVec(fGeom->Nplanes(), -999.);
    std::vector<double> energyError(fGeom->Nplanes(), -999.);

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(Event);
    auto const detProp   = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(Event, clockData);

    for(auto const& [plane, hits]: planeHits){

      unsigned int planeNumHits = hits.size();

      //Calculate the Energy for
      double energy = CalculateEnergy(clockData, detProp, hits, plane);
      std::cout << "energy :" << energy << std::endl;

      hitsize = hits.size();

      // If the energy is negative, leave it at -999
      if (energy>0)
        energyVec.at(plane) = energy;

      if (planeNumHits > bestPlaneNumHits) {
        bestPlane        = plane;
        bestPlaneNumHits = planeNumHits;
      }
      if(plane == 2){      
        fOutputTree->Fill();
        n_hit_energy.push_back(std::make_tuple(subrun, event, showernum, hitsize, hit_energy, hits_key, energy)); //save info for collection plane
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

  // function to calculate the reco energy
  double ShowerRecoEnergyfromNumElectrons::CalculateEnergy(const detinfo::DetectorClocksData& clockData,
          const detinfo::DetectorPropertiesData& detProp,
          const std::vector<art::Ptr<recob::Hit> >& hits,
          const geo::PlaneID::PlaneID_t plane){

    hit_energy.clear();
    double totalCharge = 0;
    double totalEnergy = 0;
    double correctedtotalCharge = 0;
    double nElectrons = 0;

    for (auto const& hit :hits){
      totalCharge = hit->Integral() * fCalorimetryAlg.LifetimeCorrection(clockData, detProp, hit->PeakTime()); // obtain charge and correct for lifetime
      // correct charge due to recombination
      correctedtotalCharge = totalCharge / fRecombinationFactor;
      // calculate # of electrons and the corresponding energy
      nElectrons = fCalorimetryAlg.ElectronsFromADCArea(correctedtotalCharge, plane);
      totalEnergy += (nElectrons / util::kGeVToElectrons) * 1000; // energy in MeV
      hit_energy.push_back((nElectrons / util::kGeVToElectrons) * 1000);
    }
    return totalEnergy;
  }




/*
// ModBox function to calculate recomobination
double ShowerRecoEnergyfromNumElectrons::ModBoxRecombination(const detinfo::DetectorPropertiesData& detProp, double EField){
    double rho   = detProp.Density();
    double Alpha = util::kModBoxA;
    double Beta  = util::kModBoxB/(rho * EField);
    double recombination = 0;
    recombination = std::log(Alpha + Beta * fNominalModBoxdEdx)/(Beta * fNominalModBoxdEdx);
//    std::cout << "Using the charge weighted local Efield of " << localEfield_cweighted << " kV/cm and a recombination factor of " << recombination << "." << std::endl;

    return recombination; 
}
*/

// Function to find only the largest shower
void ShowerRecoEnergyfromNumElectrons::FindLargestShower(std::vector<std::tuple<unsigned int, unsigned int, unsigned int, int, std::vector<double>, std::vector<int>, double>> recoenergy_largest_shower){
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
    TTree *recoenergy_LS = new TTree();
    art::ServiceHandle<art::TFileService> tfs;
    recoenergy_LS    = tfs->make<TTree>("Reco_LS","Reco_LS");

    recoenergy_LS->Branch("subrun", &subrun, "subrun/i");                                                                              
    recoenergy_LS->Branch("event", &event, "event/i");    
    recoenergy_LS->Branch("showernum", &showernum, "showernum/i");
    recoenergy_LS->Branch("numhits", &hitsize, "numhits/I");                                                     
    recoenergy_LS->Branch("hit_energy", &hit_energy);                                                     
    recoenergy_LS->Branch("hits_key", &hits_key);                                                     
   // recoenergy_LS->Branch("SPHits_temp", &SP_hits_temp, "SPHits_temp/I");                                                     
    recoenergy_LS->Branch("energy", &energy, "energy/d"); 

    for(unsigned int i = 0; i < recoenergy_largest_shower.size(); i++){
        subrun   = std::get<0>(recoenergy_largest_shower[i]);
        event    = std::get<1>(recoenergy_largest_shower[i]);
        showernum = std::get<2>(recoenergy_largest_shower[i]);
        hitsize   = std::get<3>(recoenergy_largest_shower[i]);
        hit_energy= std::get<4>(recoenergy_largest_shower[i]);
        hits_key  = std::get<5>(recoenergy_largest_shower[i]);
        energy    = std::get<6>(recoenergy_largest_shower[i]);
        //SP_hits_temp   = std::get<5>(recoenergy_largest_shower[i]);

        recoenergy_LS->Fill();
    }

//    m_file->Write("", TObject::kOverwrite); // save only the new version of the tree - without the arguments was duplicating original tree
}


}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerRecoEnergyfromNumElectrons)
