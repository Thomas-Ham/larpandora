//############################################################################
//### Name:        ShowerEnergyCheater                                     ###
//### Author:      Ed Tyley                                                ###
//### Date:        16.07.19                                                ###
//### Description: Cheating tool using truth for shower energy             ###
//############################################################################


//Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindManyP.h"

//LArSoft Includes
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Shower.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Tools/IShowerTool.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Algs/LArPandoraShowerCheatingAlg.h"
#include "sbndcode/RecoUtils/RecoUtils.h"                                                                                                                                              

//ROOT Inclides
#include "TFile.h"
#include "TTree.h"

//C++ Includes
#include <iostream>
#include <cmath>
#include <fstream>

//Root Includes

namespace ShowerRecoTools {

  class ShowerEnergyCheater:IShowerTool {

    public:

      ShowerEnergyCheater(const fhicl::ParameterSet& pset);

      ~ShowerEnergyCheater();

      //Calculate Cheating StartPosition 
      int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
			   art::Event& Event, 
			   reco::shower::ShowerElementHolder& ShowerEleHolder) override;

    private:

        void FindLargestShower(std::vector<std::tuple<unsigned int, unsigned int, unsigned int, int, std::vector<double>, double, int, int, int, int>> energy_largest_shower);     
        int TrueParticleID(detinfo::DetectorClocksData const& clockData, const art::Ptr<recob::Hit> hit, bool rollup_unsaved_ids);
        std::map<int,std::map<geo::PlaneID,int> > NumberofPlaneHitsPerTrack(detinfo::DetectorClocksData const& clockData, const std::vector<art::Ptr<recob::Hit> >& hits);

        //Algorithm functions
        shower::LArPandoraShowerCheatingAlg fLArPandoraShowerCheatingAlg;

        //FCL
        art::InputTag fPFParticleModuleLabel;
        art::InputTag fHitModuleLabel;
        art::InputTag fShowerModuleLabel;
        art::ServiceHandle<cheat::BackTrackerService> bt_serv;
 
        unsigned int showernum = 0;
        art::SubRunNumber_t subRunN;                                                                                             
        art::EventNumber_t EventN; 
        double trueEnergy;
        double purity;
        double completeness;
        int hitsize;
        int TrueHitsDep_WithinRecoShower;
        int TrueHitsDep_FromTrueShower;
        int total_event_pfps;
        int this_pfp;
        std::vector<double> true_hit_energy;

	std::vector<std::tuple<unsigned int, unsigned int, unsigned int, int, std::vector<double>, double, int, int, int, int>> n_hit_energy;

        std::string File_name = "Energy_files/Cathode_sample_SCE.root"; // The file to write to
        std::string Tree_name = "energy_cheater"; 
        TFile *m_file;    ///< TFile for saving event information                                                               
        TTree *fOutputTree;  
 };


    ShowerEnergyCheater::ShowerEnergyCheater(const fhicl::ParameterSet& pset) :
        IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools")),
        fLArPandoraShowerCheatingAlg(pset.get<fhicl::ParameterSet>("LArPandoraShowerCheatingAlg")),
        fPFParticleModuleLabel(pset.get<art::InputTag>("PFParticleModuleLabel","")),
        fHitModuleLabel(pset.get<art::InputTag>("HitModuleLabel")),
        fShowerModuleLabel(pset.get<art::InputTag>("ShowerModuleLabel"))
    {
        //fShowerModuleLabel = pset.get<std::vector<std::string> >("ShowerModuleLabel");                                                                                            
        
        // The TFileService lets us define a tree and takes care of writing it to file          
        art::ServiceHandle<art::TFileService> tfs;                                                                                                                       
     //   m_file      = new TFile(File_name.c_str(), "UPDATE");                          
        fOutputTree = tfs->make<TTree>(Tree_name.c_str(),Tree_name.c_str());                                                          

        //add branches                                                                      
        fOutputTree->Branch("Subrun", &subRunN, "Subrun/i");                                          
        fOutputTree->Branch("Event", &EventN, "Event/i");                                      
        fOutputTree->Branch("ShowerN", &showernum, "ShowerN/i");                                        
        fOutputTree->Branch("NHits", &hitsize, "NHits/I");                      
        fOutputTree->Branch("Hit_Energy", &true_hit_energy);                      
        fOutputTree->Branch("Energy", &trueEnergy, "Energy/d");   
        fOutputTree->Branch("TrueHitsDep_FromTrueShower", &TrueHitsDep_FromTrueShower, "TrueHitsDep_FromTrueShower/I");   
        fOutputTree->Branch("TrueHitsDep_WithinRecoShower", &TrueHitsDep_WithinRecoShower, "TrueHitsDep_WithinRecoShower/I");   
        fOutputTree->Branch("Total_Event_PFPs", &total_event_pfps, "Total_Event_PFPs/I");
        fOutputTree->Branch("This_PFP", &this_pfp, "This_PFP/I");
    }

    ShowerEnergyCheater::~ShowerEnergyCheater()
    {
       // m_file->cd();                                                                                   
       // fOutputTree->CloneTree()->Write("energy_cheater", TObject::kOverwrite); 

       bool find_largest_shower = true;                                                                                         
        if(find_largest_shower){                                                                            
        FindLargestShower(n_hit_energy);                                                                                                    
       }     

     //   m_file->Close();
    }

    int ShowerEnergyCheater::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle, 
						   art::Event& Event, 
						   reco::shower::ShowerElementHolder& ShowerEleHolder){

    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Shower Energy Cheater ~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

    showernum = ShowerEleHolder.GetShowerNumber();
		
    // get subrun number
    subRunN = Event.subRun();
    // get event number 
    EventN = Event.id().event();

    art::ServiceHandle<geo::Geometry> geom;

    //Could store these in the shower element holder and just calculate once?
    std::map<int,const simb::MCParticle*> trueParticles = fLArPandoraShowerCheatingAlg.GetTrueParticleMap();
    std::map<int,std::vector<int> > showersMothers = fLArPandoraShowerCheatingAlg.GetTrueChain(trueParticles);

    //Get the hits from the shower:
    art::Handle<std::vector<recob::PFParticle> > pfpHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, pfpHandle)){
      throw cet::exception("ShowerEnergyCheater") << "Could not get the pandora pf particles. Something is not cofingured coreectly Please give the correct pandoa module label. Stopping";
      return 1;
    }
    
    std::vector<art::Ptr<recob::PFParticle>> pfps;   
    art::fill_ptr_vector(pfps, pfpHandle); 
    total_event_pfps = pfps.size();
    this_pfp = pfparticle->Self();


    //Get the clusters
    art::Handle<std::vector<recob::Cluster> > clusHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, clusHandle)){
      throw cet::exception("ShowerEnergyCheater") << "Could not get the pandora clusters. Something is not cofingured coreectly Please give the correct pandoa module label. Stopping";
      return 1;
    }

    //Getting the Shower Information                                                                                                                                                   
    art::Handle<std::vector<recob::Shower> > showerListHandle;
   
    std::vector<art::Ptr<recob::Shower> > recoShowers;                                                                                                                                                                              
    if(Event.getByLabel(fShowerModuleLabel,showerListHandle)){                                                                                                                                                                                  
        art::fill_ptr_vector(recoShowers,showerListHandle);                                                                                                                                     
    }              

    //Association between Showers and 2d Hits                                                                                                                                          
    art::FindManyP<recob::Hit> fmh(showerListHandle, Event, fShowerModuleLabel); 

    art::ProductID showerhit_productid = fmh.at(0).front().id();
    std::cout << "showerhit_productid: " << showerhit_productid << std::endl; 
   
    art::FindManyP<recob::Cluster> fmc(pfpHandle, Event, fPFParticleModuleLabel);
    std::vector<art::Ptr<recob::Cluster> > clusters = fmc.at(pfparticle.key());

    //Get the hit association
    art::FindManyP<recob::Hit> fmhc(clusHandle, Event, fPFParticleModuleLabel);

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(Event);

    std::map<art::ProductID,std::map<int,std::map<geo::PlaneID, int> > > MCTrack_hit_map; 
    //Get all the hits                                                                                                                                                                   
    std::vector<art::Handle<std::vector<recob::Hit> > > hitHandles;                                                                                                                      
    Event.getManyByType(hitHandles);   
    //Get the number of hits associated wit each Track. This is done for every hits handle in the event.                                                                               
    for(auto& handle : hitHandles){                                                                                                                                                    
                                                                                                                                                                                           
      if(!handle.isValid()){                                                                                                                                                           
        mf::LogError("ShowerValidation") << "Bad hit handle from the all the hit handles" << std::endl;                                                                                
        continue;                                                                                                                                                                      
      }                                                                                                                                                                                
                                                                                                                                                                                           
      //RawHitFinder is a bit silly at the moment so ignore it for the time being - REMOVE                                                                                             
      if(handle.provenance()->moduleLabel() == "fasthit"){continue;}                                                                                                                   
                                                                                                                                                                                        
      //Getting the Hit Information                                                                                                                                                    
      art::Handle<std::vector<recob::Hit> > hitHandle;                                                                                                                                 
      std::vector<art::Ptr<recob::Hit> > hits_fromhandle;                                                                                                                              
      if(Event.getByLabel(handle.provenance()->moduleLabel(),hitHandle))                                                                                                                 
          {art::fill_ptr_vector(hits_fromhandle, hitHandle);}                                                                                                                              
                                                                                                                                                                                                            
      //Get a map of the number of hits per plane each track has deposited.                                                                                                            
      MCTrack_hit_map[handle.id()] = ShowerEnergyCheater::NumberofPlaneHitsPerTrack(clockData, hits_fromhandle);                                                                                 
    }                                                 

    std::pair<int,double> ShowerTrackInfo;
    const simb::MCParticle* trueParticle = nullptr;

    std::map<geo::PlaneID::PlaneID_t, std::vector<art::Ptr<recob::Hit> > > planeHits;

    for(auto const& cluster: clusters){
        
        // Get the plane
        const geo::PlaneID::PlaneID_t plane(cluster->Plane().Plane);
        //Get the hits
        std::vector<art::Ptr<recob::Hit>> hits = fmhc.at(cluster.key());

        planeHits[plane].insert(planeHits[plane].end(),hits.begin(),hits.end());
    }

    // Get energy etc for each plane
    for(auto const& [plane, hits]: planeHits){
        true_hit_energy.clear();
                
        ShowerTrackInfo = fLArPandoraShowerCheatingAlg.TrueParticleIDFromTrueChain(clockData,showersMothers,hits,plane);
        trueParticle = trueParticles[ShowerTrackInfo.first];
        int TrueHitsDep_WithinRecoShower = 0;
        int TrueHitsDep_FromTrueShower = 0;
        int ShowerTrackID = ShowerTrackInfo.first;

        for (auto const& daughter: showersMothers[ShowerTrackID]){                                                                                                                      
            for(geo::PlaneID plane_id: geom->IteratePlaneIDs()){                                                                                                                        
                //std::cout << "plane_id = " << plane_id << std::endl;
                if(plane_id.Plane == plane){
                    TrueHitsDep_FromTrueShower +=  MCTrack_hit_map[showerhit_productid][daughter][plane_id];                                                                                    
                }
            }                                                                                                                                                                            
        }   
            
        //Get the number of hits in the reco shower from the true shower.                                                                                                              
        std::map<int,std::map<geo::PlaneID,int> >  MCTrack_showerhit_map = ShowerEnergyCheater::NumberofPlaneHitsPerTrack(clockData, hits);                
        for(auto const& planehit_map : MCTrack_showerhit_map){                                                                                                                         
            if(std::find(showersMothers[ShowerTrackID].begin(), showersMothers[ShowerTrackID].end(),planehit_map.first) == showersMothers[ShowerTrackID].end()){continue;}                  
            for(auto const& hit_num : planehit_map.second){                                                                                                                            
                if(hit_num.first.Plane == plane){
                    //std::cout << "planehit_map.second: " << hit_num.first.Plane << std::endl;
                    TrueHitsDep_WithinRecoShower += hit_num.second;                                                                                                                            
                }
            }                                                                                                                                                                            
        } 

        purity = (double)TrueHitsDep_WithinRecoShower / (double) hits.size();
        completeness = (double)TrueHitsDep_WithinRecoShower / (double) TrueHitsDep_FromTrueShower;

        // Get the true energy of showering particle
        if(ShowerTrackInfo.first==-99999) {
            mf::LogError("ShowerEnergyCheater") << "True Shower Not Found. Setting hits and energy to -999.";
            hitsize    = -999;
            trueEnergy = -999;
            
        }
        else{
            trueEnergy = 0;
            hitsize = hits.size();
            // Get true energy of showering particle 
            trueEnergy = trueParticle->E() * 1000; // change to MeV
        }
    
        std::cout << "Plane: " << plane << std::endl;   
        std::cout << "Hits: " << hitsize << std::endl;
        std::cout << "TrueHitsDep_WithinRecoShower: " << TrueHitsDep_WithinRecoShower << std::endl;       
        std::cout << "TrueHitsDep_FromTrueShower: " << TrueHitsDep_FromTrueShower << std::endl;
        std::cout << "purity: " << purity << std::endl;
        std::cout << "completeness: " << completeness << std::endl;
        std::cout << "True Energy: " << trueEnergy << std::endl;

        // Get true hit energy
        double deposit = 0;
        for(auto const& h : hits){
            std::vector<sim::TrackIDE> trackIDEs = bt_serv->HitToTrackIDEs(clockData, h);

            std::vector<sim::TrackIDE> ides;
            for( auto ide : trackIDEs ){
                ides.push_back(ide);
            }
            double hit_energy = 0;
            for (std::vector<sim::TrackIDE>::iterator ideIt = ides.begin(); ideIt != ides.end(); ++ideIt){
                hit_energy += ideIt->energy;
                deposit += ideIt->energy;
            }
            true_hit_energy.push_back(hit_energy);
            std::cout << "hit energy: " << hit_energy << std::endl;
        }
        std::cout << "sum of hits energy: " << deposit << std::endl;

        // fill root file
        if(plane == 2){
            n_hit_energy.push_back(std::make_tuple(subRunN ,EventN, showernum, hitsize, true_hit_energy, trueEnergy, TrueHitsDep_FromTrueShower, TrueHitsDep_WithinRecoShower, total_event_pfps, this_pfp));
            fOutputTree->Fill();
        }
    }


		
    //Holder for the final product
    std::vector<double> ShowerEnergyCheater;
    ShowerEnergyCheater.push_back(trueEnergy);

    std::vector<double> EnergyError = {-999,-999,-999};
	 
    ShowerEleHolder.SetElement(ShowerEnergyCheater,EnergyError,"ShowerEnergyCheater");

	
    return 0;
    }


// Function to find only the largest shower
void ShowerEnergyCheater::FindLargestShower(std::vector<std::tuple<unsigned int, unsigned int, unsigned int, int, std::vector<double>, double, int, int, int, int>> energy_largest_shower){
    // Cut showers so we are only left with the biggest (most hits) from each event.
    // Feel like there should be a more ROOT-y way to do this...
    unsigned int i = 0;                                                                 
    while(i < energy_largest_shower.size()){  

        // If the shower at (i-1)'th position from the same event has fewer hits, delete it
        if(std::get<2>(energy_largest_shower[i]) != 0 && std::get<3>(energy_largest_shower[i]) > std::get<3>(energy_largest_shower[(i-1)])){   

            energy_largest_shower.erase(energy_largest_shower.begin() + (i-1));   
        }
        // Delete any remaining i'th non primary shower (i.e. the non primary shower which has fewer hits than the one at (i-1))
        else if(std::get<2>(energy_largest_shower[i]) != 0){ 
    
            energy_largest_shower.erase(energy_largest_shower.begin() + (i));
        }

        else{
            i++;
        }
    }


    // Add new tree to root file which has only the largest shower from each event
    art::ServiceHandle<art::TFileService> tfs; 
    TTree *energy_cheater_LS = new TTree();
    energy_cheater_LS = tfs->make<TTree>((Tree_name + "LS").c_str(), (Tree_name + "LS").c_str());

    energy_cheater_LS->Branch("Subrun", &subRunN, "Subrun/i");                                                                              
    energy_cheater_LS->Branch("Event", &EventN, "Event/i");    
    energy_cheater_LS->Branch("ShowerN", &showernum, "ShowerN/i");
    energy_cheater_LS->Branch("NHits", &hitsize, "NHits/I");                                                     
    energy_cheater_LS->Branch("Hit_Energy", &true_hit_energy);                                                     
    energy_cheater_LS->Branch("Energy", &trueEnergy, "Energy/d"); 
    energy_cheater_LS->Branch("TrueHitsDep_FromTrueShower", &TrueHitsDep_FromTrueShower, "TrueHitsDep_FromTrueShower/i"); 
    energy_cheater_LS->Branch("TrueHitsDep_WithinRecoShower", &TrueHitsDep_WithinRecoShower, "TrueHitsDep_WithinRecoShower/i"); 
    energy_cheater_LS->Branch("Total_Event_PFPs", &total_event_pfps, "Total_Event_PFPs/I");
    energy_cheater_LS->Branch("This_PFP", &this_pfp, "This_PFP/I");

    for(unsigned int i = 0; i < energy_largest_shower.size(); i++){
        subRunN    = std::get<0>(energy_largest_shower[i]);
        EventN     = std::get<1>(energy_largest_shower[i]);
        showernum  = std::get<2>(energy_largest_shower[i]);
        hitsize    = std::get<3>(energy_largest_shower[i]);
        true_hit_energy = std::get<4>(energy_largest_shower[i]);
        trueEnergy = std::get<5>(energy_largest_shower[i]);
        TrueHitsDep_FromTrueShower = std::get<7>(energy_largest_shower[i]);
        TrueHitsDep_WithinRecoShower = std::get<7>(energy_largest_shower[i]);
        total_event_pfps = std::get<8>(energy_largest_shower[i]);
        this_pfp = std::get<9>(energy_largest_shower[i]);

        energy_cheater_LS->Fill();
    }

    //m_file->Write("", TObject::kOverwrite); // save only the new version of the tree - without the arguments was duplicating original tree

}

std::map<int,std::map<geo::PlaneID,int> > ShowerEnergyCheater::NumberofPlaneHitsPerTrack(detinfo::DetectorClocksData const& clockData,
    const std::vector<art::Ptr<recob::Hit> >& hits){

  //art::ServiceHandle<geo::Geometry> geom;

  std::map<int, std::map<geo::PlaneID, int> > HitNum;

  //Loop over the hits and find the IDE
  //for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt) {
  for(auto const& h : hits){

    //art::Ptr<recob::Hit> hit = *hitIt;
    art::Ptr<recob::Hit> hit = h;

    geo::WireID wireid = hit->WireID();
    geo::PlaneID  PlaneID = wireid.planeID();

    std::vector<sim::TrackIDE> trackIDEs = bt_serv->HitToTrackIDEs(clockData, hit);

    std::map<int,float> hitEnergies;

    //Loop over the IDEs associated to the hit and add up energies
    for(unsigned int idIt = 0; idIt < trackIDEs.size(); ++idIt) {
      hitEnergies[TMath::Abs(trackIDEs.at(idIt).trackID)] += trackIDEs.at(idIt).energy;
    }
    //Find which track deposited the most energy.
    int   likelytrack = -9999;
    float MaxEnergy   = -9999;
    for(std::map<int,float>::iterator track_iter=hitEnergies.begin();track_iter!=hitEnergies.end();++track_iter){
      if(track_iter->second > MaxEnergy){
  MaxEnergy = track_iter->second;
  likelytrack = track_iter->first;
      }
    }

    ++HitNum[likelytrack][PlaneID];
  }
  return HitNum;
}

int ShowerEnergyCheater::TrueParticleID(detinfo::DetectorClocksData const& clockData,
                              const art::Ptr<recob::Hit> hit, bool rollup_unsaved_ids) {
  std::map<int,double> id_to_energy_map;
  std::vector<sim::TrackIDE> track_ides = bt_serv->HitToTrackIDEs(clockData, hit);
  for (unsigned int idIt = 0; idIt < track_ides.size(); ++idIt) {
    int id = track_ides.at(idIt).trackID;
    if (rollup_unsaved_ids) id = std::abs(id);
    double energy = track_ides.at(idIt).energy;
    id_to_energy_map[id]+=energy;
  }
  //Now loop over the map to find the maximum contributor
  double likely_particle_contrib_energy = -99999;
  int likely_track_id = 0;
  for (std::map<int,double>::iterator mapIt = id_to_energy_map.begin(); mapIt != id_to_energy_map.end(); mapIt++){
    double particle_contrib_energy = mapIt->second;
    if (particle_contrib_energy > likely_particle_contrib_energy){
      likely_particle_contrib_energy = particle_contrib_energy;
      likely_track_id = mapIt->first;
    }
  }
  return likely_track_id;
}




} // namespace 

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerEnergyCheater)
