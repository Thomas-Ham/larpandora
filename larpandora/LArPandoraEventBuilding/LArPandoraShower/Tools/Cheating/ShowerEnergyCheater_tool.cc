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
#include "nug4/ParticleNavigation/ParticleList.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/MCCheater/BackTrackerService.h"

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

        void FindLargestShower(std::vector<std::tuple<unsigned int, unsigned int, unsigned int, int, std::vector<double>, double, int, int, int, int, double, TVector3, double, double, std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<std::vector<double>>>> energy_largest_shower);     
        int TrueParticleID(detinfo::DetectorClocksData const& clockData, const art::Ptr<recob::Hit> hit, bool rollup_unsaved_ids);
        std::map<int,std::map<geo::PlaneID,int> > NumberofPlaneHitsPerTrack(detinfo::DetectorClocksData const& clockData, const std::vector<art::Ptr<recob::Hit> >& hits);

        //Algorithm functions
        shower::LArPandoraShowerCheatingAlg fLArPandoraShowerCheatingAlg;

        //FCL
        art::InputTag fPFParticleModuleLabel;
        art::InputTag fHitModuleLabel;
        art::InputTag fShowerModuleLabel;
        art::ServiceHandle<cheat::BackTrackerService> bt_serv;
        art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
 
        unsigned int showernum = 0;
        art::SubRunNumber_t subrun;                                                                                             
        art::EventNumber_t event; 
        double trueEnergy;
        double purity;
        double completeness;
        int hitsize;
        int TrueHitsDep_WithinRecoShower;
        int TrueHitsDep_FromTrueShower;
        int total_event_pfps;
        int this_pfp;
        double theta_XZ;
        double theta_YZ;
        double true_momentum;
        TVector3 trueDir;
        std::vector<double> true_hit_energy;
        std::vector<std::vector<double>> trackide_id;
        std::vector<std::vector<double>> trackide_energy;
        std::vector<std::vector<double>> trackide_energyfrac;

        std::vector<std::tuple<unsigned int, unsigned int, unsigned int, int, std::vector<double>, double, int, int, int, int, double, TVector3, double, double, std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<std::vector<double>>>> n_hit_energy;

        std::string Tree_name = "Truth"; 
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
        fOutputTree = tfs->make<TTree>(Tree_name.c_str(),Tree_name.c_str());                                                          

        //add branches                                                                      
        fOutputTree->Branch("subrun", &subrun, "subrun/i");                                          
        fOutputTree->Branch("event", &event, "event/i");                                      
        fOutputTree->Branch("showernum", &showernum, "showernum/i");                                        
        fOutputTree->Branch("numhits", &hitsize, "numhits/I");                      
        fOutputTree->Branch("hit_energy", &true_hit_energy);                      
        fOutputTree->Branch("energy", &trueEnergy, "energy/d");   
        fOutputTree->Branch("TrueHitsDep_FromTrueShower", &TrueHitsDep_FromTrueShower, "TrueHitsDep_FromTrueShower/I");   
        fOutputTree->Branch("TrueHitsDep_WithinRecoShower", &TrueHitsDep_WithinRecoShower, "TrueHitsDep_WithinRecoShower/I");   
        fOutputTree->Branch("Total_Event_PFPs", &total_event_pfps, "Total_Event_PFPs/I");
        fOutputTree->Branch("true_momentum", &true_momentum, "true_momentum/d");
        fOutputTree->Branch("true_direction", &trueDir);
        fOutputTree->Branch("thetaXZ", &theta_XZ, "thetaXZ/d");
        fOutputTree->Branch("thetaYZ", &theta_YZ, "thetaYZ/d");
        fOutputTree->Branch("trackide_id", &trackide_id);
        fOutputTree->Branch("trackide_energy", &trackide_energy);
        fOutputTree->Branch("trackide_energyfrac", &trackide_energyfrac);

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
    subrun = Event.subRun();
    // get event number 
    event = Event.id().event();

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

    // Get the spacepoints
    art::Handle<std::vector<recob::SpacePoint>> spHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, spHandle)){
      throw cet::exception("ShowerEnergyCheater") << "Could not get the pandora spacepoints. Something is not cofingured coreectly Please give the correct pandoa module label. Stopping";
      return 1;
    }

    art::FindManyP<recob::SpacePoint> fmsp(pfpHandle, Event, fPFParticleModuleLabel);
    std::vector<art::Ptr<recob::SpacePoint>> spacepoints = fmsp.at(pfparticle.key());

    // Get the hits association from the spacepoints
    art::FindManyP<recob::Hit> fmhsp(spHandle, Event, fPFParticleModuleLabel);

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
   
    std::vector<double> sp_hits_vec;
    sp_hits_vec.clear(); 
    for(auto const& spacepoint : spacepoints){
      //Get the hits
      std::vector<art::Ptr<recob::Hit>> hits_sp{fmhsp.at(spacepoint.key())};
      for(auto h : hits_sp){
        sp_hits_vec.push_back(h.key());
      }

      // get local Efield at SP_pos
      //TVector3 sp_pos{IShowerTool::GetLArPandoraShowerAlg().SpacePointPosition(spacepoint)}; 
      //const double localEfield_SP{IShowerTool::GetLArPandoraShowerAlg().SCECorrectEField(fNominal_Efield, sp_pos)};
    }
   

    geo::PlaneID::PlaneID_t bestPlane{std::numeric_limits<geo::PlaneID::PlaneID_t>::max()};
    unsigned int bestPlaneSubrun = 0;
    unsigned int bestPlaneEvent = 0;
    unsigned int bestPlaneShowernum = 0;
    int bestPlaneNumHits = 0;
    std::vector<double>  bestPlaneHitEnergy;
    double bestPlaneEnergy = 0;

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


            trueDir = {trueParticle->Px(),trueParticle->Py(),trueParticle->Pz()};
            true_momentum = trueParticle->P();
            
            //std::cout << trueParticle->NumberTrajectoryPoints() << std::endl; 
            //std::cout << trueParticle->Momentum(trueParticle->NumberTrajectoryPoints() - 2).T() << std::endl;
            //std::cout << trueParticle->Momentum(trueParticle->NumberTrajectoryPoints() - 2).X() << std::endl;
            //std::cout << trueParticle->Momentum(trueParticle->NumberTrajectoryPoints() - 2).Y() << std::endl;
            //std::cout << trueParticle->Momentum(trueParticle->NumberTrajectoryPoints() - 2).Z() << std::endl;

            //std::cout << true_momentum << std::endl;
            //std::cout << trueDir[0] << std::endl;
            //std::cout << trueDir[1] << std::endl;
            //std::cout << trueDir[2] << std::endl;
            
            TLorentzVector pos = trueParticle->Position();
            std::cout << "Pos: " << pos.X() << "  " << pos.Y() << "  " << pos.Z() << std::endl;

            // XZ plane
            double argXZ = trueDir[0]/true_momentum;
            if(argXZ > 1 && argXZ < 1.000001){ // stupid floating point precision causing issues 
              argXZ = 1;
            }

            if(trueDir[0] >= 0 && trueDir[2] >= 0){
              theta_XZ = 90 - (acos(argXZ)*(180/M_PI));
            }
            else if(trueDir[0] >= 0 && trueDir[2] < 0){
              theta_XZ = 90 + (acos(argXZ)*(180/M_PI));
            }
            else if(trueDir[0] < 0 && trueDir[2] > 0){
              theta_XZ = -(90 - (acos(abs(argXZ))*(180/M_PI)));
            }
            else if(trueDir[0] < 0 && trueDir[2] < 0){
              theta_XZ = -(90 + (acos(abs(argXZ))*(180/M_PI)));
            }

            // YZ plane
            if(trueDir[1] >= 0 && trueDir[2] >= 0){
              theta_YZ = 90 - (acos(trueDir[1]/true_momentum)*(180/M_PI));
            }
            else if(trueDir[1] >= 0 && trueDir[2] < 0){
              theta_YZ = 90 + (acos(trueDir[1]/true_momentum)*(180/M_PI));
            }
            else if(trueDir[1] < 0 && trueDir[2] > 0){
              theta_YZ = -(90 - (acos(abs(trueDir[1]/true_momentum))*(180/M_PI)));
            }
            else if(trueDir[1] < 0 && trueDir[2] < 0){
              theta_YZ = -(90 + (acos(abs(trueDir[1]/true_momentum))*(180/M_PI)));
            }


        }
    
        std::cout << "Plane: " << plane << std::endl;   
        std::cout << "Hits: " << hitsize << std::endl;
        std::cout << "TrueHitsDep_WithinRecoShower: " << TrueHitsDep_WithinRecoShower << std::endl;       
        std::cout << "TrueHitsDep_FromTrueShower: " << TrueHitsDep_FromTrueShower << std::endl;
        std::cout << "purity: " << purity << std::endl;
        std::cout << "completeness: " << completeness << std::endl;
        std::cout << "True Energy: " << trueEnergy << std::endl;
        std::cout << "True momentum: " << true_momentum << std::endl;
        std::cout << "True direction: (" << trueDir[0] << ", " << trueDir[1] << ", " << trueDir[2] << ")" << std::endl;
        std::cout << "ThetaXZ: " << theta_XZ << std::endl;
        std::cout << "ThetaYZ: " << theta_YZ << std::endl;

   
        // make a collection of the distinct eve ID values
        std::set<int> eveIDs; 

        // Get true hit energy
        double deposit = 0;
        //for(auto const& h_sp : sp_hits_vec){
        trackide_id.clear();
        trackide_energy.clear();
        trackide_energyfrac.clear();
        for(auto const& h : hits){
          //if( h.key() == h_sp){
            

            std::vector<sim::TrackIDE> trackIDEs = bt_serv->HitToTrackIDEs(clockData, h);
            std::vector<sim::TrackIDE> eveides = bt_serv->HitToEveTrackIDEs(clockData, h);

            std::vector<sim::TrackIDE> ides;
            for( auto ide : trackIDEs ){
                ides.push_back(ide);
            }

            std::vector<double> temp_ide_id_vec;   
            std::vector<double> temp_ide_energy_vec;   
            std::vector<double> temp_ide_energyfrac_vec;   
            temp_ide_id_vec.clear();
            temp_ide_energy_vec.clear();
            temp_ide_energyfrac_vec.clear();
           
            for (size_t t = 0; t < trackIDEs.size(); ++t) {
                // find the Eve particle for the current trackID
                int eveID = pi_serv->ParticleList().EveId(trackIDEs[t].trackID);

                std::cout  << "track id: " << trackIDEs[t].trackID << " contributed " << trackIDEs[t].energy << "/"
                           << trackIDEs[t].energyFrac << " to the current hit and has eveID: " << eveID << std::endl;
                
                temp_ide_id_vec.push_back(trackIDEs[t].trackID);
                temp_ide_energy_vec.push_back(trackIDEs[t].energy);
                temp_ide_energyfrac_vec.push_back(trackIDEs[t].energyFrac);

                /*
                const simb::MCParticle *particle = pi_serv->TrackIdToParticle_P(trackIDEs[t].trackID);
                if (particle){
                  std::cout<<"Pdgcode = "<<particle->PdgCode()<<std::endl;
                }
                */
            }
            trackide_id.push_back(temp_ide_id_vec);
            trackide_energy.push_back(temp_ide_energy_vec);
            trackide_energyfrac.push_back(temp_ide_energyfrac_vec);

            /*
            for (size_t e = 0; e < eveides.size(); ++e) {
              std::cout << "eve id: " << eveides[e].trackID << " contributed " << eveides[e].energy << "/"
                << eveides[e].energyFrac << " to the current hit" << std::endl;
            }
            */




            double hit_energy = 0;
            for (std::vector<sim::TrackIDE>::iterator ideIt = ides.begin(); ideIt != ides.end(); ++ideIt){
                hit_energy += ideIt->energy;
                deposit += ideIt->energy;
            }
            true_hit_energy.push_back(hit_energy);
            //if(h->WireID().Wire == 1131 || h->WireID().Wire == 1163 || h->WireID().Wire == 1161){
              std::cout << "WireID: "  << h->WireID() << "  PeakTime: " << h->PeakTime() <<  "  hit energy: " << std::setprecision (15) << hit_energy << std::endl;
            //}
          //} 
        //}
        }
        std::cout << "sum of hits energy: " << deposit << std::endl;
    
        if(hitsize > bestPlaneNumHits){
              bestPlane        = plane;
              bestPlaneNumHits = hitsize;
              bestPlaneSubrun  = subrun;
              bestPlaneEvent   = event;
              bestPlaneShowernum = showernum;
              bestPlaneEnergy  = trueEnergy;
              bestPlaneHitEnergy = true_hit_energy;
        }


        // fill root file
        
        if(plane == 1){
            n_hit_energy.push_back(std::make_tuple(subrun ,event, showernum, hitsize, true_hit_energy, trueEnergy, TrueHitsDep_FromTrueShower, TrueHitsDep_WithinRecoShower, total_event_pfps, this_pfp, true_momentum, trueDir, theta_XZ, theta_YZ, trackide_id, trackide_energy, trackide_energyfrac));
            fOutputTree->Fill();
        }
        
    }


    // Only set the best plane if it has some hits in it
    if (bestPlane < geom->Nplanes()){
      // Need to cast as an int for legacy default of -999
      // have to define a new variable as we pass-by-reference when filling
      mf::LogVerbatim("myShowerESTAREnergy") << "best plane: " << bestPlane << 
                                                "  best plane subrun: " << bestPlaneSubrun <<
                                                "  best plane event: " << bestPlaneEvent <<
                                                "  best plane showernum: " << bestPlaneShowernum <<
                                                "  best plane number of hit: " << bestPlaneNumHits <<
                                                "  best energy: " << bestPlaneEnergy << std::endl;
                                                //"  best hit energy: " << bestPlaneHitEnergy << std::endl;
      
      //n_hit_energy.push_back(std::make_tuple(bestPlaneSubrun, bestPlaneEvent, bestPlaneShowernum, bestPlaneNumHits, bestPlaneHitEnergy, bestPlaneEnergy, TrueHitsDep_FromTrueShower, TrueHitsDep_WithinRecoShower, total_event_pfps, this_pfp, true_momentum, trueDir, theta_XZ, theta_YZ));
      //fOutputTree->Fill();
    }

    
    //Holder for the final product
    std::vector<double> ShowerEnergyCheater;
    ShowerEnergyCheater.push_back(trueEnergy);

    std::vector<double> EnergyError = {-999,-999,-999};
   
    ShowerEleHolder.SetElement(ShowerEnergyCheater,EnergyError,"ShowerEnergyCheater");

  
    return 0;
    }


// Function to find only the largest shower
void ShowerEnergyCheater::FindLargestShower(std::vector<std::tuple<unsigned int, unsigned int, unsigned int, int, std::vector<double>, double, int, int, int, int, double, TVector3, double, double, std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<std::vector<double>>>> energy_largest_shower){
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
    energy_cheater_LS = tfs->make<TTree>((Tree_name + "_LS").c_str(), (Tree_name + "_LS").c_str());

    energy_cheater_LS->Branch("subrun", &subrun, "subrun/i");                                                                              
    energy_cheater_LS->Branch("event", &event, "event/i");    
    energy_cheater_LS->Branch("showernum", &showernum, "showernum/i");
    energy_cheater_LS->Branch("numhits", &hitsize, "numhits/I");                                                     
    energy_cheater_LS->Branch("hit_energy", &true_hit_energy);                                                     
    energy_cheater_LS->Branch("energy", &trueEnergy, "energy/d"); 
    energy_cheater_LS->Branch("TrueHitsDep_FromTrueShower", &TrueHitsDep_FromTrueShower, "TrueHitsDep_FromTrueShower/i"); 
    energy_cheater_LS->Branch("TrueHitsDep_WithinRecoShower", &TrueHitsDep_WithinRecoShower, "TrueHitsDep_WithinRecoShower/i"); 
    energy_cheater_LS->Branch("Total_Event_PFPs", &total_event_pfps, "Total_Event_PFPs/I");
    energy_cheater_LS->Branch("This_PFP", &this_pfp, "This_PFP/I");
    energy_cheater_LS->Branch("true_momentum", &true_momentum, "true_momentum/d");
    energy_cheater_LS->Branch("true_direction", &trueDir);
    energy_cheater_LS->Branch("thetaXZ", &theta_XZ, "thetaXZ/d");
    energy_cheater_LS->Branch("thetaYZ", &theta_YZ, "thetaYZ/d");
    energy_cheater_LS->Branch("trackide_id", &trackide_id);
    energy_cheater_LS->Branch("trackide_energy", &trackide_energy);
    energy_cheater_LS->Branch("trackide_energyfrac", &trackide_energyfrac);

    for(unsigned int i = 0; i < energy_largest_shower.size(); i++){
        subrun    = std::get<0>(energy_largest_shower[i]);
        event     = std::get<1>(energy_largest_shower[i]);
        showernum  = std::get<2>(energy_largest_shower[i]);
        hitsize    = std::get<3>(energy_largest_shower[i]);
        true_hit_energy = std::get<4>(energy_largest_shower[i]);
        trueEnergy = std::get<5>(energy_largest_shower[i]);
        TrueHitsDep_FromTrueShower = std::get<7>(energy_largest_shower[i]);
        TrueHitsDep_WithinRecoShower = std::get<7>(energy_largest_shower[i]);
        total_event_pfps = std::get<8>(energy_largest_shower[i]);
        this_pfp = std::get<9>(energy_largest_shower[i]);
        true_momentum = std::get<10>(energy_largest_shower[i]);
        trueDir = std::get<11>(energy_largest_shower[i]);
        theta_XZ = std::get<12>(energy_largest_shower[i]);
        theta_YZ = std::get<13>(energy_largest_shower[i]);
        trackide_id = std::get<14>(energy_largest_shower[i]);
        trackide_energy = std::get<15>(energy_largest_shower[i]);
        trackide_energyfrac = std::get<16>(energy_largest_shower[i]);

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
