//############################################################################
//### Name:        ShowerRecoEnergyestar                        ###
//### Author:      Tom Ham			                           ###
//### Date:        01/04/2020                                              ###
//### Description: Tool for finding the Energy of the shower by going      ###
//###              from number of hits -> number of electrons -> energy.   ###
//###              Derived from the linear energy algorithm, written for   ###
//###              the EMShower_module.cc                                  ###
//############################################################################


//Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileDirectory.h"

//LArSoft Includes
#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Tools/IShowerTool.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"

// ROOT includes
#include "TFile.h"
#include "TTree.h" 
#include "TH1D.h" 
#include "TGraph2D.h"
#include "TGraph.h"

//C++ Includes
#include <iostream>
#include <vector>
#include <tuple>

namespace ShowerRecoTools {

  class ShowerRecoEnergyestar:IShowerTool {

    public:

    ShowerRecoEnergyestar(const fhicl::ParameterSet& pset);
    
    ~ShowerRecoEnergyestar(); 
        
    //Physics Function. Calculate the shower Energy.
    int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
    art::Event& Event,
    reco::shower::ShowerElementHolder& ShowerElementHolder
    ) override;

    private:
    
    double CalculateEnergy(const detinfo::DetectorClocksData& clockData,
            const detinfo::DetectorPropertiesData& detProp, 
            std::vector<art::Ptr<recob::Hit> >& hits,
            geo::View_t& view);
    void FindLargestShower(std::vector<std::tuple<unsigned int, unsigned int, unsigned int, int, double, int>> recoenergy_largest_shower);    
    void FindPrimaryShower(std::vector<std::tuple<unsigned int, unsigned int, unsigned int, int, double, std::vector<double>, int, unsigned int, unsigned int>> recoenergy_largest_shower);

    art::InputTag fPFParticleModuleLabel;

    //Services
    detinfo::DetectorProperties const* detProp = nullptr;
    art::ServiceHandle<geo::Geometry> fGeom;
    calo::CalorimetryAlg              fCalorimetryAlg;    

    //fcl params 
    std::string fShowerEnergyestar; 
    bool fSCECorrectEField;
    double fNominalModBoxdEdx;
    double fNominalRecombinationFactor;

    // Declare variables etc.
    double Energy                = 0;
    std::vector<double> Hit_Energy;
    double nominal_Efield        = 0;
    double localEfield           = 0;
    double localEfield_cweighted = 0;
    

    unsigned int showernum = 0;
    art::SubRunNumber_t subRunN;
    art::EventNumber_t EventN;
    int hitsize;
    int SP_hits;
    int SP_hits_temp;
    int pdgcode;
    int total_event_pfps;
    int this_pfp;

    std::vector<std::tuple<int, double>> hit_key_EField_sp;

    // vec to store subrun #, event #, shower #, # of hits, energy, hit energy, pdg code, Total, PFPs, this PFP
    std::vector<std::tuple<unsigned int, unsigned int, unsigned int, int, double, std::vector<double>, int, unsigned int, unsigned int>> n_hit_energy; // more useful when making plots

    bool write_to_file = false;
    std::string File_name = "Energy_files/Cathode_sample_SCE.root"; // The file to write to
    //std::string Tree_name = "recoenergy_SCE_correction_ModBox_perhit"; 
    std::string Tree_name = "recoenergy"; 

    TFile *m_file;    ///< TFile for saving event information  
    TTree *fOutputTree;

    int n_space_points = 0;
    int n_no_space_points = 0;

    // Read in estar lookup curve
    //TFile *estar = new TFile("/sbnd/app/users/tham/estar/estar_lookup_curve_wo_density.root");
    TFile *estar = new TFile("/sbnd/app/users/tham/estar/estar_lookup_curve_wo_density.root");
    TGraph2D *t_estar = (TGraph2D*)estar->Get("g");


}; // class

    ShowerRecoEnergyestar::ShowerRecoEnergyestar(const fhicl::ParameterSet& pset) :
    IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools")),
    fPFParticleModuleLabel(pset.get<art::InputTag>("PFParticleModuleLabel","")),
//    detprop(lar::providerFrom<detinfo::DetectorPropertiesService>()),
    fCalorimetryAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg")),
    fShowerEnergyestar(pset.get<std::string>("ShowerEnergyestar")),
    fSCECorrectEField(pset.get<bool>("SCECorrectEField")),
    fNominalModBoxdEdx(pset.get<double>("NominalModBoxdEdx")),
    fNominalRecombinationFactor(pset.get<double>("NominalRecombinationFactor"))
    {

        art::ServiceHandle<art::TFileService> tfs;
   //     if(write_to_file){
            // The TFileService lets us define a tree and takes care of writing it to file
     //       m_file      = new TFile(File_name.c_str(), "UPDATE"); 
            fOutputTree = tfs->make<TTree>(Tree_name.c_str(),Tree_name.c_str());

            //add branches                                                                          
            fOutputTree->Branch("Subrun", &subRunN, "Subrun/i");
            fOutputTree->Branch("Event", &EventN, "Event/i");
            fOutputTree->Branch("ShowerN", &showernum, "ShowerN/i");
            fOutputTree->Branch("NHits", &hitsize, "NHits/I");
            fOutputTree->Branch("Energy", &Energy, "Event/d");
            fOutputTree->Branch("Hit_Energy", &Hit_Energy);
            fOutputTree->Branch("PDGcode", &pdgcode, "PDGcode/I");
            fOutputTree->Branch("Total_Event_PFPs", &total_event_pfps, "total_event_pfps/I");
            fOutputTree->Branch("This_PFP", &this_pfp, "This_PFP/I");
   //   }
    }

    ShowerRecoEnergyestar::~ShowerRecoEnergyestar()
    {
        if(write_to_file){
            //store output tree                                                                                            
            m_file->cd();                                                  
            fOutputTree->CloneTree()->Write(Tree_name.c_str(), TObject::kOverwrite);

            // Find only the largest shower from each event and write to the output file
            // To find the largest shower need to compare with the previous one and since the tool iterates over them one at a time 
            // we need to run at the end once we've looked at all the showers (hence we're in the destructor since it runs at the end).
            bool find_largest_shower = true;
            if(find_largest_shower){
               // FindLargestShower(n_hit_energy);
            }
            m_file->Close();
        }
 //       fOutputTree->CloneTree()->Write(Tree_name.c_str(),TObject::kOverwrite);
          FindPrimaryShower(n_hit_energy);
//        FindLargestShower(n_hit_energy);
    }

    int ShowerRecoEnergyestar::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
    art::Event& Event,
    reco::shower::ShowerElementHolder& ShowerEleHolder
    ){

    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Shower Reco Energy Tool ~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    
    pdgcode = pfparticle->PdgCode();
    //std::cout << pfparticle->PdgCode() << std::endl;
    // get shower number per event
    showernum = ShowerEleHolder.GetShowerNumber();

    // get subrun number
    subRunN = Event.subRun();

    // get event number 
    EventN = Event.id().event();

    //std::cout << pfparticle->Self() << std::endl;
    std::cout << "subrun: " << subRunN << "  showernum: " << showernum << "  Event: " << EventN << std::endl;

    //ShowerEleHolder.PrintElements();
    // Get the number of planes
    unsigned int numPlanes = fGeom->Nplanes();

    //Holder for the final product
    std::vector<double> ShowerRecoEnergyestar(numPlanes, -999);
    std::vector<double> ShowerRecoEnergyestarError(numPlanes, -999);

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(Event);
    auto const detProp   = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(Event, clockData);

    // Get the assocated pfParicle vertex PFParticles
    art::Handle<std::vector<recob::PFParticle> > pfpHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, pfpHandle)){
        throw cet::exception("ShowerRecoEnergyestar") << "Could not get the pandora pf particles. Something is not configured correctly Please give the correct pandora module label. Stopping";
    return 1;
    }
    std::vector<art::Ptr<recob::PFParticle>> pfps;   
    art::fill_ptr_vector(pfps, pfpHandle); 
    
    total_event_pfps = pfps.size();
    this_pfp = pfparticle->Self();
    std::cout << "Total PFPs from event: " << total_event_pfps << "  This is PFP: " << this_pfp << std::endl;

    std::map<geo::View_t, std::vector<art::Ptr<recob::Hit>>> view_hits;

    //Get the clusters
    art::Handle<std::vector<recob::Cluster> > clusHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, clusHandle)){
      throw cet::exception("ShowerRecoEnergyestar") << "Could not get the pandora clusters. Something is not configured correctly Please give the correct pandora module label. Stopping";
      return 1;
    }
    art::FindManyP<recob::Cluster> fmc(pfpHandle, Event, fPFParticleModuleLabel);
    std::vector<art::Ptr<recob::Cluster> > clusters = fmc.at(pfparticle.key());


    //Get the sapcepoints
    art::Handle<std::vector<recob::SpacePoint> > spHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, spHandle)){
        throw cet::exception("thShowerSpacePoints") << "Could not get the pandora space points. Something is not cofingured coreectly Please give the correct pandoa module     label. Stopping";
    return 1;
    }
    art::FindManyP<recob::SpacePoint> fmsp(pfpHandle, Event, fPFParticleModuleLabel);
    std::vector<art::Ptr<recob::SpacePoint> > spacepoints = fmsp.at(pfparticle.key());

    // Get E-field -  Electric Field in the drift region in KV/cm
    nominal_Efield  = detProp.Efield(); // Nominal E-field 

    // Get hit association for spacepoints
    SP_hits = 0;
    hit_key_EField_sp.clear();
    art::FindManyP<recob::Hit> fmhsp(spHandle, Event, fPFParticleModuleLabel);
    for(auto const& spacepoint : spacepoints){
        //Get the hits
        std::vector<art::Ptr<recob::Hit> > hits = fmhsp.at(spacepoint.key());
        SP_hits += hits.size();

        // get local Efield at SP_pos
        TVector3 SP_pos = IShowerTool::GetLArPandoraShowerAlg().SpacePointPosition(spacepoint); 
        double localEfield_SP = IShowerTool::GetLArPandoraShowerAlg().SCECorrectEField(nominal_Efield, SP_pos);

        for(auto const& h : hits){
            hit_key_EField_sp.push_back(std::make_tuple(h.key(), localEfield_SP));
        }

                
    }
//    std::cout << "hits associated with SP: " << SP_hits << std::endl;
 
 //   for(auto i = hit_key_EField_sp.begin(); i != hit_key_EField_sp.end(); i++){
 //       std::cout << std::get<0>(*i) << "   " << std::get<1>(*i) << std::endl;
 //   }

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
            mf::LogWarning("ShowerRecoEnergyestar") << "Shower centre calculation doesn't look to be sensible. Reconstruction is probably dodgy." << std::endl;
        }

    }

    //Get the hit association for clusters
    art::FindManyP<recob::Hit> fmhc(clusHandle, Event, fPFParticleModuleLabel);

    std::vector<std::vector<art::Ptr<recob::Hit> > > trackHits;
    trackHits.resize(numPlanes);
    
    //Loop over the clusters in the plane and get the hits
    for(auto const& cluster: clusters){
        //std::cout << "cluster key: " << cluster.key() << std::endl;
        //Get the hits
        std::vector<art::Ptr<recob::Hit> > hits = fmhc.at(cluster.key());
        std::cout << "hits associated with cluster: " << hits.size() << std::endl;

        if(hits.size() == 0){
            mf::LogWarning("ShowerRecoEnergyestar") << "No hit for the cluster. This suggest the find many is wrong."<< std::endl;
            continue;
        }

        //Get the view.
        geo::View_t view = cluster->View();

        view_hits[view].insert(view_hits[view].end(),hits.begin(),hits.end());
    }




    std::map<unsigned int, double > view_energies;
    std::vector<art::Ptr<recob::Hit>> hits;
    //Accounting for events crossing the cathode.
    for(auto const& view_hit_iter: view_hits){
 	
        hits = view_hit_iter.second;
        geo::View_t view = view_hit_iter.first;
 	
        //Calculate the Energy
        Energy = CalculateEnergy(clockData, detProp, hits,view);
        
        hitsize = hits.size();
		
        // Print out the energy for each plane
        std::cout << "Subrun: " << subRunN 
                  << " Event: " << EventN 
                  << " ShowerNum: " << showernum 
                  << " View: "<< view 
                  << "  hitsize: " << hitsize  
                  << " energy: " << Energy << std::endl;

       
        unsigned int viewNum = view;
        view_energies[viewNum] = Energy;
        
        ShowerRecoEnergyestar.at(view) = Energy;   
 
        ShowerEleHolder.SetElement(ShowerRecoEnergyestar,ShowerRecoEnergyestarError,fShowerEnergyestar);
        if(view == 2){   
        //if(write_to_file){
            fOutputTree->Fill();
        //}
        n_hit_energy.push_back(std::make_tuple(subRunN, EventN, showernum, hitsize, Energy, Hit_Energy, pdgcode, total_event_pfps, this_pfp)); //save info for collection plane
        }

    }

    //TODO think of a better way of doing this
    for (unsigned int plane=0; plane<numPlanes; ++plane) {
       
        try{
            Energy = view_energies.at(plane);
            if (Energy<0){
                mf::LogWarning("ShowerLinearEnergy") << "Negative shower energy: "<<Energy << ". Setting the energy to -999." << std::endl;
                Energy=-999;
             
             if(plane == 2){
                SP_hits_temp = 0;
                for(auto const& h : hits){
                    auto it = std::find_if(hit_key_EField_sp.begin(), hit_key_EField_sp.end(), [&](const std::tuple<unsigned int, double>& e){
                    return std::get<0>(e) == h.key();
                    });

                    if(it != hit_key_EField_sp.end()){
                        SP_hits_temp += 1;
                    }
                }
                std::cout << SP_hits_temp << std::endl;
                //n_hit_energy.push_back(std::make_tuple(subRunN ,EventN, showernum, hitsize, Energy, Hit_Energy, pdgcode, total_event_pfps, this_pfp)); //save info for collection plane
                //fOutputTree->Fill();
            }
            }  
        } 

        catch(...){
            mf::LogWarning("ShowerLinearEnergy") <<"No energy calculation for plane "<< plane << ". Setting the energy to -999." << std::endl;
            // if there's no calculation, set the energy to -999.
            Energy = -999;
            if(plane == 2){
                //n_hit_energy.push_back(std::make_tuple(subRunN, EventN, showernum, hitsize, Energy, Hit_Energy, pdgcode, total_event_pfps, this_pfp)); //save info for collection plane
                //fOutputTree->Fill();
            }
        }              	
  
    }

    if(ShowerRecoEnergyestar.size() == 0){
        throw cet::exception("ShowerLinearEnergy") << "Energy Vector is empty";
        return 1;
        }

	 
      return 0;

}


// function to calculate the reco energy
double ShowerRecoEnergyestar::CalculateEnergy(const detinfo::DetectorClocksData& clockData,
       const detinfo::DetectorPropertiesData& detProp,
       std::vector<art::Ptr<recob::Hit> >& hits, geo::View_t& view){
 
    double totalCharge = 0;
    double totalEnergy = 0;
    double nElectrons = 0;
   
    double total_electrons = 0; 
    double efield = 0; 
    //double energy_fit = 0;
   
    Hit_Energy.clear(); 
    // Loop over the hits
    for(auto const& h : hits){
        //std::cout << h.key() << std::endl; 
        if(fSCECorrectEField){
            auto it = std::find_if(hit_key_EField_sp.begin(), hit_key_EField_sp.end(), [&](const std::tuple<unsigned int, double>& e){
                return std::get<0>(e) == h.key();
            });
    
            // do a per hit SCE correction since we have associated space points
            if(it != hit_key_EField_sp.end()){
                //std::cout << "hit has SP " << std::endl;
                n_space_points += 1;
                efield = std::get<1>(*it);
            }

            // no associated space points :( so correct these hits using the shower centre position 
            else{
                //std::cout << "hit has no SP " << std::endl;    
                n_no_space_points += 1;
                efield = localEfield_cweighted;

            }
        }

        else{
            efield = 0.5; //nominal efield
            //std::cout << "Using the nominal Efield of " << nominal_Efield << "kV/cm and the nominal recombination factor of " << recombination << "." << std::endl;
        }

    //std::cout << "Using recombination factor of: " << recombination << std::endl;
    totalCharge = (h)->Integral() * fCalorimetryAlg.LifetimeCorrection(clockData, detProp, (h)->PeakTime()); // obtain charge and correct for lifetime and recombination
    nElectrons = fCalorimetryAlg.ElectronsFromADCArea(totalCharge, view);
    total_electrons += fCalorimetryAlg.ElectronsFromADCArea(totalCharge, view);
    //std::cout << nElectrons << "        " << std::endl; //t_estar->Interpolate(nElectrons, efield) << std::endl;
    totalEnergy += t_estar->Interpolate(nElectrons, efield);
    //std::cout << "hit_energy: " << t_estar->Interpolate(nElectrons, efield) << std::endl; 
    //energy_fit += -0.0503288 + (3.37871e-05 * nElectrons);
    //totalEnergy += (0.0117063 + (3.2847e-05 * nElectrons));
    //totalEnergy += t_estar->Eval(nElectrons);
    Hit_Energy.push_back(t_estar->Interpolate(nElectrons, efield));
    //std::cout << "hit_energy: " << t_estar->Interpolate(nElectrons, efield) << std::endl;
    //Hit_Energy.push_back(0.0117063 + (3.2847e-05 * nElectrons));

    }
    std::cout << "EField: " << efield << std::endl;
    return totalEnergy;
 
}

// Function to cut events where we have > 1 'major' shower
void ShowerRecoEnergyestar::FindPrimaryShower(std::vector<std::tuple<unsigned int, unsigned int, unsigned int, int, double, std::vector<double>, int, unsigned int, unsigned int>> recoenergy_largest_shower){
/*
    std::vector<std::tuple<unsigned int, unsigned int>> multiple_shower_vec;
    for(unsigned int i = 0; i < recoenergy_largest_shower.size(); i++){
        // Get all the cases where there's 2 or more showers
        if(std::get<2>(recoenergy_largest_shower[i]) == 1){
            // save the subrun and event number
            multiple_shower_vec.push_back(std::make_tuple(std::get<0>(recoenergy_largest_shower[i]), std::get<1>(recoenergy_largest_shower[i])));
        }
    }

    
    std::vector<std::tuple<unsigned int, unsigned int>> delete_vec;
    for(unsigned int j = 0; j < multiple_shower_vec.size(); j++){
        double max_hits = 0;
        for(unsigned int i = 0; i < recoenergy_largest_shower.size(); i++){
            if(std::get<0>(recoenergy_largest_shower[i]) == std::get<0>(multiple_shower_vec[j]) &&
              std::get<1>(recoenergy_largest_shower[i]) == std::get<1>(multiple_shower_vec[j])){
                // Assume first shower is (one of) the largest..
                if(std::get<2>(recoenergy_largest_shower[i]) == 0){max_hits = std::get<3>(recoenergy_largest_shower[i]);}
                else if(std::get<2>(recoenergy_largest_shower[i]) != 0 &&
                        std::get<3>(recoenergy_largest_shower[i]) > (0.05 * max_hits)){
                        delete_vec.push_back(std::make_tuple(std::get<0>(multiple_shower_vec[j]), std::get<1>(multiple_shower_vec[j])));
                }
                std::cout << max_hits << std::endl;
            }
        }
    }

    for(unsigned int j = 0; j < delete_vec.size(); j++){
        for(unsigned int i = 0; i < recoenergy_largest_shower.size();){
            if(std::get<0>(recoenergy_largest_shower[i]) == std::get<0>(delete_vec[j]) &&
               std::get<1>(recoenergy_largest_shower[i]) == std::get<1>(delete_vec[j])){
                
                recoenergy_largest_shower.erase(recoenergy_largest_shower.begin() + i);
            }
            else{i++;}
        }
    }
 */
    // Add new tree to root file which has only the largest shower from each event
    //TTree *recoenergy_cut = new TTree((Tree_name + "cut").c_str(), (Tree_name + "cut").c_str());
    TTree *recoenergy_cut = new TTree();
    art::ServiceHandle<art::TFileService> tfs;
    recoenergy_cut    = tfs->make<TTree>("recoenergy_cut","Reco Energy cut");

    recoenergy_cut->Branch("Subrun", &subRunN, "Subrun/i");                                                                              
    recoenergy_cut->Branch("Event", &EventN, "Event/i");    
    recoenergy_cut->Branch("ShowerN", &showernum, "ShowerN/i");
    recoenergy_cut->Branch("NHits", &hitsize, "NHits/I");                                                     
   // recoenergy_LS->Branch("SPHits_temp", &SP_hits_temp, "SPHits_temp/I");                                                     
    recoenergy_cut->Branch("Energy", &Energy, "Energy/d"); 
    recoenergy_cut->Branch("Hit_Energy", &Hit_Energy); 
    recoenergy_cut->Branch("pdgcode", &pdgcode, "pdgcode/i"); 
    recoenergy_cut->Branch("TotalEventPFPs", &total_event_pfps, "total_event_pfps/i"); 
    recoenergy_cut->Branch("This_PFP", &this_pfp, "this_pfp/i"); 

    for(unsigned int i = 0; i < recoenergy_largest_shower.size(); i++){
        subRunN          = std::get<0>(recoenergy_largest_shower[i]);
        EventN           = std::get<1>(recoenergy_largest_shower[i]);
        showernum        = std::get<2>(recoenergy_largest_shower[i]);
        hitsize          = std::get<3>(recoenergy_largest_shower[i]);
        Energy           = std::get<4>(recoenergy_largest_shower[i]);
        Hit_Energy       = std::get<5>(recoenergy_largest_shower[i]);
        pdgcode          = std::get<6>(recoenergy_largest_shower[i]);
        total_event_pfps = std::get<7>(recoenergy_largest_shower[i]);
        this_pfp         = std::get<8>(recoenergy_largest_shower[i]);

        recoenergy_cut->Fill();
    }
       
//}
// Function to find only the largest shower
//void ShowerRecoEnergyestar::FindLargestShower(std::vector<std::tuple<unsigned int, unsigned int, unsigned int, int, double>> recoenergy_largest_shower){


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
 //   art::ServiceHandle<art::TFileService> tfs;
    recoenergy_LS    = tfs->make<TTree>("recoenergy_LS","Reco Energy Largest Shower");

    recoenergy_LS->Branch("Subrun", &subRunN, "Subrun/i");                                                                              
    recoenergy_LS->Branch("Event", &EventN, "Event/i");    
    recoenergy_LS->Branch("ShowerN", &showernum, "ShowerN/i");
    recoenergy_LS->Branch("NHits", &hitsize, "NHits/i");                                                     
   // recoenergy_LS->Branch("SPHits_temp", &SP_hits_temp, "SPHits_temp/I");                                                     
    recoenergy_LS->Branch("Energy", &Energy, "Energy/d"); 
    recoenergy_LS->Branch("Hit_Energy", &Hit_Energy); 
    recoenergy_LS->Branch("PDGcode", &pdgcode, "PDGcode/I"); 
    recoenergy_LS->Branch("TotalEventPFPs", &total_event_pfps, "total_event_pfps/i"); 
    recoenergy_LS->Branch("This_PFP", &this_pfp, "this_pfp/i"); 

    for(unsigned int i = 0; i < recoenergy_largest_shower.size(); i++){
        subRunN          = std::get<0>(recoenergy_largest_shower[i]);
        EventN           = std::get<1>(recoenergy_largest_shower[i]);
        showernum        = std::get<2>(recoenergy_largest_shower[i]);
        hitsize          = std::get<3>(recoenergy_largest_shower[i]);
        Energy           = std::get<4>(recoenergy_largest_shower[i]);
        Hit_Energy       = std::get<5>(recoenergy_largest_shower[i]);
        pdgcode          = std::get<6>(recoenergy_largest_shower[i]);
        total_event_pfps = std::get<7>(recoenergy_largest_shower[i]);
        this_pfp         = std::get<8>(recoenergy_largest_shower[i]);

        recoenergy_LS->Fill();
    }
//    m_file->Write("", TObject::kOverwrite); // save only the new version of the tree - without the arguments was duplicating original tree
}


} // namespace

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerRecoEnergyestar)










