//############################################################################
//### Name:        myShowerESTAREnergy                                       ###
//### Author:      Tom Ham                                                 ###
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
#include "lardataobj/AnalysisBase/T0.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"

// ROOT includes
#include "TFile.h"
#include "TGraph2D.h"

namespace ShowerRecoTools {

  class myShowerESTAREnergy:IShowerTool {

  public:

  myShowerESTAREnergy(const fhicl::ParameterSet& pset);
    
  ~myShowerESTAREnergy(); 
        
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

  void FindLargestShower(std::vector<std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, std::vector<double>, double,
      int, int, std::vector<double>, // sphits, nosphits, Goodness of fit
      std::vector<double>, std::vector<double>, std::vector<double>,        // Peak time, RMS, PeakAmplitude
      std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>        // SummedADC, Integral, Multiplicity, charge
      >> energy_largest_shower);     

  art::InputTag fPFParticleLabel;

  //Services
  art::ServiceHandle<geo::Geometry> fGeom;
  calo::CalorimetryAlg              fCalorimetryAlg;    

  //fcl params 
  std::string fShowerEnergyOutputLabel;
  std::string fShowerBestPlaneOutputLabel; 
  bool fSCECorrectEField;

  double fNominal_Efield;          // Nominal E-field in the detector
  double fLocal_Efield;            // Local E-field at the geometric centre of the shower 
  double fLocal_Efield_cweighted;  // Local E-field at the charge weighted centre of the shower

  std::vector<std::tuple<int, double>> sp_hits_tuple; // tuple holding the key and local E-field for all the sp hits 

  std::string fname = "/sbnd/app/users/tham/estar/ESTAREnergyLookupCurve.root";
  std::string ESTARPath;
  std::string ESTAR_TGraph_name;
  TGraph2D *t_estar;

  std::vector<double> hit_energy;
  std::vector<double> hit_goodness_of_fit;
  std::vector<double> hit_peak_time;
  std::vector<double> hit_rms;
  std::vector<double> hit_peak_amplitude;
  std::vector<double> hit_summed_adc;
  std::vector<double> hit_integral;
  std::vector<double> hit_multiplicity;
  std::vector<double> hit_charge;
  unsigned int hitsize;
  double energy;
  unsigned int showernum;
  art::SubRunNumber_t subrun;                                                                                             
  art::EventNumber_t event;
  
  int sphits{0}, nosphits{0};

  TTree *fOutputTree;
  std::vector<std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, std::vector<double>, double,
    int, int, std::vector<double>,
    std::vector<double>, std::vector<double>, std::vector<double>,        // Peak time, RMS, PeakAmplitude
    std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>        // SummedADC, Integral, Multiplicity, charge
    >> energy_vec;

  }; // class

  myShowerESTAREnergy::myShowerESTAREnergy(const fhicl::ParameterSet& pset) :
    IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools")),
    fPFParticleLabel(pset.get<art::InputTag>("PFParticleLabel")),
    fCalorimetryAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg")),
    fShowerEnergyOutputLabel(pset.get<std::string>("ShowerEnergyOutputLabel")),
    fShowerBestPlaneOutputLabel(pset.get<std::string>("ShowerBestPlaneOutputLabel")),
    fSCECorrectEField(pset.get<bool>("SCECorrectEField")),
    fNominal_Efield(0),
    fLocal_Efield(0),
    fLocal_Efield_cweighted(0),
    ESTARPath(pset.get<std::string>("ESTARFname")),
    ESTAR_TGraph_name(pset.get<std::string>("ESTARTGraph"))
  {

    /*
    //Get the ESTAR lookup curve file name
    cet::search_path sp("FW_SEARCH_PATH");
    if (!sp.find_file(ESTARPath, fname)) {
      std::cout << ESTARPath << std::endl;
      throw cet::exception("myShowerESTAREnergy") << "Could not find the ESTAR lookup curve file.";
    }
    */
    
    TFile fin(fname.c_str(), "READ");
    if (!fin.IsOpen()) {
      throw cet::exception("myShowerESTAREnergy") << "Could not read the ESTAR file. Stopping";
    }

    // Get the TGraph.
    t_estar = dynamic_cast<TGraph2D*>(fin.Get(ESTAR_TGraph_name.c_str()));
    if(!t_estar){
      throw cet::exception("myShowerESTAREnergy") << "Could not read the ESTAR TGraph";
    }

    // save info 
    art::ServiceHandle<art::TFileService> tfs; 
    fOutputTree = tfs->make<TTree>("Reco","Reco");  
    fOutputTree->Branch("subrun", &subrun);    
    fOutputTree->Branch("event", &event);    
    fOutputTree->Branch("showernum", &showernum);    
    fOutputTree->Branch("numhits", &hitsize);    
    fOutputTree->Branch("hit_energy", &hit_energy);    
    fOutputTree->Branch("energy", &energy);    
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

  }

  myShowerESTAREnergy::~myShowerESTAREnergy()
  {
    bool find_largest_shower = true;                                                                                         
    if(find_largest_shower){                                                                            
      FindLargestShower(energy_vec);                                                                                                    
    }  
  }

  int myShowerESTAREnergy::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
        art::Event& Event,
        reco::shower::ShowerElementHolder& ShowerEleHolder
  ){

    mf::LogInfo("myShowerESTAREnergy") << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Shower ESTAR Energy Tool ~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    if(fSCECorrectEField){mf::LogVerbatim("myShowerESTAREnergy") << "Applying SCE correction" << std::endl;}

    showernum = ShowerEleHolder.GetShowerNumber();
    subrun = Event.subRun();
    event = Event.id().event();

    // Get the associated pfParicle vertex PFParticles
    auto const pfpHandle{Event.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleLabel)};

    // Get the clusters
    auto const clusHandle{Event.getValidHandle<std::vector<recob::Cluster>>(fPFParticleLabel)};

    const art::FindManyP<recob::Cluster>& fmc{ShowerEleHolder.GetFindManyP<recob::Cluster>(pfpHandle, Event, fPFParticleLabel)};
    std::vector<art::Ptr<recob::Cluster>> clusters{fmc.at(pfparticle.key())};

    // Get the hits association from the clusters
    const art::FindManyP<recob::Hit>& fmhc{ShowerEleHolder.GetFindManyP<recob::Hit>(clusHandle, Event, fPFParticleLabel)};
    
    std::map<geo::PlaneID::PlaneID_t, std::vector<art::Ptr<recob::Hit>>> planeHits;

    // Get the spacepoints
    auto const spHandle{Event.getValidHandle<std::vector<recob::SpacePoint>>(fPFParticleLabel)};

    const art::FindManyP<recob::SpacePoint>& fmsp{ShowerEleHolder.GetFindManyP<recob::SpacePoint>(pfpHandle, Event, fPFParticleLabel)};
    std::vector<art::Ptr<recob::SpacePoint>> spacepoints{fmsp.at(pfparticle.key())};

    // Get the hits association from the spacepoints
    const art::FindManyP<recob::Hit>& fmhsp{ShowerEleHolder.GetFindManyP<recob::Hit>(spHandle, Event, fPFParticleLabel)};

    //Get T0-PFParticle association
    double pfpT0Time = 0;
    const art::FindManyP<anab::T0>& fmpfpt0 = ShowerEleHolder.GetFindManyP<anab::T0>(pfpHandle, Event, fPFParticleLabel);
    std::vector<art::Ptr<anab::T0> > pfpT0Vec = fmpfpt0.at(pfparticle.key());
    if (pfpT0Vec.size()==1) {
     pfpT0Time = pfpT0Vec.front()->Time();
    }
    std::cout << "T0 time: " << pfpT0Time << std::endl;


    auto const clockData{art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(Event)};
    auto const detProp{art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(Event, clockData)};
    
    // Get E-field -  Electric Field in the drift region in KV/cm
    fNominal_Efield = detProp.Efield(); // Nominal E-field 
    std::cout << "nom efield: " << fNominal_Efield << std::endl;
   
    // Not all hits from the clusters will match to Space points. If that's the case let's use some sort of shower centre
    // Get shower centre 
    const TVector3 showercentre{IShowerTool::GetLArPandoraShowerAlg().ShowerCentre(spacepoints)};
    
    // Get the charge weighted shower centre - sum(charge*position)/sum(charge)
    // Charge is obtained from the lifetime corrected hits and hits are obtained from the spacepoints 
    const TVector3 chargecentre{IShowerTool::GetLArPandoraShowerAlg().ShowerCentre(clockData, detProp, spacepoints, fmhsp)};
    std::cout << "charge centre: " << chargecentre.X() << "  " << chargecentre.Y() << "  " << chargecentre.Z() << std::endl;

    // Get space charge corrected E-field
    if(fSCECorrectEField){
      // Check the shower centre exists and is sensible
      // Think the charge weighted centre can take nan values which causes the E-field calculation to fail
      if(!isnan(showercentre.Mag()) && !isnan(chargecentre.Mag())){ 
        fLocal_Efield = IShowerTool::GetLArPandoraShowerAlg().SCECorrectEField(fNominal_Efield, showercentre);
        fLocal_Efield_cweighted = IShowerTool::GetLArPandoraShowerAlg().SCECorrectEField(fNominal_Efield, chargecentre);
      }
      else{
        mf::LogWarning("myShowerESTAREnergy") << "Shower centre calculation doesn't look to be sensible. Reconstruction is probably dodgy." << std::endl;
      }
    }

    std::cout << "cweighted efield: " << fLocal_Efield_cweighted << std::endl;
    // Loop over the clusters in the plane and get the hits
    for(auto const& cluster: clusters){
      std::cout << "cluster hits: " << cluster->NHits() << std::endl;       
      std::cout << "cluster charge: " << cluster->Integral() << std::endl;       
      std::cout << "cluster charge: " << fCalorimetryAlg.ElectronsFromADCArea(cluster->Integral(), cluster->Plane().Plane) << std::endl;       
      // Get the hits
      const std::vector<art::Ptr<recob::Hit>> hits{fmhc.at(cluster.key())};
      // Get the plane.
      const geo::PlaneID::PlaneID_t plane{cluster->Plane().Plane};
       
      planeHits[plane].insert(planeHits[plane].end(),hits.begin(),hits.end());
    }

    // Loop over the spacepoints and get the hits
    sp_hits_tuple.clear();
    for(auto const& spacepoint : spacepoints){
      //Get the hits
      const std::vector<art::Ptr<recob::Hit>> hits{fmhsp.at(spacepoint.key())};

      // get local Efield at SP_pos
      TVector3 sp_pos{IShowerTool::GetLArPandoraShowerAlg().SpacePointPosition(spacepoint)}; 
      const double localEfield_SP{IShowerTool::GetLArPandoraShowerAlg().SCECorrectEField(fNominal_Efield, sp_pos)};

      // Store the localEfield of all the hits with associated SP's so we can compare with the hits from the clusters later
      for(auto const& h : hits){
        sp_hits_tuple.push_back(std::make_tuple(h.key(), localEfield_SP));
      }           
    }

    // Calculate the energy for each plane && best plane
    hit_energy.clear();
    hit_goodness_of_fit.clear();
    hit_peak_time.clear();
    hit_rms.clear();
    hit_peak_amplitude.clear();
    hit_summed_adc.clear();
    hit_integral.clear();
    hit_multiplicity.clear();
    hit_charge.clear();
    geo::PlaneID::PlaneID_t bestPlane{std::numeric_limits<geo::PlaneID::PlaneID_t>::max()};
    unsigned int bestPlaneSubrun = 0;
    unsigned int bestPlaneEvent = 0;
    unsigned int bestPlaneShowernum = 0;
    unsigned int bestPlaneNumHits = 0;
    std::vector<double>  bestPlaneHitEnergy;
    double bestPlaneEnergy = 0;

    //Holder for the final product
    std::vector<double> myShowerESTAREnergy(fGeom->Nplanes(), -999);
    std::vector<double> myShowerESTAREnergyError(fGeom->Nplanes(), -999);

    for(auto const& [plane, hits] : planeHits){

      unsigned int planeNumHits = hits.size();

      //Calculate the Energy
      energy = CalculateEnergy(clockData, detProp, hits, plane);

      if(energy > std::numeric_limits<double>::epsilon()){    
        myShowerESTAREnergy.at(plane) = energy;
      }

      if(planeNumHits > bestPlaneNumHits){
        bestPlane        = plane;
        bestPlaneNumHits = planeNumHits;
        bestPlaneSubrun  = subrun;
        bestPlaneEvent   = event;
        bestPlaneShowernum = showernum;
        bestPlaneEnergy  = energy;
        bestPlaneHitEnergy = hit_energy;
      }
        
      mf::LogVerbatim("myShowerESTAREnergy") 
                  << "plane: " << plane
                  << "  hits: " << hits.size()
                  << "  energy: " << energy << std::endl;  
      // fill root file
      
      if(plane == 2){
        energy_vec.push_back(std::make_tuple(subrun, event, showernum, planeNumHits, hit_energy, energy,
        sphits, nosphits, hit_goodness_of_fit,
        hit_peak_time, hit_rms, hit_peak_amplitude,
        hit_summed_adc, hit_integral, hit_multiplicity, hit_charge));
        fOutputTree->Fill();
      }
      
    }


    ShowerEleHolder.SetElement(myShowerESTAREnergy,myShowerESTAREnergyError,fShowerEnergyOutputLabel);

    // Only set the best plane if it has some hits in it
    if (bestPlane < fGeom->Nplanes()){
      // Need to cast as an int for legacy default of -999
      // have to define a new variable as we pass-by-reference when filling
      int bestPlaneVal(bestPlane);
      ShowerEleHolder.SetElement(bestPlaneVal, fShowerBestPlaneOutputLabel);
      mf::LogVerbatim("myShowerESTAREnergy") << "best plane: " << bestPlane << 
                                                "  best plane subrun: " << bestPlaneSubrun <<
                                                "  best plane event: " << bestPlaneEvent <<
                                                "  best plane shoernum: " << bestPlaneShowernum <<
                                                "  best plane number of hit: " << bestPlaneNumHits <<
                                                "  best energy: " << bestPlaneEnergy << std::endl;
      
      //energy_vec.push_back(std::make_tuple(bestPlaneSubrun, bestPlaneEvent, bestPlaneShowernum, bestPlaneNumHits, bestPlaneHitEnergy, bestPlaneEnergy));
      //fOutputTree->Fill();
    }
      
    return 0;
  }


  // function to calculate the reco energy
  double myShowerESTAREnergy::CalculateEnergy(const detinfo::DetectorClocksData& clockData,
         const detinfo::DetectorPropertiesData& detProp,
         const std::vector<art::Ptr<recob::Hit> >& hits, const geo::PlaneID::PlaneID_t plane){
   
    double totalEnergy{0};
    double efield{0};
    double hit_charge_sum{0};
    //bool SP = false; 
    sphits = 0; nosphits = 0;
    // Loop over the hits
    for(auto const& h : hits){
      // If we're correcting for SCE, check if the hit is one that was we found to have an associated SP
      if(fSCECorrectEField){
        const auto it = std::find_if(sp_hits_tuple.begin(), sp_hits_tuple.end(), [&](const std::tuple<unsigned int, double>& e){
        return std::get<0>(e) == h.key();
        });
      
        // do a per hit SCE correction since we have associated space points
        if(it != sp_hits_tuple.end()){
          efield = std::get<1>(*it);
          //efield = fNominal_Efield;
          ++sphits;
          //SP = true;
        }

        // no associated space points so correct these hits using the shower centre position 
        else{
          efield = fLocal_Efield_cweighted;
          ++nosphits;
          //SP = false;
          continue; 
        }
      }

      else{
        efield = fNominal_Efield; 
      }

      const double totalCharge{(h)->Integral() * fCalorimetryAlg.LifetimeCorrection(clockData, detProp, (h)->PeakTime())}; // obtain charge and correct for lifetime
      //const double totalCharge{(h)->Integral()}; // obtain charge and correct for lifetime
      const double nElectrons{fCalorimetryAlg.ElectronsFromADCArea(totalCharge, plane)};
      if(nElectrons !=0 && t_estar->Interpolate(nElectrons, efield) <= 0){
        mf::LogWarning("myShowerESTAREnergy") << "Interpolation is returning an energy of 0 (or negative) even though we do have some collected charge. " 
                                            << "Probably due to an out of bounds entry for the charge or the efield." << std::endl; 
      }
      // Make sure we aren't counting energy which is somehow negative
      if(t_estar->Interpolate(nElectrons, efield) > 0){
        //std::cout << "hit E: " << t_estar->Interpolate(nElectrons, efield) << "   SP: " << SP << std::endl;
        totalEnergy += t_estar->Interpolate(nElectrons, efield);
        hit_charge_sum += nElectrons;
        if(plane == 2){
          hit_energy.push_back(t_estar->Interpolate(nElectrons, efield));
          hit_goodness_of_fit.push_back(h->GoodnessOfFit());
          hit_peak_time.push_back(h->PeakTime());
          hit_rms.push_back(h->RMS());
          hit_peak_amplitude.push_back(h->PeakAmplitude());
          hit_summed_adc.push_back(h->SummedADC());
          hit_integral.push_back(h->Integral());
          hit_multiplicity.push_back(h->Multiplicity());
          hit_charge.push_back(nElectrons);

          //std::cout << t_estar->Interpolate(nElectrons, efield) << std::endl;
          //if(h->WireID().Wire == 1131 || h->WireID().Wire == 1163 || h->WireID().Wire == 1164){
          
          //if(subrun == 49 && event == 100){  // && t_estar->Interpolate(nElectrons, efield)  > 9.7402 && t_estar->Interpolate(nElectrons, efield) < 9.7403) {
            //std::cout << "subrun: " << subrun << "  Event: " << event << std::endl;

            std::cout << "channel: " << h->Channel() << "   WireID: " << h->WireID() << "  PeakTime: "<< h->PeakTime() << "  RMS: " << h->RMS() << 
                         "  nElectrons: " << nElectrons << "  Hit Energy: " << t_estar->Interpolate(nElectrons, efield) << " GoodnessofFit: " << h->GoodnessOfFit() <<
                        "  Start tick: " << h->StartTick() << "   End tick: " << h->EndTick() << "  Integral: " << h->Integral() << " PeakAmp: " << h->PeakAmplitude() << std::endl;
            //std::cout << t_estar->Interpolate(nElectrons, efield) << std::endl;
            //std::cin.get();
          //}
          
        }
      }
    }
    std::cout << "hit_charge_sum: " << hit_charge_sum << std::endl;
    std::cout << "hit_charge_sum Energy: " << t_estar->Interpolate(hit_charge_sum, efield) << std::endl;
    //std::cout << "sphit: " << sphits << std::endl;
    //std::cout << "nosphit: " << nosphits << std::endl;

    return totalEnergy; 
  }

 // Function to find only the largest shower
void myShowerESTAREnergy::FindLargestShower(std::vector<std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, std::vector<double>, double,
    int, int, std::vector<double>,
    std::vector<double>, std::vector<double>, std::vector<double>,        // Peak time, RMS, PeakAmplitude
    std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>        // SummedADC, Integral, Multiplicity
    >> energy_largest_shower){
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
    TTree *energy_LS = new TTree();
    energy_LS = tfs->make<TTree>(("Reco_LS"), ("Reco_LS"));

    energy_LS->Branch("subrun", &subrun, "subrun/i");                                                                              
    energy_LS->Branch("event", &event, "event/i");    
    energy_LS->Branch("showernum", &showernum, "showernum/i");
    energy_LS->Branch("numhits", &hitsize, "numhits/i");                                                     
    energy_LS->Branch("hit_energy", &hit_energy);                                                     
    energy_LS->Branch("energy", &energy, "energy/d"); 
    energy_LS->Branch("sphits", &sphits, "sphits/i"); 
    energy_LS->Branch("nosphits", &nosphits, "nosphits/i"); 
    energy_LS->Branch("hit_goodness_of_fit", &hit_goodness_of_fit);                                                     
    energy_LS->Branch("hit_peak_time", &hit_peak_time);                                                     
    energy_LS->Branch("hit_rms", &hit_rms);                                                     
    energy_LS->Branch("hit_peak_amplitude", &hit_peak_amplitude);                                                     
    energy_LS->Branch("hit_summed_adc", &hit_summed_adc);                                                     
    energy_LS->Branch("hit_integral", &hit_integral);                                                     
    energy_LS->Branch("hit_multiplicity", &hit_multiplicity);                                                     
    energy_LS->Branch("hit_charge", &hit_charge);                                                     

    for(unsigned int i = 0; i < energy_largest_shower.size(); i++){
        subrun    = std::get<0>(energy_largest_shower[i]);
        event     = std::get<1>(energy_largest_shower[i]);
        showernum  = std::get<2>(energy_largest_shower[i]);
        hitsize    = std::get<3>(energy_largest_shower[i]);
        hit_energy = std::get<4>(energy_largest_shower[i]);
        energy = std::get<5>(energy_largest_shower[i]);
        sphits = std::get<6>(energy_largest_shower[i]);
        nosphits = std::get<7>(energy_largest_shower[i]);
        hit_goodness_of_fit = std::get<8>(energy_largest_shower[i]);
        hit_peak_time = std::get<9>(energy_largest_shower[i]);
        hit_rms = std::get<10>(energy_largest_shower[i]);
        hit_peak_amplitude = std::get<11>(energy_largest_shower[i]);
        hit_summed_adc = std::get<12>(energy_largest_shower[i]);
        hit_integral = std::get<13>(energy_largest_shower[i]);
        hit_multiplicity = std::get<14>(energy_largest_shower[i]);
        hit_charge = std::get<15>(energy_largest_shower[i]);

        energy_LS->Fill();
    }

    //m_file->Write("", TObject::kOverwrite); // save only the new version of the tree - without the arguments was duplicating original tree

  } 




} // namespace



DEFINE_ART_CLASS_TOOL(ShowerRecoTools::myShowerESTAREnergy)
