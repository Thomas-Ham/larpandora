//############################################################################
//### Name:        IShowerTool                                             ###
//### Author:      Dominic Barker                                          ###
//### Date:        13.05.19                                                ###
//### Description: Generic Tool for finding the shower energy. Used in     ###
//###              LArPandoraModularShowerCreation_module.cc               ###
//############################################################################

#ifndef IShowerTool_H
#define IShowerTool_H

//Framwork Includes
#include "art/Framework/Core/ProducesCollector.h"
#include "art/Framework/Principal/Event.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "fhiclcpp/ParameterSet.h"

//LArSoft Includes
#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Algs/LArPandoraShowerAlg.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Algs/ShowerElementHolder.hh"
#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Algs/ShowerProducedPtrsHolder.hh"

namespace ShowerRecoTools {
  class IShowerTool {

  public:
    IShowerTool(const fhicl::ParameterSet& pset)
      : fLArPandoraShowerAlg(pset.get<fhicl::ParameterSet>("LArPandoraShowerAlg"))
      , fRunEventDisplay(pset.get<bool>("EnableEventDisplay")){};

    virtual ~IShowerTool() noexcept = default;

    //Generic Elemnt Finder. Used to calculate thing about the shower.
    virtual int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
                                 art::Event& Event,
                                 reco::shower::ShowerElementHolder& ShowerEleHolder) = 0;

    //Main function that runs the shower tool.  This includes running the derived function
    //that calculates the shower element and also runs the event display if requested
    int
    RunShowerTool(const art::Ptr<recob::PFParticle>& pfparticle,
                  art::Event& Event,
                  reco::shower::ShowerElementHolder& ShowerEleHolder,
                  std::string evd_display_name_append = "")
    {

      int calculation_status = CalculateElement(pfparticle, Event, ShowerEleHolder);
      if (calculation_status != 0) return calculation_status;
      if (fRunEventDisplay) {
        IShowerTool::GetLArPandoraShowerAlg().DebugEVD(
          pfparticle, Event, ShowerEleHolder, evd_display_name_append);
      }
      return calculation_status;
    }

    //Function to initialise the producer i.e produces<std::vector<recob::Vertex> >(); commands go here.
    virtual void
    InitialiseProducers()
    {}

    //Set the point looking back at the producer module show we can make things in the module
    void
    SetPtr(art::ProducesCollector* collector)
    {
      collectorPtr = collector;
    }

    //Initialises the unique ptr holder so that the tool can access it behind the scenes.
    void
    InitaliseProducerPtr(reco::shower::ShowerProducedPtrsHolder& uniqueproducerPtrs)
    {
      UniquePtrs = &uniqueproducerPtrs;
    }

    //End function so the user can add associations
    virtual int
    AddAssociations(const art::Ptr<recob::PFParticle>& pfpPtr,
                    art::Event& Event,
                    reco::shower::ShowerElementHolder& ShowerEleHolder)
    {
      return 0;
    }

  protected:
    const shower::LArPandoraShowerAlg&
    GetLArPandoraShowerAlg() const
    {
      return fLArPandoraShowerAlg;
    };

  private:
    //ptr to the holder of all the unique ptrs.
    reco::shower::ShowerProducedPtrsHolder* UniquePtrs;

    //Algorithm functions
    shower::LArPandoraShowerAlg fLArPandoraShowerAlg;

    //Flags
    bool fRunEventDisplay;

    art::ProducesCollector* collectorPtr;

  protected:
    //Function to return the art:ptr for the corrsponding index iter. This allows the user the make associations
    template <class T>
    art::Ptr<T>
    GetProducedElementPtr(std::string Name,
                          reco::shower::ShowerElementHolder& ShowerEleHolder,
                          int iter = -1)
    {

      //Check the element has been set
      bool check_element = ShowerEleHolder.CheckElement(Name);
      if (!check_element) {
        throw cet::exception("IShowerTool") << "tried to get a element that does not exist. Failed "
                                               "at making the art ptr for Element: "
                                            << Name << std::endl;
      }

      //Check the unique ptr has been set.
      bool check_ptr = UniquePtrs->CheckUniqueProduerPtr(Name);
      if (!check_ptr) {
        throw cet::exception("IShowerTool")
          << "tried to get a ptr that does not exist. Failed at making the art ptr for Element"
          << Name;
      }

      //Check if the user has defined an index if not just use the current shower index/
      int index;
      if (iter != -1) { index = iter; }
      else {
        index = ShowerEleHolder.GetShowerNumber();
      }

      //Make the ptr
      return UniquePtrs->GetArtPtr<T>(Name, index);
    }

    //Function so that the user can add products to the art event. This will set up the unique ptrs and the ptr makers required.
    //Example: InitialiseProduct<std::vector<recob<vertex>>("MyVertex")
    template <class T>
    void
    InitialiseProduct(std::string Name, std::string InstanceName = "")
    {

      if (collectorPtr == nullptr) {
        mf::LogWarning("IShowerTool") << "The art::ProducesCollector ptr has not been set";
        return;
      }

      collectorPtr->produces<T>(InstanceName);
      UniquePtrs->SetShowerUniqueProduerPtr(type<T>(), Name, InstanceName);
    }

    //Function so that the user can add assocations to the event.
    //Example: AddSingle<art::Assn<recob::Vertex,recob::shower>((art::Ptr<recob::Vertex>) Vertex, (art::Prt<recob::shower>) Shower), "myassn")
    template <class T, class A, class B>
    void
    AddSingle(A& a, B& b, std::string Name)
    {
      UniquePtrs->AddSingle<T>(a, b, Name);
    }

    //Function to get the size of the vector, for the event,  that is held in the unique producer ptr that will be put in the event.
    int
    GetVectorPtrSize(std::string Name)
    {
      return UniquePtrs->GetVectorPtrSize(Name);
    }

    void
    PrintPtrs()
    {
      UniquePtrs->PrintPtrs();
    }

    void
    PrintPtr(std::string Name)
    {
      UniquePtrs->PrintPtr(Name);
    }
  };
}

#endif
