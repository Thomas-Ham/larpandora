//############################################################################
//### Name:        ShowerIncrementalTrackHitFinder                         ###
//### Author:      Dom Barker                                              ###
//### Date:        13.05.19                                                ###
//### Description: Tool to incrementally add space points to the initial   ###
//###              track space points until the residuals blow up          ###
//############################################################################

//Framework Includes
#include "art/Utilities/ToolMacros.h"

//LArSoft Includes
#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Tools/IShowerTool.h"

//Root Includes
#include "TGraph2D.h"
#include "TPrincipal.h"

namespace ShowerRecoTools {

  class ShowerIncrementalTrackHitFinder : public IShowerTool {

  public:
    ShowerIncrementalTrackHitFinder(const fhicl::ParameterSet& pset);

    //Generic Direction Finder
    int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
                         art::Event& Event,
                         reco::shower::ShowerElementHolder& ShowerEleHolder) override;

  private:
    std::vector<art::Ptr<recob::SpacePoint>> RunIncrementalSpacePointFinder(
      const art::Event& Event,
      std::vector<art::Ptr<recob::SpacePoint>> const& sps,
      const art::FindManyP<recob::Hit>& fmh);

    void PruneFrontOfSPSPool(std::vector<art::Ptr<recob::SpacePoint>>& sps_pool,
                             std::vector<art::Ptr<recob::SpacePoint>> const& initial_track);

    void PruneTrack(std::vector<art::Ptr<recob::SpacePoint>>& initial_track);

    void AddSpacePointsToSegment(std::vector<art::Ptr<recob::SpacePoint>>& segment,
                                 std::vector<art::Ptr<recob::SpacePoint>>& sps_pool,
                                 size_t num_sps_to_take);

    bool IsSegmentValid(std::vector<art::Ptr<recob::SpacePoint>> const& segment);

    bool IncrementallyFitSegment(const detinfo::DetectorClocksData& clockData,
                                 const detinfo::DetectorPropertiesData& detProp,
                                 std::vector<art::Ptr<recob::SpacePoint>>& segment,
                                 std::vector<art::Ptr<recob::SpacePoint>>& sps_pool,
                                 const art::FindManyP<recob::Hit>& fmh,
                                 double current_residual);

    double FitSegmentAndCalculateResidual(const detinfo::DetectorClocksData& clockData,
                                          const detinfo::DetectorPropertiesData& detProp,
                                          std::vector<art::Ptr<recob::SpacePoint>>& segment,
                                          const art::FindManyP<recob::Hit>& fmh);

    double FitSegmentAndCalculateResidual(const detinfo::DetectorClocksData& clockData,
                                          const detinfo::DetectorPropertiesData& detProp,
                                          std::vector<art::Ptr<recob::SpacePoint>>& segment,
                                          const art::FindManyP<recob::Hit>& fmh,
                                          int& max_residual_point);

    bool RecursivelyReplaceLastSpacePointAndRefit(
      const detinfo::DetectorClocksData& clockData,
      const detinfo::DetectorPropertiesData& detProp,
      std::vector<art::Ptr<recob::SpacePoint>>& segment,
      std::vector<art::Ptr<recob::SpacePoint>>& reduced_sps_pool,
      const art::FindManyP<recob::Hit>& fmh,
      double current_residual);

    bool
    IsResidualOK(double new_residual, double current_residual)
    {
      return new_residual - current_residual < fMaxResidualDiff;
    };
    bool
    IsResidualOK(double new_residual, double current_residual, size_t no_sps)
    {
      return (new_residual - current_residual < fMaxResidualDiff &&
              new_residual / no_sps < fMaxAverageResidual);
    };
    bool
    IsResidualOK(double residual, size_t no_sps)
    {
      return residual / no_sps < fMaxAverageResidual;
    }

    double CalculateResidual(std::vector<art::Ptr<recob::SpacePoint>>& sps,
                             TVector3& PCAEigenvector,
                             TVector3& TrackPosition);

    double CalculateResidual(std::vector<art::Ptr<recob::SpacePoint>>& sps,
                             TVector3& PCAEigenvector,
                             TVector3& TrackPosition,
                             int& max_residual_point);

    //Function to calculate the shower direction using a charge weight 3D PCA calculation.
    TVector3 ShowerPCAVector(std::vector<art::Ptr<recob::SpacePoint>>& sps);

    TVector3 ShowerPCAVector(const detinfo::DetectorClocksData& clockData,
                             const detinfo::DetectorPropertiesData& detProp,
                             const std::vector<art::Ptr<recob::SpacePoint>>& sps,
                             const art::FindManyP<recob::Hit>& fmh);

    std::vector<art::Ptr<recob::SpacePoint>> CreateFakeShowerTrajectory(TVector3 start_position,
                                                                        TVector3 start_direction);
    std::vector<art::Ptr<recob::SpacePoint>> CreateFakeSPLine(TVector3 start_position,
                                                              TVector3 start_direction,
                                                              int npoints);
    void RunTestOfIncrementalSpacePointFinder(const art::Event& Event,
                                              const art::FindManyP<recob::Hit>& dud_fmh);

    void MakeTrackSeed(const detinfo::DetectorClocksData& clockData,
                       const detinfo::DetectorPropertiesData& detProp,
                       std::vector<art::Ptr<recob::SpacePoint>>& segment,
                       const art::FindManyP<recob::Hit>& fmh);

    //Services
    art::InputTag fPFParticleLabel;
    int fVerbose;
    bool fUseShowerDirection;
    bool fChargeWeighted;
    bool fForwardHitsOnly;
    float fMaxResidualDiff;
    float fMaxAverageResidual;
    int fStartFitSize;
    int fNMissPoints;
    float fTrackMaxAdjacentSPDistance;
    bool fRunTest;
    bool fMakeTrackSeed;
    float fStartDistanceCut;
    float fDistanceCut;
    std::string fShowerStartPositionInputLabel;
    std::string fShowerDirectionInputLabel;
    std::string fInitialTrackHitsOutputLabel;
    std::string fInitialTrackSpacePointsOutputLabel;
  };

  ShowerIncrementalTrackHitFinder::ShowerIncrementalTrackHitFinder(const fhicl::ParameterSet& pset)
    : IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools"))
    , fPFParticleLabel(pset.get<art::InputTag>("PFParticleLabel"))
    , fVerbose(pset.get<int>("Verbose"))
    , fUseShowerDirection(pset.get<bool>("UseShowerDirection"))
    , fChargeWeighted(pset.get<bool>("ChargeWeighted"))
    , fForwardHitsOnly(pset.get<bool>("ForwardHitsOnly"))
    , fMaxResidualDiff(pset.get<float>("MaxResidualDiff"))
    , fMaxAverageResidual(pset.get<float>("MaxAverageResidual"))
    , fStartFitSize(pset.get<int>("StartFitSize"))
    , fNMissPoints(pset.get<int>("NMissPoints"))
    , fTrackMaxAdjacentSPDistance(pset.get<float>("TrackMaxAdjacentSPDistance"))
    , fRunTest(0)
    , fMakeTrackSeed(pset.get<bool>("MakeTrackSeed"))
    , fStartDistanceCut(pset.get<float>("StartDistanceCut"))
    , fDistanceCut(pset.get<float>("DistanceCut"))
    , fShowerStartPositionInputLabel(pset.get<std::string>("ShowerStartPositionInputLabel"))
    , fShowerDirectionInputLabel(pset.get<std::string>("ShowerDirectionInputLabel"))
    ,

    fInitialTrackHitsOutputLabel(pset.get<std::string>("InitialTrackHitsOutputLabel"))
    , fInitialTrackSpacePointsOutputLabel(
        pset.get<std::string>("InitialTrackSpacePointsOutputLabel"))
  {
    if (fStartFitSize == 0) {
      throw cet::exception("ShowerIncrementalTrackHitFinder")
        << "We cannot make a track if you don't gives us at leats one hit. Change fStartFitSize "
           "please to something sensible";
    }
  }

  int
  ShowerIncrementalTrackHitFinder::CalculateElement(
    const art::Ptr<recob::PFParticle>& pfparticle,
    art::Event& Event,
    reco::shower::ShowerElementHolder& ShowerEleHolder)
  {

    //This is all based on the shower vertex being known. If it is not lets not do the track
    if (!ShowerEleHolder.CheckElement(fShowerStartPositionInputLabel)) {
      if (fVerbose)
        mf::LogError("ShowerIncrementalTrackHitFinder")
          << "Start position not set, returning " << std::endl;
      return 1;
    }

    // Get the assocated pfParicle Handle
    auto const pfpHandle = Event.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleLabel);

    // Get the spacepoint - PFParticle assn
    const art::FindManyP<recob::SpacePoint>& fmspp =
      ShowerEleHolder.GetFindManyP<recob::SpacePoint>(pfpHandle, Event, fPFParticleLabel);

    // Get the spacepoints
    auto const spHandle = Event.getValidHandle<std::vector<recob::SpacePoint>>(fPFParticleLabel);

    // Get the hits associated with the space points
    const art::FindManyP<recob::Hit>& fmh =
      ShowerEleHolder.GetFindManyP<recob::Hit>(spHandle, Event, fPFParticleLabel);

    // Get the SpacePoints
    std::vector<art::Ptr<recob::SpacePoint>> spacePoints = fmspp.at(pfparticle.key());

    //We cannot progress with no spacepoints.
    if (spacePoints.empty()) {
      if (fVerbose)
        mf::LogError("ShowerIncrementalTrackHitFinder")
          << "No space points, returning " << std::endl;
      return 1;
    }

    TVector3 ShowerStartPosition = {-999, -999, -999};
    ShowerEleHolder.GetElement(fShowerStartPositionInputLabel, ShowerStartPosition);

    //Decide if the you want to use the direction of the shower or make one.
    if (fUseShowerDirection) {

      if (!ShowerEleHolder.CheckElement(fShowerDirectionInputLabel)) {
        if (fVerbose)
          mf::LogError("ShowerIncrementalTrackHitFinder")
            << "Direction not set, returning " << std::endl;
        return 1;
      }

      TVector3 ShowerDirection = {-999, -999, -999};
      ShowerEleHolder.GetElement(fShowerDirectionInputLabel, ShowerDirection);

      //Order the spacepoints
      IShowerTool::GetLArPandoraShowerAlg().OrderShowerSpacePoints(
        spacePoints, ShowerStartPosition, ShowerDirection);
      //Remove the back hits if requird.
      if (fForwardHitsOnly) {
        int back_sps = 0;
        for (auto spacePoint : spacePoints) {
          double proj = IShowerTool::GetLArPandoraShowerAlg().SpacePointProjection(
            spacePoint, ShowerStartPosition, ShowerDirection);
          if (proj < 0) { ++back_sps; }
          if (proj > 0) { break; }
        }
        spacePoints.erase(spacePoints.begin(), spacePoints.begin() + back_sps);
      }
    }
    else {
      //Order the spacepoint using the magnitude away from the vertex
      IShowerTool::GetLArPandoraShowerAlg().OrderShowerSpacePoints(spacePoints,
                                                                   ShowerStartPosition);
    }

    //Remove the first x spacepoints
    int frontsp = 0;
    for (auto spacePoint : spacePoints) {
      double dist =
        (IShowerTool::GetLArPandoraShowerAlg().SpacePointPosition(spacePoint) - ShowerStartPosition)
          .Mag();
      if (dist > fStartDistanceCut) { break; }
      ++frontsp;
    }
    spacePoints.erase(spacePoints.begin(), spacePoints.begin() + frontsp);

    //Bin anything above x cm
    int sp_iter = 0;
    for (auto spacePoint : spacePoints) {
      double dist =
        (IShowerTool::GetLArPandoraShowerAlg().SpacePointPosition(spacePoint) - ShowerStartPosition)
          .Mag();
      if (dist > fDistanceCut) { break; }
      ++sp_iter;
    }
    spacePoints.erase(spacePoints.begin() + sp_iter, spacePoints.end());

    if (spacePoints.size() < 3) {
      if (fVerbose)
        mf::LogError("ShowerIncrementalTrackHitFinder")
          << "Not enough spacepoints bailing" << std::endl;
      return 1;
    }

    //Create fake hits and test the algorithm
    if (fRunTest) RunTestOfIncrementalSpacePointFinder(Event, fmh);

    //Actually runt he algorithm.
    std::vector<art::Ptr<recob::SpacePoint>> track_sps =
      RunIncrementalSpacePointFinder(Event, spacePoints, fmh);

    // Get the hits associated to the space points and seperate them by planes
    std::vector<art::Ptr<recob::Hit>> trackHits;
    for (auto const& spacePoint : track_sps) {
      std::vector<art::Ptr<recob::Hit>> hits = fmh.at(spacePoint.key());
      for (auto const& hit : hits) {
        trackHits.push_back(hit);
      }
    }

    //Add to the holder
    ShowerEleHolder.SetElement(trackHits, fInitialTrackHitsOutputLabel);
    ShowerEleHolder.SetElement(track_sps, fInitialTrackSpacePointsOutputLabel);

    return 0;
  }

  TVector3
  ShowerIncrementalTrackHitFinder::ShowerPCAVector(std::vector<art::Ptr<recob::SpacePoint>>& sps)
  {

    //Initialise the the PCA.
    TPrincipal* pca = new TPrincipal(3, "");

    //Normalise the spacepoints, charge weight and add to the PCA.
    for (auto& sp : sps) {

      TVector3 sp_position = IShowerTool::GetLArPandoraShowerAlg().SpacePointPosition(sp);

      double sp_coord[3];
      sp_coord[0] = sp_position.X();
      sp_coord[1] = sp_position.Y();
      sp_coord[2] = sp_position.Z();

      //Add to the PCA
      pca->AddRow(sp_coord);
    }

    //Evaluate the PCA
    pca->MakePrincipals();

    //Get the Eigenvectors.
    const TMatrixD* Eigenvectors = pca->GetEigenVectors();

    TVector3 Eigenvector = {(*Eigenvectors)[0][0], (*Eigenvectors)[1][0], (*Eigenvectors)[2][0]};

    delete pca;

    return Eigenvector;
  }

  //Function to calculate the shower direction using a charge weight 3D PCA calculation.
  TVector3
  ShowerIncrementalTrackHitFinder::ShowerPCAVector(
    const detinfo::DetectorClocksData& clockData,
    const detinfo::DetectorPropertiesData& detProp,
    const std::vector<art::Ptr<recob::SpacePoint>>& sps,
    const art::FindManyP<recob::Hit>& fmh)
  {

    //Initialise the the PCA.
    TPrincipal* pca = new TPrincipal(3, "");

    float TotalCharge = 0;

    //Normalise the spacepoints, charge weight and add to the PCA.
    for (auto& sp : sps) {

      TVector3 sp_position = IShowerTool::GetLArPandoraShowerAlg().SpacePointPosition(sp);

      float wht = 1;

      if (fChargeWeighted) {

        //Get the charge.
        float Charge = IShowerTool::GetLArPandoraShowerAlg().SpacePointCharge(sp, fmh);
        //        std::cout << "Charge: " << Charge << std::endl;

        //Get the time of the spacepoint
        float Time = IShowerTool::GetLArPandoraShowerAlg().SpacePointTime(sp, fmh);

        //Correct for the lifetime at the moment.
        Charge *= std::exp((sampling_rate(clockData) * Time) / (detProp.ElectronLifetime() * 1e3));
        //        std::cout << "Charge: "<< Charge << std::endl;

        //Charge Weight
        wht *= std::sqrt(Charge / TotalCharge);
      }

      double sp_coord[3];
      sp_coord[0] = sp_position.X() * wht;
      sp_coord[1] = sp_position.Y() * wht;
      sp_coord[2] = sp_position.Z() * wht;

      //Add to the PCA
      pca->AddRow(sp_coord);
    }

    //Evaluate the PCA
    pca->MakePrincipals();

    //Get the Eigenvectors.
    const TMatrixD* Eigenvectors = pca->GetEigenVectors();

    TVector3 Eigenvector = {(*Eigenvectors)[0][0], (*Eigenvectors)[1][0], (*Eigenvectors)[2][0]};

    delete pca;

    return Eigenvector;
  }

  //Function to remove the spacepoint with the highest residual until we have a track which matches the
  //residual criteria.
  void
  ShowerIncrementalTrackHitFinder::MakeTrackSeed(const detinfo::DetectorClocksData& clockData,
                                                 const detinfo::DetectorPropertiesData& detProp,
                                                 std::vector<art::Ptr<recob::SpacePoint>>& segment,
                                                 const art::FindManyP<recob::Hit>& fmh)
  {

    bool ok = true;

    int maxresidual_point = 0;

    //Check the residual
    double residual =
      FitSegmentAndCalculateResidual(clockData, detProp, segment, fmh, maxresidual_point);

    //Is it okay
    ok = IsResidualOK(residual, segment.size());

    //Remove points until we can fit a track.
    while (!ok && segment.size() != 1) {

      //Remove the point with the highest residual
      for (auto sp = segment.begin(); sp != segment.end(); ++sp) {
        if (sp->key() == (unsigned)maxresidual_point) {
          segment.erase(sp);
          break;
        }
      }

      //Check the residual
      double residual =
        FitSegmentAndCalculateResidual(clockData, detProp, segment, fmh, maxresidual_point);

      //Is it okay
      ok = IsResidualOK(residual, segment.size());
    }
  }

  std::vector<art::Ptr<recob::SpacePoint>>
  ShowerIncrementalTrackHitFinder::RunIncrementalSpacePointFinder(
    const art::Event& Event,
    std::vector<art::Ptr<recob::SpacePoint>> const& sps,
    const art::FindManyP<recob::Hit>& fmh)
  {

    auto const clockData =
      art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(Event);
    auto const detProp =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(Event, clockData);

    //Create space point pool (yes we are copying the input vector because we're going to twiddle with it
    std::vector<art::Ptr<recob::SpacePoint>> sps_pool = sps;
    std::vector<art::Ptr<recob::SpacePoint>> initial_track;
    std::vector<art::Ptr<recob::SpacePoint>> track_segment_copy;

    while (sps_pool.size() > 0) {
      //PruneFrontOfSPSPool(sps_pool, initial_track);

      std::vector<art::Ptr<recob::SpacePoint>> track_segment;
      AddSpacePointsToSegment(track_segment, sps_pool, (size_t)(fStartFitSize));
      if (!IsSegmentValid(track_segment)) {
        //Clear the pool and lets leave this place
        sps_pool.clear();
        break;
      }

      //Lets really try to make the initial track seed.
      if (fMakeTrackSeed && sps_pool.size() + fStartFitSize == sps.size()) {
        MakeTrackSeed(clockData, detProp, track_segment, fmh);
        if (track_segment.empty()) break;

        track_segment_copy = track_segment;
      }

      //A sleight of hand coming up.  We are going to move the last sp from the segment back into the pool so
      //that it makes kick starting the recursion easier (sneaky)
      //TODO defend against segments that are too small for this to work (I dunno who is running the alg with
      //fStartFitMinSize==0 but whatever
      sps_pool.insert(sps_pool.begin(), track_segment.back());
      track_segment.pop_back();
      double current_residual = 0;
      size_t initial_segment_size = track_segment.size();

      IncrementallyFitSegment(clockData, detProp, track_segment, sps_pool, fmh, current_residual);

      //Check if the track has grown in size at all
      if (initial_segment_size == track_segment.size()) {
        //The incremental fitter could not grow th track at all.  SAD!
        //Clear the pool and let's get out of here
        sps_pool.clear();
        break;
      }
      else {
        //We did some good fitting and everyone is really happy with it
        //Let's store all of the hits in the final space point vector
        AddSpacePointsToSegment(initial_track, track_segment, track_segment.size());
      }
    }

    //If we have failed then no worry we have the seed. We shall just give those points.
    if (fMakeTrackSeed && initial_track.empty()) initial_track = track_segment_copy;

    //Runt the algorithm that attepmts to remove hits too far away from the track.
    PruneTrack(initial_track);

    return initial_track;
  }

  void
  ShowerIncrementalTrackHitFinder::PruneFrontOfSPSPool(
    std::vector<art::Ptr<recob::SpacePoint>>& sps_pool,
    std::vector<art::Ptr<recob::SpacePoint>> const& initial_track)
  {

    //If the initial track is empty then there is no pruning to do
    if (initial_track.empty()) return;
    double distance = IShowerTool::GetLArPandoraShowerAlg().DistanceBetweenSpacePoints(
      initial_track.back(), sps_pool.front());
    while (distance > 1 && sps_pool.size() > 0) {
      sps_pool.erase(sps_pool.begin());
      distance = IShowerTool::GetLArPandoraShowerAlg().DistanceBetweenSpacePoints(
        initial_track.back(), sps_pool.front());
    }
    return;
  }

  void
  ShowerIncrementalTrackHitFinder::PruneTrack(
    std::vector<art::Ptr<recob::SpacePoint>>& initial_track)
  {

    if (initial_track.empty()) return;
    std::vector<art::Ptr<recob::SpacePoint>>::iterator sps_it = initial_track.begin();
    while (sps_it != std::next(initial_track.end(), -1)) {
      std::vector<art::Ptr<recob::SpacePoint>>::iterator next_sps_it = std::next(sps_it, 1);
      double distance =
        IShowerTool::GetLArPandoraShowerAlg().DistanceBetweenSpacePoints(*sps_it, *next_sps_it);
      if (distance > fTrackMaxAdjacentSPDistance) { initial_track.erase(next_sps_it); }
      else {
        sps_it++;
      }
    }
    return;
  }

  void
  ShowerIncrementalTrackHitFinder::AddSpacePointsToSegment(
    std::vector<art::Ptr<recob::SpacePoint>>& segment,
    std::vector<art::Ptr<recob::SpacePoint>>& sps_pool,
    size_t num_sps_to_take)
  {
    size_t new_segment_size = segment.size() + num_sps_to_take;
    while (segment.size() < new_segment_size && sps_pool.size() > 0) {
      segment.push_back(sps_pool[0]);
      sps_pool.erase(sps_pool.begin());
    }
    return;
  }

  bool
  ShowerIncrementalTrackHitFinder::IsSegmentValid(
    std::vector<art::Ptr<recob::SpacePoint>> const& segment)
  {
    bool ok = true;
    if (segment.size() < (size_t)(fStartFitSize)) return !ok;

    return ok;
  }

  bool
  ShowerIncrementalTrackHitFinder::IncrementallyFitSegment(
    const detinfo::DetectorClocksData& clockData,
    const detinfo::DetectorPropertiesData& detProp,
    std::vector<art::Ptr<recob::SpacePoint>>& segment,
    std::vector<art::Ptr<recob::SpacePoint>>& sps_pool,
    const art::FindManyP<recob::Hit>& fmh,
    double current_residual)
  {

    bool ok = true;
    //Firstly, are there any space points left???
    if (sps_pool.empty()) return !ok;
    //Fit the current line
    current_residual = FitSegmentAndCalculateResidual(clockData, detProp, segment, fmh);
    //Take a space point from the pool and plonk it onto the seggieweggie
    AddSpacePointsToSegment(segment, sps_pool, 1);
    //Fit again
    double residual = FitSegmentAndCalculateResidual(clockData, detProp, segment, fmh);

    ok = IsResidualOK(residual, current_residual, segment.size());
    if (!ok) {
      //Create a sub pool of space points to pass to the refitter
      std::vector<art::Ptr<recob::SpacePoint>> sub_sps_pool;
      AddSpacePointsToSegment(sub_sps_pool, sps_pool, fNMissPoints);
      //We'll need an additional copy of this pool, as we will need the space points if we have to start a new
      //segment later, but all of the funtionality drains the pools during use
      std::vector<art::Ptr<recob::SpacePoint>> sub_sps_pool_cache = sub_sps_pool;
      //The most recently added SP to the segment is bad but it will get thrown away by RecursivelyReplaceLastSpacePointAndRefit
      //It's possible that we will need it if we end up forming an entirely new line from scratch, so
      //add the bad SP to the front of the cache
      sub_sps_pool_cache.insert(sub_sps_pool_cache.begin(), segment.back());
      ok = RecursivelyReplaceLastSpacePointAndRefit(
        clockData, detProp, segment, sub_sps_pool, fmh, current_residual);
      if (ok) {
        //The refitting may have dropped a couple of points but it managed to find a point that kept the residual
        //at a sensible value.
        //Add the remaining SPS in the reduced pool back t othe start of the larger pool
        while (sub_sps_pool.size() > 0) {
          sps_pool.insert(sps_pool.begin(), sub_sps_pool.back());
          sub_sps_pool.pop_back();
        }
        //We'll need the latest residual now that we've managed to refit the track
        residual = FitSegmentAndCalculateResidual(clockData, detProp, segment, fmh);
      }
      else {
        //All of the space points in the reduced pool could not sensibly refit the track.  The reduced pool will be
        //empty so move all of the cached space points back into the main pool
        //        std::cout<<"The refitting was NOT a success, dumping  all " << sub_sps_pool_cache.size() << "  sps back into the pool" << std::endl;
        while (sub_sps_pool_cache.size() > 0) {
          sps_pool.insert(sps_pool.begin(), sub_sps_pool_cache.back());
          sub_sps_pool_cache.pop_back();
        }
        //The bad point is still on the segment, so remove it
        segment.pop_back();
        return !ok;
      }
    }

    //Update the residual
    current_residual = residual;

    //Round and round we go
    //NOBODY GETS OFF MR BONES WILD RIDE
    return IncrementallyFitSegment(clockData, detProp, segment, sps_pool, fmh, current_residual);
  }

  double
  ShowerIncrementalTrackHitFinder::FitSegmentAndCalculateResidual(
    const detinfo::DetectorClocksData& clockData,
    const detinfo::DetectorPropertiesData& detProp,
    std::vector<art::Ptr<recob::SpacePoint>>& segment,
    const art::FindManyP<recob::Hit>& fmh)
  {

    TVector3 primary_axis;
    if (fChargeWeighted)
      primary_axis = ShowerPCAVector(clockData, detProp, segment, fmh);
    else
      primary_axis = ShowerPCAVector(segment);

    TVector3 segment_centre;
    if (fChargeWeighted)
      segment_centre =
        IShowerTool::GetLArPandoraShowerAlg().ShowerCentre(clockData, detProp, segment, fmh);
    else
      segment_centre = IShowerTool::GetLArPandoraShowerAlg().ShowerCentre(segment);

    double residual = CalculateResidual(segment, primary_axis, segment_centre);

    return residual;
  }

  double
  ShowerIncrementalTrackHitFinder::FitSegmentAndCalculateResidual(
    const detinfo::DetectorClocksData& clockData,
    const detinfo::DetectorPropertiesData& detProp,
    std::vector<art::Ptr<recob::SpacePoint>>& segment,
    const art::FindManyP<recob::Hit>& fmh,
    int& max_residual_point)
  {

    TVector3 primary_axis;
    if (fChargeWeighted)
      primary_axis = ShowerPCAVector(clockData, detProp, segment, fmh);
    else
      primary_axis = ShowerPCAVector(segment);

    TVector3 segment_centre;
    if (fChargeWeighted)
      segment_centre =
        IShowerTool::GetLArPandoraShowerAlg().ShowerCentre(clockData, detProp, segment, fmh);
    else
      segment_centre = IShowerTool::GetLArPandoraShowerAlg().ShowerCentre(segment);

    double residual = CalculateResidual(segment, primary_axis, segment_centre, max_residual_point);

    return residual;
  }

  bool
  ShowerIncrementalTrackHitFinder::RecursivelyReplaceLastSpacePointAndRefit(
    const detinfo::DetectorClocksData& clockData,
    const detinfo::DetectorPropertiesData& detProp,
    std::vector<art::Ptr<recob::SpacePoint>>& segment,
    std::vector<art::Ptr<recob::SpacePoint>>& reduced_sps_pool,
    const art::FindManyP<recob::Hit>& fmh,
    double current_residual)
  {

    bool ok = true;
    //If the pool is empty, then there is nothing to do (sad)
    if (reduced_sps_pool.empty()) return !ok;
    //Drop the last space point
    segment.pop_back();
    //Add one point
    AddSpacePointsToSegment(segment, reduced_sps_pool, 1);
    double residual = FitSegmentAndCalculateResidual(clockData, detProp, segment, fmh);

    ok = IsResidualOK(residual, current_residual, segment.size());
    //    std::cout<<"recursive refit: isok " << ok << "  res: " << residual << "  curr res: " << current_residual << std::endl;
    if (ok) return ok;
    return RecursivelyReplaceLastSpacePointAndRefit(
      clockData, detProp, segment, reduced_sps_pool, fmh, current_residual);
  }

  double
  ShowerIncrementalTrackHitFinder::CalculateResidual(std::vector<art::Ptr<recob::SpacePoint>>& sps,
                                                     TVector3& PCAEigenvector,
                                                     TVector3& TrackPosition)
  {

    double Residual = 0;

    for (auto const& sp : sps) {

      //Get the relative position of the spacepoint
      TVector3 pos = IShowerTool::GetLArPandoraShowerAlg().SpacePointPosition(sp) - TrackPosition;

      //Gen the perpendicular distance
      double len = pos.Dot(PCAEigenvector);
      double perp = (pos - len * PCAEigenvector).Mag();

      Residual += perp;
    }
    return Residual;
  }

  double
  ShowerIncrementalTrackHitFinder::CalculateResidual(std::vector<art::Ptr<recob::SpacePoint>>& sps,
                                                     TVector3& PCAEigenvector,
                                                     TVector3& TrackPosition,
                                                     int& max_residual_point)
  {

    double Residual = 0;
    double max_residual = -999;

    for (auto const& sp : sps) {

      //Get the relative position of the spacepoint
      TVector3 pos = IShowerTool::GetLArPandoraShowerAlg().SpacePointPosition(sp) - TrackPosition;

      //Gen the perpendicular distance
      double len = pos.Dot(PCAEigenvector);
      double perp = (pos - len * PCAEigenvector).Mag();

      Residual += perp;

      if (perp > max_residual) {
        max_residual = perp;
        max_residual_point = sp.key();
      }
    }
    return Residual;
  }

  std::vector<art::Ptr<recob::SpacePoint>>
  ShowerIncrementalTrackHitFinder::CreateFakeShowerTrajectory(TVector3 start_position,
                                                              TVector3 start_direction)
  {
    std::vector<art::Ptr<recob::SpacePoint>> fake_sps;
    std::vector<art::Ptr<recob::SpacePoint>> segment_a =
      CreateFakeSPLine(start_position, start_direction, 20);
    fake_sps.insert(std::end(fake_sps), std::begin(segment_a), std::end(segment_a));

    //make a new segment:
    TVector3 sp_position =
      IShowerTool::GetLArPandoraShowerAlg().SpacePointPosition(fake_sps.back());
    TVector3 direction = start_direction;
    direction.RotateX(10. * 3.142 / 180.);
    std::vector<art::Ptr<recob::SpacePoint>> segment_b =
      CreateFakeSPLine(sp_position, direction, 10);
    fake_sps.insert(std::end(fake_sps), std::begin(segment_b), std::end(segment_b));

    //Now make three branches that come from the end of the segment
    TVector3 branching_position =
      IShowerTool::GetLArPandoraShowerAlg().SpacePointPosition(fake_sps.back());

    TVector3 direction_branch_a = direction;
    direction_branch_a.RotateZ(15. * 3.142 / 180.);
    std::vector<art::Ptr<recob::SpacePoint>> branch_a =
      CreateFakeSPLine(branching_position, direction_branch_a, 6);
    fake_sps.insert(std::end(fake_sps), std::begin(branch_a), std::end(branch_a));

    TVector3 direction_branch_b = direction;
    direction_branch_b.RotateY(20. * 3.142 / 180.);
    std::vector<art::Ptr<recob::SpacePoint>> branch_b =
      CreateFakeSPLine(branching_position, direction_branch_b, 10);
    fake_sps.insert(std::end(fake_sps), std::begin(branch_b), std::end(branch_b));

    TVector3 direction_branch_c = direction;
    direction_branch_c.RotateX(3. * 3.142 / 180.);
    std::vector<art::Ptr<recob::SpacePoint>> branch_c =
      CreateFakeSPLine(branching_position, direction_branch_c, 20);
    fake_sps.insert(std::end(fake_sps), std::begin(branch_c), std::end(branch_c));

    return fake_sps;
  }

  std::vector<art::Ptr<recob::SpacePoint>>
  ShowerIncrementalTrackHitFinder::CreateFakeSPLine(TVector3 start_position,
                                                    TVector3 start_direction,
                                                    int npoints)
  {
    std::vector<art::Ptr<recob::SpacePoint>> fake_sps;
    art::ProductID prod_id(std::string("totally_genuine"));
    size_t current_id = 500000;

    double step_length = 0.2;
    for (double i_point = 0; i_point < npoints; i_point++) {
      TVector3 new_position = start_position + i_point * step_length * start_direction;
      Double32_t xyz[3] = {new_position.X(), new_position.Y(), new_position.Z()};
      Double32_t err[3] = {0., 0., 0.};
      recob::SpacePoint* sp = new recob::SpacePoint(xyz, err, 0, 1);
      fake_sps.emplace_back(art::Ptr<recob::SpacePoint>(prod_id, sp, current_id++));
    }
    return fake_sps;
  }

  void
  ShowerIncrementalTrackHitFinder::RunTestOfIncrementalSpacePointFinder(
    const art::Event& Event,
    const art::FindManyP<recob::Hit>& dud_fmh)
  {
    TVector3 start_position(50, 50, 50);
    TVector3 start_direction(0, 0, 1);
    std::vector<art::Ptr<recob::SpacePoint>> fake_sps =
      CreateFakeShowerTrajectory(start_position, start_direction);

    IShowerTool::GetLArPandoraShowerAlg().OrderShowerSpacePoints(fake_sps, start_position);

    std::vector<art::Ptr<recob::SpacePoint>> track_sps =
      RunIncrementalSpacePointFinder(Event, fake_sps, dud_fmh);

    TGraph2D graph_sps;
    for (size_t i_sp = 0; i_sp < fake_sps.size(); i_sp++) {
      graph_sps.SetPoint(graph_sps.GetN(),
                         fake_sps[i_sp]->XYZ()[0],
                         fake_sps[i_sp]->XYZ()[1],
                         fake_sps[i_sp]->XYZ()[2]);
    }
    TGraph2D graph_track_sps;
    for (size_t i_sp = 0; i_sp < track_sps.size(); i_sp++) {
      graph_track_sps.SetPoint(graph_track_sps.GetN(),
                               track_sps[i_sp]->XYZ()[0],
                               track_sps[i_sp]->XYZ()[1],
                               track_sps[i_sp]->XYZ()[2]);
    }

    art::ServiceHandle<art::TFileService> tfs;

    TCanvas* canvas = tfs->make<TCanvas>("test_inc_can", "test_inc_can");
    canvas->SetName("test_inc_can");
    graph_sps.SetMarkerStyle(8);
    graph_sps.SetMarkerColor(1);
    graph_sps.SetFillColor(1);
    graph_sps.Draw("p");

    graph_track_sps.SetMarkerStyle(8);
    graph_track_sps.SetMarkerColor(2);
    graph_track_sps.SetFillColor(2);
    graph_track_sps.Draw("samep");
    canvas->Write();

    fRunTest = false;
    return;
  }
}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerIncrementalTrackHitFinder)
