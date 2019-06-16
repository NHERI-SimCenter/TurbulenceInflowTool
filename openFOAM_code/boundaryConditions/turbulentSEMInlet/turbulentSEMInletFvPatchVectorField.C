/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "turbulentSEMInletFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "momentOfInertia.H"
#include "globalIndex.H"
#include "IFstream.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::label Foam::turbulentSEMInletFvPatchVectorField::seedIterMax_ = 1000;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::pointToPointPlanarInterpolation&
Foam::turbulentSEMInletFvPatchVectorField::patchMapper() const
{
    // Initialise interpolation (2D planar interpolation by triangulation)
    if (mapperPtr_.empty())
    {

        // Reread values and interpolate
        fileName samplePointsFile
        (
            this->db().time().path()
           /this->db().time().caseConstant()
           /"boundaryData"
           /this->patch().name()
           /"points"
        );

        pointField samplePoints((IFstream(samplePointsFile)()));

        if (debug)
        {
            InfoInFunction
                << " Read " << samplePoints.size() << " sample points from "
                << samplePointsFile << endl;
        }


        // tbd: run-time selection
        bool nearestOnly =
        (
           !mapMethod_.empty()
         && mapMethod_ != "planarInterpolation"
        );

        // Allocate the interpolator
        mapperPtr_.reset
        (
            new pointToPointPlanarInterpolation
            (
                samplePoints,
                this->patch().patch().faceCentres(),
                perturb_,
                nearestOnly
            )
        );
    }

    return *mapperPtr_;
}


void Foam::turbulentSEMInletFvPatchVectorField::initialisePatch()
{
    const vectorField nf(patch().nf());

    // Patch normal points into domain
    patchNormal_ = -gAverage(nf);

    // Check that patch is planar
    scalar error = max(magSqr(patchNormal_ + nf));

    if (error > SMALL)
    {
        WarningInFunction
            << "Patch " << patch().name() << " is not planar"
            << endl;
    }

    patchNormal_ /= mag(patchNormal_) + ROOTVSMALL;

    // Decompose the patch faces into triangles to enable point search

    const polyPatch& patch = this->patch().patch();
    const pointField& points = patch.points();

    // Triangulate the patch faces and create addressing
    DynamicList<label> triToFace(2*patch.size());
    DynamicList<scalar> triMagSf(2*patch.size());
    DynamicList<face> triFace(2*patch.size());
    DynamicList<face> tris(5);

    // Set zero value at the start of the tri area list
    triMagSf.append(0.0);

    forAll(patch, faceI)
    {
        const face& f = patch[faceI];

        tris.clear();
        f.triangles(points, tris);

        forAll(tris, i)
        {
            triToFace.append(faceI);
            triFace.append(tris[i]);
            triMagSf.append(tris[i].mag(points));
        }
    }

    forAll(sumTriMagSf_, i)
    {
        sumTriMagSf_[i] = 0.0;
    }

    sumTriMagSf_[Pstream::myProcNo() + 1] = sum(triMagSf);

    Pstream::listCombineGather(sumTriMagSf_, maxEqOp<scalar>());
    Pstream::listCombineScatter(sumTriMagSf_);

    for (label i = 1; i < triMagSf.size(); i++)
    {
        triMagSf[i] += triMagSf[i-1];
    }

    // Transfer to persistent storage
    triFace_.transfer(triFace);
    triToFace_.transfer(triToFace);
    triCumulativeMagSf_.transfer(triMagSf);

    // Convert sumTriMagSf_ into cumulative sum of areas per proc
    for (label i = 1; i < sumTriMagSf_.size(); i++)
    {
        sumTriMagSf_[i] += sumTriMagSf_[i-1];
    }

    // Global patch area (over all processors)
    patchArea_ = sumTriMagSf_.last();

    // Local patch bounds (this processor)
    patchBounds_ = boundBox(patch.localPoints(), false);
    patchBounds_.inflate(0.1);

    // Determine if all eddies spawned from a single processor
    singleProc_ = patch.size() == returnReduce(patch.size(), sumOp<label>());
    reduce(singleProc_, orOp<bool>());
}


void Foam::turbulentSEMInletFvPatchVectorField::initialiseEddyBox()
{
    const scalarField& magSf = patch().magSf();

    // Maximum extent across all processors
    maxSigmaX_ = cmptMax(vector(fabs(gMax(Lu_)&patchNormal_),fabs(gMax(Lv_)&patchNormal_),fabs(gMax(Lw_)&patchNormal_)));

    // Eddy box volume
    V0_ = 2*gSum(magSf)*maxSigmaX_;

    if (Pstream::master())
    {
        Info<< "Patch: " << patch().patch().name() << " eddy box:" << nl
            << "    volume    : " << V0_ << nl
            << "    maxSigmaX : " << maxSigmaX_ << nl
            << endl;
    }
}


Foam::pointIndexHit Foam::turbulentSEMInletFvPatchVectorField::setNewPosition
(
    const bool global
)
{
    // Initialise to miss
    pointIndexHit pos(false, vector::max, -1);

    const polyPatch& patch = this->patch().patch();
    const pointField& points = patch.points();

    if (global)
    {
        scalar areaFraction = patchArea_*rndGen_.globalScalar01();

        // Determine which processor to use
        label procI = 0;

        forAllReverse(sumTriMagSf_, i)
        {
            if (areaFraction >= sumTriMagSf_[i])
            {
                procI = i;
                break;
            }
        }

        if (Pstream::myProcNo() == procI)
        {
            // Find corresponding decomposed face triangle
            label triI = 0;
            scalar offset = sumTriMagSf_[procI];
            forAllReverse(triCumulativeMagSf_, i)
            {
                if (areaFraction > triCumulativeMagSf_[i] + offset)
                {
                    triI = i;
                    break;
                }
            }

            // Find random point in triangle
            const face& tf = triFace_[triI];
            const triPointRef tri(points[tf[0]], points[tf[1]], points[tf[2]]);

            pos.setHit();
            pos.setIndex(triToFace_[triI]);
            pos.rawPoint() = tri.randomPoint(rndGen_);
        }
    }
    else
    {
        // Find corresponding decomposed face triangle on local processor
        label triI = 0;
        scalar maxAreaLimit = triCumulativeMagSf_.last();
        scalar areaFraction = maxAreaLimit*rndGen_.scalar01();

        forAllReverse(triCumulativeMagSf_, i)
        {
            if (areaFraction > triCumulativeMagSf_[i])
            {
                triI = i;
                break;
            }
        }

        // Find random point in triangle
        const face& tf = triFace_[triI];
        const triPointRef tri(points[tf[0]], points[tf[1]], points[tf[2]]);

        pos.setHit();
        pos.setIndex(triToFace_[triI]);
        pos.rawPoint() = tri.randomPoint(rndGen_);
    }

    return pos;
}


void Foam::turbulentSEMInletFvPatchVectorField::initialiseEddies()
{
    DynamicList<eddy> eddies(size());

    // Initialise eddy properties
    scalar sumVolEddy = 0;
    scalar sumVolEddyAllProc = 0;
    scalar sigmaRatio = sigmaToLengthscaleRatio();

    while (sumVolEddyAllProc/V0_ < eddyDensity_)
    {
        // Get new parallel consistent position
        pointIndexHit pos(setNewPosition(true));
        label faceI = pos.index();

        // Note: only 1 processor will pick up this face
        if (faceI != -1)
        {
            eddy e
            (
                velocityShape_,
                faceI,
                pos.hitPoint(),
                rndGen_.scalarAB(-maxSigmaX_, maxSigmaX_),
                U_[faceI]&patchNormal_,
                sigmaRatio*Lu_[faceI],
                sigmaRatio*Lv_[faceI],
                sigmaRatio*Lw_[faceI],
                R_[faceI],
                rndGen_
            );

            eddies.append(e);
            sumVolEddy += e.volume();
        }

        sumVolEddyAllProc = returnReduce(sumVolEddy, sumOp<scalar>());
    }

    eddies_.transfer(eddies);
    nEddy_ = eddies_.size();

    if (debug)
    {
        Pout<< "Patch:" << patch().patch().name();

        if (Pstream::parRun())
        {
            Pout<< " processor:" << Pstream::myProcNo();
        }

        Pout<< " seeded:" << nEddy_ << " eddies" << endl;
    }

    reduce(nEddy_, sumOp<label>());

    if (nEddy_ > 0)
    {
        Info<< "Turbulent SEM patch: " << patch().name()
            << " seeded " << nEddy_ << " eddies with total volume "
            << sumVolEddyAllProc
            << endl;
    }
    else
    {
        WarningInFunction
            << "Patch: " << patch().patch().name()
            << " on field " << internalField().name()
            << ": No eddies seeded - please check your set-up" << endl;
    }
}


void Foam::turbulentSEMInletFvPatchVectorField::convectEddies
(
    const scalar deltaT
)
{
    // Note: all operations applied to local processor only

    label nRecycled = 0;
    scalar sigmaRatio = sigmaToLengthscaleRatio();

    forAll(eddies_, eddyI)
    {
        eddy& e = eddies_[eddyI];
        e.move(deltaT);

        const scalar position0 = e.x();

        // Check to see if eddy has exited downstream box plane
        if (fabs(position0) > maxSigmaX_)
        {
            pointIndexHit pos(setNewPosition(false));
            label faceI = pos.index();

            scalar Umax = gMax(mag(U_&patchNormal_));
            scalar Uran = rndGen_.scalar01()*Umax;

            label iter = 0;

            while (Uran > mag(U_[faceI]&patchNormal_) && iter++ < seedIterMax_)
            {
                pointIndexHit pos(setNewPosition(false));
                faceI = pos.index();
                Uran = rndGen_.scalar01()*Umax;
            }

            e = eddy
            (
                velocityShape_,
                faceI,
                pos.hitPoint(),
                deltaT*(U_[faceI]&patchNormal_)-maxSigmaX_,
                U_[faceI]&patchNormal_,
                sigmaRatio*Lu_[faceI],
                sigmaRatio*Lv_[faceI],
                sigmaRatio*Lw_[faceI],
                R_[faceI],
                rndGen_
            );

            nRecycled++;
        }
    }

    reduce(nRecycled, sumOp<label>());

    if (nRecycled > 0)
    {
        Info<< "Patch: " << patch().patch().name() << " recycled "
            << nRecycled << " eddies" << endl;
    }
}


Foam::vector Foam::turbulentSEMInletFvPatchVectorField::uDashEddy
(
    const List<eddy>& eddies,
    const point& patchFaceCf
) const
{
    vector uDash(vector::zero);

    forAll(eddies, k)
    {
        const eddy& e = eddies[k];
        uDash += e.uDash(patchFaceCf, patchNormal_);
    }

    return uDash;
}


void Foam::turbulentSEMInletFvPatchVectorField::calcOverlappingProcEddies
(
    List<List<eddy>>& overlappingEddies
) const
{
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag + 1;

    List<boundBox> patchBBs(Pstream::nProcs());
    patchBBs[Pstream::myProcNo()] = patchBounds_;
    Pstream::gatherList(patchBBs);
    Pstream::scatterList(patchBBs);

    // Per processor indices into all segments to send
    List<DynamicList<label>> dynSendMap(Pstream::nProcs());

    forAll(eddies_, i)
    {
        // Collect overlapping eddies
        const eddy& e = eddies_[i];

        // Eddy bounds
        point x = e.position(patchNormal_);
        boundBox ebb = e.bounds();
        ebb.min() += x;
        ebb.max() += x;

        forAll(patchBBs, procI)
        {
            // Not including intersection with local patch
            if (procI != Pstream::myProcNo())
            {
                if (ebb.overlaps(patchBBs[procI]))
                {
                    dynSendMap[procI].append(i);
                }
            }
        }
    }

    labelListList sendMap(Pstream::nProcs());

    forAll(sendMap, procI)
    {
        sendMap[procI].transfer(dynSendMap[procI]);
    }

    // Send the number of eddies for local processors to receive
    labelListList sendSizes(Pstream::nProcs());
    sendSizes[Pstream::myProcNo()].setSize(Pstream::nProcs());

    forAll(sendMap, procI)
    {
        sendSizes[Pstream::myProcNo()][procI] = sendMap[procI].size();
    }

    Pstream::gatherList(sendSizes);
    Pstream::scatterList(sendSizes);

    // Determine order of receiving
    labelListList constructMap(Pstream::nProcs());

    // Local segment first
    constructMap[Pstream::myProcNo()] = identity
    (
        sendMap[Pstream::myProcNo()].size()
    );

    label segmentI = constructMap[Pstream::myProcNo()].size();

    forAll(constructMap, procI)
    {
        if (procI != Pstream::myProcNo())
        {
            // What I need to receive is what other processor is sending to me
            label nRecv = sendSizes[procI][Pstream::myProcNo()];
            constructMap[procI].setSize(nRecv);

            for (label i = 0; i < nRecv; i++)
            {
                constructMap[procI][i] = segmentI++;
            }
        }
    }

    mapDistribute map(segmentI, xferMove(sendMap), xferMove(constructMap));

    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    for (label domain = 0; domain < Pstream::nProcs(); domain++)
    {
        const labelList& sendElems = map.subMap()[domain];

        if (domain != Pstream::myProcNo() && sendElems.size())
        {
            List<eddy> subEddies(UIndirectList<eddy>(eddies_, sendElems));

            UOPstream toDomain(domain, pBufs);

            toDomain<< subEddies;
        }
    }

    // Start receiving
    pBufs.finishedSends();

    // Consume
    for (label domain = 0; domain < Pstream::nProcs(); domain++)
    {
        const labelList& recvElems = map.constructMap()[domain];

        if (domain != Pstream::myProcNo() && recvElems.size())
        {
            UIPstream str(domain, pBufs);
            {
                str >> overlappingEddies[domain];
            }
        }
    }

    // Restore tag
    UPstream::msgType() = oldTag;

}

Foam::scalar
Foam::turbulentSEMInletFvPatchVectorField::sigmaToLengthscaleRatio() const
{
    if (velocityShape_ == "tent")
    {
        return 4.0/3.0;
    }
    else if (velocityShape_ == "step")
    {
        return 1.0;
    }
    else if (velocityShape_ == "gaussian")
    {
        const scalar pi = constant::mathematical::pi;
        return sqrt(9.0/pi);
    }
    else
    {
        Info << "velocity shape " << velocityShape_ << "does not exist (ERROR)" << endl;
        return 0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbulentSEMInletFvPatchVectorField::
turbulentSEMInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    curTimeIndex_(-1),

    interpolateU_(false),
    interpolateR_(false),
    interpolateL_(false),
    U_(),
    R_(),
    Lu_(),
    Lv_(),
    Lw_(),

    patchBounds_(boundBox::invertedBox),
    patchArea_(-1),
    patchNormal_(vector::zero),
    triFace_(),
    triToFace_(),
    triCumulativeMagSf_(),
    sumTriMagSf_(Pstream::nProcs() + 1, 0.0),
    V0_(0),

    eddyDensity_(1),
    nEddy_(0),
    eddies_(0),
    velocityShape_("gaussian"),
    rndGen_((Pstream::myProcNo()+1)*time(NULL)),
    maxSigmaX_(0),
    singleProc_(false),

    perturb_(1e-5),
    mapMethod_("planarInterpolation"),
    mapperPtr_(nullptr)
{}


Foam::turbulentSEMInletFvPatchVectorField::
turbulentSEMInletFvPatchVectorField
(
    const turbulentSEMInletFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    curTimeIndex_(-1),

    interpolateU_(ptf.interpolateU_),
    interpolateR_(ptf.interpolateR_),
    interpolateL_(ptf.interpolateL_),
    U_(ptf.U_, mapper),
    R_(ptf.R_, mapper),
    Lu_(ptf.Lu_, mapper),
    Lv_(ptf.Lv_, mapper),
    Lw_(ptf.Lw_, mapper),

    patchBounds_(ptf.patchBounds_),
    patchArea_(ptf.patchArea_),
    patchNormal_(ptf.patchNormal_),
    triFace_(ptf.triFace_),
    triToFace_(ptf.triToFace_),
    triCumulativeMagSf_(ptf.triCumulativeMagSf_),
    sumTriMagSf_(ptf.sumTriMagSf_),
    V0_(ptf.V0_),

    eddyDensity_(ptf.eddyDensity_),
    nEddy_(ptf.nEddy_),
    eddies_(ptf.eddies_),
    velocityShape_(ptf.velocityShape_),
    rndGen_(ptf.rndGen_),
    maxSigmaX_(ptf.maxSigmaX_),
    singleProc_(ptf.singleProc_),

    perturb_(ptf.perturb_),
    mapMethod_(ptf.mapMethod_),
    mapperPtr_(nullptr)
{}


Foam::turbulentSEMInletFvPatchVectorField::
turbulentSEMInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    curTimeIndex_(-1),

    interpolateU_(false),
    interpolateR_(false),
    interpolateL_(false),
    U_(interpolateOrRead<vector>("U", dict, interpolateU_)),
    R_(interpolateOrRead<symmTensor>("R", dict, interpolateR_)),
    Lu_(interpolateOrRead<vector>("Lu", dict, interpolateL_)),
    Lv_(interpolateOrRead<vector>("Lv", dict, interpolateL_)),
    Lw_(interpolateOrRead<vector>("Lw", dict, interpolateL_)),

    patchBounds_(boundBox::invertedBox),
    patchArea_(-1),
    patchNormal_(vector::zero),
    triFace_(),
    triToFace_(),
    triCumulativeMagSf_(),
    sumTriMagSf_(Pstream::nProcs() + 1, 0.0),
    V0_(0),

    eddyDensity_(dict.lookupOrDefault<scalar>("eddyDensity", 1)),
    nEddy_(0),
    eddies_(),
    velocityShape_(dict.lookupOrDefault<word>("velocityShape", "gaussian")),
    rndGen_(0),
    maxSigmaX_(0),
    singleProc_(false),

    perturb_(dict.lookupOrDefault<scalar>("perturb", 1e-5)),
    mapMethod_(dict.lookup("mapMethod")),
    mapperPtr_(nullptr)
{
}


Foam::turbulentSEMInletFvPatchVectorField::
turbulentSEMInletFvPatchVectorField
(
    const turbulentSEMInletFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    curTimeIndex_(-1),

    interpolateU_(ptf.interpolateU_),
    interpolateR_(ptf.interpolateR_),
    interpolateL_(ptf.interpolateL_),
    U_(ptf.U_),
    R_(ptf.R_),
    Lu_(ptf.Lu_),
    Lv_(ptf.Lv_),
    Lw_(ptf.Lw_),

    patchBounds_(ptf.patchBounds_),
    patchArea_(ptf.patchArea_),
    patchNormal_(ptf.patchNormal_),
    triFace_(ptf.triFace_),
    triToFace_(ptf.triToFace_),
    triCumulativeMagSf_(ptf.triCumulativeMagSf_),
    sumTriMagSf_(ptf.sumTriMagSf_),
    V0_(ptf.V0_),

    eddyDensity_(ptf.eddyDensity_),
    nEddy_(ptf.nEddy_),
    eddies_(ptf.eddies_),
    velocityShape_(ptf.velocityShape_),
    rndGen_(ptf.rndGen_),
    maxSigmaX_(ptf.maxSigmaX_),
    singleProc_(ptf.singleProc_),

    perturb_(ptf.perturb_),
    mapMethod_(ptf.mapMethod_),
    mapperPtr_(nullptr)
{}


Foam::turbulentSEMInletFvPatchVectorField::
turbulentSEMInletFvPatchVectorField
(
    const turbulentSEMInletFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf),
    curTimeIndex_(-1),

    interpolateU_(ptf.interpolateU_),
    interpolateR_(ptf.interpolateR_),
    interpolateL_(ptf.interpolateL_),
    U_(ptf.U_),
    R_(ptf.R_),
    Lu_(ptf.Lu_),
    Lv_(ptf.Lv_),
    Lw_(ptf.Lw_),

    patchBounds_(ptf.patchBounds_),
    patchArea_(ptf.patchArea_),
    patchNormal_(ptf.patchNormal_),
    triFace_(ptf.triFace_),
    triToFace_(ptf.triToFace_),
    triCumulativeMagSf_(ptf.triCumulativeMagSf_),
    sumTriMagSf_(ptf.sumTriMagSf_),
    V0_(ptf.V0_),

    eddyDensity_(ptf.eddyDensity_),
    nEddy_(ptf.nEddy_),
    eddies_(ptf.eddies_),
    velocityShape_(ptf.velocityShape_),
    rndGen_(ptf.rndGen_),
    maxSigmaX_(ptf.maxSigmaX_),
    singleProc_(ptf.singleProc_),

    perturb_(ptf.perturb_),
    mapMethod_(ptf.mapMethod_),
    mapperPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::turbulentSEMInletFvPatchVectorField::~
turbulentSEMInletFvPatchVectorField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::turbulentSEMInletFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (curTimeIndex_ == -1)
    {
        initialisePatch();
        initialiseEddyBox();
        initialiseEddies();
    }

    if (curTimeIndex_ != db().time().timeIndex())
    {
        if (debug)
        {
            label n = eddies_.size();
            Info<< "Number of eddies: " << returnReduce(n, sumOp<label>())
                << endl;
        }

        const scalar deltaT = db().time().deltaTValue();

        // Move eddies using mean velocity
        convectEddies(deltaT);

        // Set velocity
        vectorField& U = *this;
        U = U_;

        const pointField& Cf = patch().Cf();

        // Apply normalisation coefficient
        const scalar c = Foam::sqrt(V0_)/Foam::sqrt(scalar(nEddy_));

        // In parallel, need to collect all eddies that will interact with
        // local faces

        if (singleProc_ || !Pstream::parRun())
        {
            forAll(U, faceI)
            {
                U[faceI] += c*uDashEddy(eddies_, Cf[faceI]);
            }
        }
        else
        {
            // Process local eddy contributions
            forAll(U, faceI)
            {
                U[faceI] += c*uDashEddy(eddies_, Cf[faceI]);
            }

            // Add contributions from overlapping eddies
            List<List<eddy>> overlappingEddies(Pstream::nProcs());
            calcOverlappingProcEddies(overlappingEddies);

            forAll(overlappingEddies, procI)
            {
                const List<eddy>& eddies = overlappingEddies[procI];

                if (eddies.size())
                {
                    //Pout<< "Applying " << eddies.size()
                    //    << " eddies from processor " << procI << endl;

                    forAll(U, faceI)
                    {
                        U[faceI] += c*uDashEddy(eddies, Cf[faceI]);
                    }
                }
            }
        }

        // Re-scale to ensure correct flow rate
        scalar fCorr = gSum(U_&patch().Sf())/gSum(U&patch().Sf());

        U *= fCorr;

        if (Pstream::master())
        {
            Info<< "mass flow correction coefficient: " << fCorr << endl;
        }

        curTimeIndex_ = db().time().timeIndex();
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::turbulentSEMInletFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    writeEntry("value", os);

    os.writeKeyword("eddyDensity") << eddyDensity_ << token::END_STATEMENT << nl;
    os.writeKeyword("velocityShape") << velocityShape_ << token::END_STATEMENT << nl;
    os.writeKeyword("perturb") << perturb_ << token::END_STATEMENT << nl;

    if (!interpolateR_)
    {
        R_.writeEntry("R", os);
    }

    if (!interpolateL_)
    {
        Lu_.writeEntry("Lu", os);
        Lv_.writeEntry("Lv", os);
        Lw_.writeEntry("Lw", os);
    }

    if (!interpolateU_)
    {
        U_.writeEntry("U", os);
    }

    if (!mapMethod_.empty())
    {
        os.writeKeyword("mapMethod") << mapMethod_ << token::END_STATEMENT << nl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       turbulentSEMInletFvPatchVectorField
   );
}


// ************************************************************************* //
