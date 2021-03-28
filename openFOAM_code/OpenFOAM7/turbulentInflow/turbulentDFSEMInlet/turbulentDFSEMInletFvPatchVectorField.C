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

#include "turbulentDFSEMInletFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "momentOfInertia.H"
#include "IFstream.H"
#include "OFstream.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::label Foam::turbulentDFSEMInletFvPatchVectorField::seedIterMax_ = 1000;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::pointToPointPlanarInterpolation&
Foam::turbulentDFSEMInletFvPatchVectorField::patchMapper() const
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


void Foam::turbulentDFSEMInletFvPatchVectorField::initialisePatch()
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

    boundBox globalBounds = boundBox(patch.localPoints(), true);

    patchSpanY_ = globalBounds.span().y();
    patchSpanZ_ = globalBounds.span().z();

    patchBounds_.inflate(0.1);

    // Determine if all eddies spawned from a single processor
    singleProc_ = patch.size() == returnReduce(patch.size(), sumOp<label>());
    reduce(singleProc_, orOp<bool>());
}


void Foam::turbulentDFSEMInletFvPatchVectorField::initialiseEddyBox()
{
    const scalarField& magSf = patch().magSf();

    //const scalarField cellDx(Foam::sqrt(magSf));
    const scalarField cellDx(max(Foam::sqrt(magSf), 2/patch().deltaCoeffs()));

    // Inialise eddy box extents
    forAll(*this, faceI)
    {
        scalar& s = sigmax_[faceI];

        // Length scale in x direction (based on eq. 14)
        s = mag(L_[faceI]);
        s = min(s, kappa_*delta_);

        // Allow eddies to be smaller than the mesh scale as suggested by
        // the reference?
        // s = min(s, nCellPerEddy_*cellDx[faceI]);
        s = max(s, nCellPerEddy_*cellDx[faceI]);
    }

    // Maximum extent across all processors
    maxSigmaX_ = gMax(sigmax_);

    // Eddy box volume
    v0_ = 2*gSum(magSf)*maxSigmaX_;

    {
        Info<< "Patch: " << patch().patch().name() << " eddy box:" << nl
            << "    volume    : " << v0_ << nl
            << "    maxSigmaX : " << maxSigmaX_ << nl
            << endl;
    }
}


Foam::pointIndexHit Foam::turbulentDFSEMInletFvPatchVectorField::setNewPosition
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


void Foam::turbulentDFSEMInletFvPatchVectorField::initialiseEddies()
{
    // Initialise eddy properties
    scalar sumVolEddy = 0;
    scalar sumVolEddyAllProc = 0;

    if (nEddyLocal_ && !isCleanRestart_)
    {
        eddies_.setSize(nEddyLocal_);

        forAll(eddies_, k)
        {    
            label faceI = eddyLabel_[k];

            dfeddy e
            (
                faceI,
                eddyPosition_[k],
                eddyDistance_[k],
                R_[faceI],
                eddyScale_[k],
                eddyIntensity_[k]
            );

            eddies_[k] = e;

            sumVolEddy += e.volume();
        }

        sumVolEddyAllProc = returnReduce(sumVolEddy, sumOp<scalar>());
    }
    else
    {
        DynamicList<dfeddy> eddies(size());

        while (sumVolEddyAllProc/v0_ < density_)
        {
            bool search = true;
            label iter = 0;

            while (search && iter++ < seedIterMax_)
            {
                // Get new parallel consistent position
                pointIndexHit pos(setNewPosition(true));
                label faceI = pos.index();

                // Note: only 1 processor will pick up this face
                if (faceI != -1)
                {
                    dfeddy e
                    (
                        faceI,
                        pos.hitPoint(),
                        rndGen_.scalarAB(-maxSigmaX_, maxSigmaX_),
                        L_[faceI],
                        R_[faceI],
                        rndGen_
                    );

                    // If eddy valid, patchFaceI is non-zero
                    if (e.patchFaceI() != -1)
                    {
                        eddies.append(e);
                        sumVolEddy += e.volume();
                        search = false;
                    }
                }
                // else eddy on remote processor

                reduce(search, andOp<bool>());
            }

            sumVolEddyAllProc = returnReduce(sumVolEddy, sumOp<scalar>());
        }

        eddies_.transfer(eddies);

        nEddyLocal_ = eddies_.size();

        eddyLabel_.setSize(nEddyLocal_);
        eddyPosition_.setSize(nEddyLocal_);
        eddyDistance_.setSize(nEddyLocal_);
        eddyScale_.setSize(nEddyLocal_);
        eddyIntensity_.setSize(nEddyLocal_);

        forAll(eddies_, k)
        {
            const dfeddy& e = eddies_[k];
            eddyLabel_[k] = e.patchFaceI();
            eddyPosition_[k] = e.position0();
            eddyDistance_[k] = e.x();
            eddyScale_[k] = e.sigma();
            eddyIntensity_[k] = e.alpha();
        }
    }

    nEddyGlobal_ = nEddyLocal_;
    reduce(nEddyGlobal_, sumOp<label>());

    if (debug)
    {
        Pout<< "Patch:" << patch().patch().name();

        if (Pstream::parRun())
        {
            Pout<< " processor:" << Pstream::myProcNo();
        }

        Pout<< " seeded:" << nEddyGlobal_ << " eddies" << endl;
    }

    if (nEddyGlobal_ > 0)
    {
        Info<< "Turbulent DFDFSEM patch: " << patch().name()
            << " seeded " << nEddyGlobal_ << " eddies with total volume "
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


void Foam::turbulentDFSEMInletFvPatchVectorField::convectEddies
(
    const scalar deltaT
)
{
    // Note: all operations applied to local processor only

    label nRecycled = 0;

    forAll(eddies_, eddyI)
    {
        dfeddy& e = eddies_[eddyI];

        e.move(deltaT*U_[e.patchFaceI()]);

        eddyDistance_[eddyI] = e.x();

        const scalar position0 = e.x();

        // Check to see if eddy has exited downstream box plane
        if (position0 > maxSigmaX_)
        {
            bool search = true;
            label iter = 0;

            while (search && iter++ < seedIterMax_)
            {
               // Spawn new eddy with new random properties (intensity etc)
               pointIndexHit pos(setNewPosition(false));
               label faceI = pos.index();

               scalar Urand = rndGen_.scalar01()*UMax_;
               label iter0 = 0;

               while (Urand > U_[faceI] && iter0++ < seedIterMax_)
               {
                   pointIndexHit pos(setNewPosition(false));
                   faceI = pos.index();
                   Urand = rndGen_.scalar01()*UMax_;
               }

               e = dfeddy
                    (
                        faceI,
                        pos.hitPoint(),
                        -maxSigmaX_ + rndGen_.scalar01()*deltaT*U_[faceI],
                        L_[faceI],
                        R_[faceI],
                        rndGen_
                    );

                eddyLabel_[eddyI] = e.patchFaceI();
                eddyPosition_[eddyI] = e.position0();
                eddyDistance_[eddyI] = e.x();
                eddyScale_[eddyI] = e.sigma();
                eddyIntensity_[eddyI] = e.alpha();

                if (e.patchFaceI() != -1)
                {
                    search = false;
                }
            }

            nRecycled++;
        }
    }

    reduce(nRecycled, sumOp<label>());

    if (debug && nRecycled > 0)
    {
        Info<< "Patch: " << patch().patch().name() << " recycled "
            << nRecycled << " eddies" << endl;
    }
}


Foam::vectorField Foam::turbulentDFSEMInletFvPatchVectorField::uDashEddy
(
    const List<dfeddy>& eddies,
    const pointField& Cf
) const
{
    vectorField uDash(Cf.size(), vector::zero);

    forAll(eddies, k)
    {
        const dfeddy& e = eddies[k];
        uDash += e.uDash(Cf, patchNormal_);

        if (periodicInY_)
        {
            const vector yOffSet = vector(0, patchSpanY_, 0);
            uDash += e.uDash(Cf+yOffSet, patchNormal_);
            uDash += e.uDash(Cf-yOffSet, patchNormal_);
        }

        if (periodicInZ_)
        {
            const vector zOffSet = vector(0, 0, patchSpanZ_);
            uDash += e.uDash(Cf+zOffSet, patchNormal_);
            uDash += e.uDash(Cf-zOffSet, patchNormal_);
        }

        if (periodicInY_&&periodicInZ_)
        {
            const vector yOffSet = vector(0, patchSpanY_, 0);
            const vector zOffSet = vector(0, 0, patchSpanZ_);

            uDash += e.uDash(Cf+yOffSet+zOffSet, patchNormal_);
            uDash += e.uDash(Cf+yOffSet-zOffSet, patchNormal_);
            uDash += e.uDash(Cf-yOffSet+zOffSet, patchNormal_);
            uDash += e.uDash(Cf-yOffSet-zOffSet, patchNormal_);
        }
    }

    return uDash;
}


void Foam::turbulentDFSEMInletFvPatchVectorField::calcOverlappingProcEddies
(
    List<List<dfeddy>>& overlappingEddies
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

    if (!periodicInY_ && !periodicInZ_)
    {
        forAll(eddies_, i)
        {
            // Collect overlapping eddies
            const dfeddy& e = eddies_[i];

            // Vorton bounds
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
    }
    else if (periodicInY_ && !periodicInZ_)
    {
        forAll(eddies_, i)
        {
            // Collect overlapping eddies
            const dfeddy& e = eddies_[i];

            // Vorton bounds
            point x = e.position(patchNormal_);

            List<boundBox> lebb(3, e.bounds());

            lebb[0].min() += x;
            lebb[0].max() += x;

            lebb[1].min() += x+vector(0, patchSpanY_, 0);
            lebb[1].max() += x+vector(0, patchSpanY_, 0);

            lebb[2].min() += x-vector(0, patchSpanY_, 0);
            lebb[2].max() += x-vector(0, patchSpanY_, 0);

            forAll(patchBBs, procI)
            {
                // Not including intersection with local patch
                if (procI != Pstream::myProcNo())
                {
                    forAll(lebb, indi)
                    {
                        if (lebb[indi].overlaps(patchBBs[procI]))
                        {
                            dynSendMap[procI].append(i);
                            break;
                        }
                    }
                }      
            }
        }
    }
    else if (!periodicInY_ && periodicInZ_)
    {
        forAll(eddies_, i)
        {
            // Collect overlapping eddies
            const dfeddy& e = eddies_[i];

            // Vorton bounds
            point x = e.position(patchNormal_);

            List<boundBox> lebb(3, e.bounds());

            lebb[0].min() += x;
            lebb[0].max() += x;

            lebb[1].min() += x+vector(0, 0, patchSpanZ_);
            lebb[1].max() += x+vector(0, 0, patchSpanZ_);

            lebb[2].min() += x-vector(0, 0, patchSpanZ_);
            lebb[2].max() += x-vector(0, 0, patchSpanZ_);

            forAll(patchBBs, procI)
            {
                // Not including intersection with local patch
                if (procI != Pstream::myProcNo())
                {
                    forAll(lebb, indi)
                    {
                        if (lebb[indi].overlaps(patchBBs[procI]))
                        {
                            dynSendMap[procI].append(i);
                            break;
                        }
                    }
                }      
            }
        }
    }
    else
    {
        forAll(eddies_, i)
        {
            // Collect overlapping eddies
            const dfeddy& e = eddies_[i];

            // Vorton bounds
            point x = e.position(patchNormal_);

            List<boundBox> lebb(9, e.bounds());

            lebb[0].min() += x;
            lebb[0].max() += x;

            lebb[1].min() += x+vector(0, patchSpanY_, 0);
            lebb[1].max() += x+vector(0, patchSpanY_, 0);

            lebb[2].min() += x-vector(0, patchSpanY_, 0);
            lebb[2].max() += x-vector(0, patchSpanY_, 0);

            lebb[3].min() += x+vector(0, 0, patchSpanZ_);
            lebb[3].max() += x+vector(0, 0, patchSpanZ_);

            lebb[4].min() += x-vector(0, 0, patchSpanZ_);
            lebb[4].max() += x-vector(0, 0, patchSpanZ_);

            lebb[5].min() += x+vector(0, patchSpanY_, patchSpanZ_);
            lebb[5].max() += x+vector(0, patchSpanY_, patchSpanZ_);

            lebb[6].min() += x-vector(0, patchSpanY_, patchSpanZ_);
            lebb[6].max() += x-vector(0, patchSpanY_, patchSpanZ_);

            lebb[7].min() += x+vector(0, patchSpanY_,-patchSpanZ_);
            lebb[7].max() += x+vector(0, patchSpanY_,-patchSpanZ_);

            lebb[8].min() += x-vector(0, patchSpanY_,-patchSpanZ_);
            lebb[8].max() += x-vector(0, patchSpanY_,-patchSpanZ_);

            forAll(patchBBs, procI)
            {
                // Not including intersection with local patch
                if (procI != Pstream::myProcNo())
                {
                    forAll(lebb, indi)
                    {
                        if (lebb[indi].overlaps(patchBBs[procI]))
                        {
                            dynSendMap[procI].append(i);
                            break;
                        }
                    }
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

    mapDistribute map(segmentI, std::move(sendMap), std::move(constructMap));

    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    for (label domain = 0; domain < Pstream::nProcs(); domain++)
    {
        const labelList& sendElems = map.subMap()[domain];

        if (domain != Pstream::myProcNo() && sendElems.size())
        {
            List<dfeddy> subEddies(UIndirectList<dfeddy>(eddies_, sendElems));

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


void Foam::turbulentDFSEMInletFvPatchVectorField::initialiseParameters()
{   
    const vectorField Cf(patch().Cf());
    
    forAll(U_, label)
    {
        bool isPositiveDefinite(true);
        tensor Lund(tensor::zero);
        
        if (U_[label] < 0)
        {
            Pout << "error: the patch-normal velocity magnitude at the point " << Cf[label]
                 << " is no larger than 0, please modify the input parameters for U" << endl;
        }
        
        if (L_[label] < 0)
        {
            Pout << "error: the integral length scale at the point " << Cf[label]
                 << " is no larger than 0, please modify the input parameters for L" << endl;
        }
        
        if (R_[label].component(symmTensor::XX) < 0)
        {
            isPositiveDefinite = false;
        }
        else
        {
            Lund.component(tensor::XX) = sqrt(R_[label].component(symmTensor::XX));
            Lund.component(tensor::YX) = R_[label].component(symmTensor::XY)/Lund.component(tensor::XX);
            Lund.component(tensor::ZX) = R_[label].component(symmTensor::XZ)/Lund.component(tensor::XX);
            
            const scalar sqrLundYY = R_[label].component(symmTensor::YY)-sqr(Lund.component(tensor::YX));
            
            if (sqrLundYY < 0)
            {
                isPositiveDefinite = false;
            }
            else
            {
                Lund.component(tensor::YY) = sqrt(sqrLundYY);
                Lund.component(tensor::ZY) = (R_[label].component(symmTensor::YZ)-Lund.component(tensor::YX)*Lund.component(tensor::ZX))/Lund.component(tensor::YY);
                
                const scalar sqrLundZZ = R_[label].component(symmTensor::ZZ)-sqr(Lund.component(tensor::ZX))-sqr(Lund.component(tensor::ZY));
                
                if (sqrLundZZ < 0)
                {
                    isPositiveDefinite = false;
                }
                else
                {
                    Lund.component(tensor::ZZ) = sqrt(sqrLundZZ);
                }
            }
        }
        
        if (!isPositiveDefinite)
        {
            Pout << "error: the Reynolds stress tensor at the point " << Cf[label]
                 << " is not positive definite, please modify the input parameters for R" << endl;
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbulentDFSEMInletFvPatchVectorField::
turbulentDFSEMInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    delta_(Zero),
    density_(Zero),
    kappa_(Zero),

    perturb_(1e-5),
    mapMethod_("nearestCell"),
    mapperPtr_(nullptr),
    interpolateR_(false),
    interpolateL_(false),
    interpolateU_(false),
    R_(),
    L_(),
    U_(),
    UMean_(Zero),
    UMax_(Zero),

    patchArea_(-1),
    triFace_(),
    triToFace_(),
    triCumulativeMagSf_(),
    sumTriMagSf_(Pstream::nProcs() + 1, Zero),
    periodicInY_(false),
    periodicInZ_(false),
    patchSpanY_(0),
    patchSpanZ_(0),

    eddies_(),
    nEddyGlobal_(Zero),
    nEddyLocal_(Zero),
    eddyLabel_(),
    eddyPosition_(),
    eddyDistance_(),
    eddyIntensity_(),

    nCellPerEddy_(5),
    patchNormal_(Zero),
    v0_(Zero),
    rndGen_((Pstream::myProcNo()+1)*time(NULL)),
    sigmax_(size(), Zero),
    maxSigmaX_(Zero),
    curTimeIndex_(-1),
    patchBounds_(boundBox::invertedBox),
    singleProc_(false),
    isCleanRestart_(false),
    isRestart_(false)
{}


Foam::turbulentDFSEMInletFvPatchVectorField::
turbulentDFSEMInletFvPatchVectorField
(
    const turbulentDFSEMInletFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    delta_(ptf.delta_),
    density_(ptf.density_),
    kappa_(ptf.kappa_),

    perturb_(ptf.perturb_),
    mapMethod_(ptf.mapMethod_),
    mapperPtr_(nullptr),
    interpolateR_(ptf.interpolateR_),
    interpolateL_(ptf.interpolateL_),
    interpolateU_(ptf.interpolateU_),
    R_(mapper(ptf.R_)),
    L_(mapper(ptf.L_)),
    U_(mapper(ptf.U_)),
    UMean_(ptf.UMean_),
    UMax_(ptf.UMax_),

    patchArea_(ptf.patchArea_),
    triFace_(ptf.triFace_),
    triToFace_(ptf.triToFace_),
    triCumulativeMagSf_(ptf.triCumulativeMagSf_),
    sumTriMagSf_(ptf.sumTriMagSf_),
    periodicInY_(ptf.periodicInY_),
    periodicInZ_(ptf.periodicInZ_),
    patchSpanY_(ptf.patchSpanY_),
    patchSpanZ_(ptf.patchSpanZ_),

    eddies_(),
    nEddyGlobal_(Zero),
    nEddyLocal_(Zero),
    eddyLabel_(),
    eddyPosition_(),
    eddyDistance_(),
    eddyScale_(),
    eddyIntensity_(),

    nCellPerEddy_(ptf.nCellPerEddy_),
    patchNormal_(ptf.patchNormal_),
    v0_(ptf.v0_),
    rndGen_(ptf.rndGen_),
    sigmax_(mapper(sigmax_)),
    maxSigmaX_(ptf.maxSigmaX_),
    curTimeIndex_(ptf.curTimeIndex_),
    patchBounds_(ptf.patchBounds_),
    singleProc_(ptf.singleProc_),
    isCleanRestart_(ptf.isCleanRestart_),
    isRestart_(ptf.isRestart_)
{}


Foam::turbulentDFSEMInletFvPatchVectorField::
turbulentDFSEMInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    delta_(readScalar(dict.lookup("delta"))),
    density_(dict.lookupOrDefault<scalar>("density", 1.0)),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),

    perturb_(dict.lookupOrDefault<scalar>("perturb", 1e-5)),
    mapMethod_(dict.lookupOrDefault<word>("mapMethod", "nearestCell")),
    mapperPtr_(nullptr),
    interpolateR_(dict.lookupOrDefault<bool>("interpolateR", false)),
    interpolateL_(dict.lookupOrDefault<bool>("interpolateL", false)),
    interpolateU_(dict.lookupOrDefault<bool>("interpolateU", false)),
    R_(interpolateOrRead<symmTensor>("R", dict, interpolateR_)),
    L_(interpolateOrRead<scalar>("L", dict, interpolateL_)),
    U_(interpolateOrRead<scalar>("U", dict, interpolateU_)),
    UMean_(0),
    UMax_(0),

    patchArea_(-1),
    triFace_(),
    triToFace_(),
    triCumulativeMagSf_(),
    sumTriMagSf_(Pstream::nProcs() + 1, Zero),
    periodicInY_(dict.lookupOrDefault<bool>("periodicInY", false)),
    periodicInZ_(dict.lookupOrDefault<bool>("periodicInZ", false)),
    patchSpanY_(0),
    patchSpanZ_(0),

    eddies_(),
    nEddyGlobal_(dict.lookupOrDefault<label>("nEddy", 0)),
    nEddyLocal_(dict.lookupOrDefault<label>("nEddyLocal", 0)),
    eddyLabel_(),
    eddyPosition_(),
    eddyDistance_(),
    eddyScale_(),
    eddyIntensity_(),

    nCellPerEddy_(dict.lookupOrDefault<label>("nCellPerEddy", 5)),
    patchNormal_(Zero),
    v0_(Zero),
    rndGen_((Pstream::myProcNo()+1)*time(NULL)),
    sigmax_(size(), Zero),
    maxSigmaX_(Zero),
    curTimeIndex_(-1),
    patchBounds_(boundBox::invertedBox),
    singleProc_(false),
    isCleanRestart_(dict.lookupOrDefault<bool>("cleanRestart", false)),
    isRestart_(false)
{
    // Set UMean as patch area average value
    UMean_ = gSum(U_*patch().magSf())/(gSum(patch().magSf()) + ROOTVSMALL);
    UMax_ = gMax(U_);

    if (nEddyLocal_ && !isCleanRestart_)
    {
        isRestart_ = true;

        ITstream& is = dict.lookup("eddyLabel");
        is >> static_cast<List<label>&>(eddyLabel_);

        is = dict.lookup("eddyPosition");
        is >> static_cast<List<vector>&>(eddyPosition_);

        is = dict.lookup("eddyDistance");
        is >> static_cast<List<scalar>&>(eddyDistance_);

        is = dict.lookup("eddyScale");
        is >> static_cast<List<vector>&>(eddyScale_);

        is = dict.lookup("eddyIntensity");
        is >> static_cast<List<vector>&>(eddyIntensity_);
    }
}


Foam::turbulentDFSEMInletFvPatchVectorField::
turbulentDFSEMInletFvPatchVectorField
(
    const turbulentDFSEMInletFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    delta_(ptf.delta_),
    density_(ptf.density_),
    kappa_(ptf.kappa_),

    perturb_(ptf.perturb_),
    mapMethod_(ptf.mapMethod_),
    mapperPtr_(nullptr),
    interpolateR_(ptf.interpolateR_),
    interpolateL_(ptf.interpolateL_),
    interpolateU_(ptf.interpolateU_),
    R_(ptf.R_),
    L_(ptf.L_),
    U_(ptf.U_),
    UMean_(ptf.UMean_),
    UMax_(ptf.UMax_),

    patchArea_(ptf.patchArea_),
    triFace_(ptf.triFace_),
    triToFace_(ptf.triToFace_),
    triCumulativeMagSf_(ptf.triCumulativeMagSf_),
    sumTriMagSf_(ptf.sumTriMagSf_),
    periodicInY_(ptf.periodicInY_),
    periodicInZ_(ptf.periodicInZ_),
    patchSpanY_(ptf.patchSpanY_),
    patchSpanZ_(ptf.patchSpanZ_),

    eddies_(),
    nEddyGlobal_(Zero),
    nEddyLocal_(Zero),
    eddyLabel_(),
    eddyPosition_(),
    eddyDistance_(),
    eddyScale_(),
    eddyIntensity_(),

    nCellPerEddy_(ptf.nCellPerEddy_),
    patchNormal_(ptf.patchNormal_),
    v0_(ptf.v0_),
    rndGen_(ptf.rndGen_),
    sigmax_(ptf.sigmax_),
    maxSigmaX_(ptf.maxSigmaX_),
    curTimeIndex_(ptf.curTimeIndex_),
    patchBounds_(ptf.patchBounds_),
    singleProc_(ptf.singleProc_),
    isCleanRestart_(ptf.isCleanRestart_),
    isRestart_(ptf.isRestart_)
{}


Foam::turbulentDFSEMInletFvPatchVectorField::
turbulentDFSEMInletFvPatchVectorField
(
    const turbulentDFSEMInletFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    delta_(ptf.delta_),
    density_(ptf.density_),
    kappa_(ptf.kappa_),

    perturb_(ptf.perturb_),
    mapMethod_(ptf.mapMethod_),
    mapperPtr_(nullptr),
    interpolateR_(ptf.interpolateR_),
    interpolateL_(ptf.interpolateL_),
    interpolateU_(ptf.interpolateU_),
    R_(ptf.R_),
    L_(ptf.L_),
    U_(ptf.U_),
    UMean_(ptf.UMean_),
    UMax_(ptf.UMax_),

    patchArea_(ptf.patchArea_),
    triFace_(ptf.triFace_),
    triToFace_(ptf.triToFace_),
    triCumulativeMagSf_(ptf.triCumulativeMagSf_),
    sumTriMagSf_(ptf.sumTriMagSf_),
    periodicInY_(ptf.periodicInY_),
    periodicInZ_(ptf.periodicInZ_),
    patchSpanY_(ptf.patchSpanY_),
    patchSpanZ_(ptf.patchSpanZ_),

    eddies_(),
    nEddyGlobal_(Zero),
    nEddyLocal_(Zero),
    eddyLabel_(),
    eddyPosition_(),
    eddyDistance_(),
    eddyScale_(),
    eddyIntensity_(),

    nCellPerEddy_(ptf.nCellPerEddy_),
    patchNormal_(ptf.patchNormal_),
    v0_(ptf.v0_),
    rndGen_(ptf.rndGen_),
    sigmax_(ptf.sigmax_),
    maxSigmaX_(ptf.maxSigmaX_),
    curTimeIndex_(ptf.curTimeIndex_),
    patchBounds_(ptf.patchBounds_),
    singleProc_(ptf.singleProc_),
    isCleanRestart_(ptf.isCleanRestart_),
    isRestart_(ptf.isRestart_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::turbulentDFSEMInletFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<vector>::autoMap(m);

    // Clear interpolator
    mapperPtr_.clear();

    m(U_, U_);
    m(R_, R_);
    m(L_, L_);
}


void Foam::turbulentDFSEMInletFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<vector>::rmap(ptf, addr);

    const turbulentDFSEMInletFvPatchVectorField& dfsemptf =
        refCast<const turbulentDFSEMInletFvPatchVectorField>(ptf);

    R_.rmap(dfsemptf.R_, addr);
    L_.rmap(dfsemptf.L_, addr);
    U_.rmap(dfsemptf.U_, addr);

    // Clear interpolator
    mapperPtr_.clear();
}


void Foam::turbulentDFSEMInletFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (curTimeIndex_ == -1)
    {
        initialisePatch();
        initialiseParameters();
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
        //U = UMean_;
        U = U_*patchNormal_;

        const pointField& Cf = patch().Cf();

        // Apply second part of normalisation coefficient
        const scalar c = Foam::sqrt(v0_)/Foam::sqrt(scalar(nEddyGlobal_));

        // In parallel, need to collect all eddies that will interact with
        // local faces

        if (singleProc_ || !Pstream::parRun())
        {
            U += c*uDashEddy(eddies_, Cf);
        }
        else
        {
            // Process local eddy contributions
            U += c*uDashEddy(eddies_, Cf);

            // Add contributions from overlapping eddies
            List<List<dfeddy>> overlappingEddies(Pstream::nProcs());
            calcOverlappingProcEddies(overlappingEddies);

            forAll(overlappingEddies, procI)
            {
                const List<dfeddy>& eddies = overlappingEddies[procI];

                if (eddies.size())
                {
                    //Pout<< "Applying " << eddies.size()
                    //    << " eddies from processor " << procI << endl;

                    U += c*uDashEddy(eddies, Cf);
                }
            }
        }

        // Re-scale to ensure correct flow rate
        scalar fCorr = gSum(UMean_ * patch().magSf())/gSum(U & -patch().Sf());

        if (Pstream::master())
        {
            Info<< "mass flow correction coefficient: " << fCorr << endl;
        }

        U *= fCorr;

        if (debug)
        {
            Info<< "Patch:" << patch().patch().name()
                << " min/max(U):" << gMin(U) << ", " << gMax(U) << endl;
        }

        curTimeIndex_ = db().time().timeIndex();
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::turbulentDFSEMInletFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    writeEntry(os, "value", *this);

    writeEntryIfDifferent<bool>(os, "periodicInY", false, periodicInY_);
    writeEntryIfDifferent<bool>(os, "periodicInZ", false, periodicInZ_);

    writeEntry(os, "U", U_);
    writeEntry(os, "R", R_);
    writeEntry(os, "L", L_);

    os.writeKeyword("delta") << delta_ << token::END_STATEMENT << nl;

    writeEntryIfDifferent<scalar>(os, "density", 1.0, density_);
    writeEntryIfDifferent<scalar>(os, "kappa", 0.41, kappa_);
    writeEntryIfDifferent<scalar>(os, "perturb", 1e-5, perturb_);
    writeEntryIfDifferent<label>(os, "nCellPerEddy", 5, nCellPerEddy_);

    writeEntryIfDifferent<label>(os, "nEddy", 0, nEddyGlobal_);
    writeEntryIfDifferent<label>(os, "nEddyLocal", 0, nEddyLocal_);

    if (nEddyLocal_)
    {
        writeEntry(os, "eddyLabel", eddyLabel_);
        writeEntry(os, "eddyPosition", eddyPosition_);
        writeEntry(os, "eddyDistance", eddyDistance_);
        writeEntry(os, "eddyScale", eddyScale_);
        writeEntry(os, "eddyIntensity", eddyIntensity_);
    }

    if (!mapMethod_.empty())
    {
        writeEntryIfDifferent<word>
        (
            os,
            "mapMethod",
            "nearestCell",
            mapMethod_
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       turbulentDFSEMInletFvPatchVectorField
   );
}


// ************************************************************************* //
