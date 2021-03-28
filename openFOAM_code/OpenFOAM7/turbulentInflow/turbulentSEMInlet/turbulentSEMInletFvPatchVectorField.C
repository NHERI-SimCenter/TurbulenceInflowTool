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
#include "IFstream.H"
#include "OFstream.H"
#include "IOmanip.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::label Foam::turbulentSEMInletFvPatchVectorField::seedIterMax_ = 1000;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::turbulentSEMInletFvPatchVectorField::createFiles()
{
    if (nOutputFace_ > 0)
    {
        fileName baseDir = db().time().path();

        word outputPrefix = "postProcessing";
        word prefix = "turbulentSEMInlet";

        if (Pstream::parRun())
        {
            baseDir = baseDir/".."/outputPrefix;
        }
        else
        {
            baseDir = baseDir/outputPrefix;
        }

        const word startTimeName = db().time().timeName(db().time().startTime().value());

        label start = 0;

        for (label i=0; i<Pstream::myProcNo(); i++)
        {
            start = start + patchSize_[i];
        }

        forAll(outputFaceIndices_, i)
        {
            if (!filePtrs_.set(i))
            {
                fileName outputDir(baseDir/prefix/startTimeName);

                mkDir(outputDir);

                word fName(db().time().timeName(outputFaceIndices_[i]+start));

                // Check if file already exists
                IFstream is(outputDir/(fName + ".dat"));

                filePtrs_.set(i, new OFstream(outputDir/(fName + ".dat")));

                initStream(filePtrs_[i]);

                writeFileHeader(i);
            }
        }
    }
}

void Foam::turbulentSEMInletFvPatchVectorField::writeFileHeader(const label i)
{
    label start = 0;

    for (label i=0; i<Pstream::myProcNo(); i++)
    {
        start = start + patchSize_[i];
    }

    writeHeaderValue(filePtrs_[i], "face index", outputFaceIndices_[i]+start);
    writeHeaderValue(filePtrs_[i], "face centroid", patch().Cf()[outputFaceIndices_[i]]);
    writeHeaderValue(filePtrs_[i], "mean velocity magnitude", U_[outputFaceIndices_[i]]);
    writeHeaderValue(filePtrs_[i], "Reynolds stress tensor", R_[outputFaceIndices_[i]]);
    writeHeaderValue(filePtrs_[i], "integral length scales", L_[outputFaceIndices_[i]]);
    writeHeaderValue(filePtrs_[i], "synthetic method", "synthetic eddy");
    writeHeaderValue(filePtrs_[i], "eddy shape function", eddyType_);

    writeCommented(filePtrs_[i], "Time");

    writeTabbed(filePtrs_[i], "ux");
    writeTabbed(filePtrs_[i], "uy");
    writeTabbed(filePtrs_[i], "uz");

    filePtrs_[i] << endl;
}

void Foam::turbulentSEMInletFvPatchVectorField::writeValues(const label i, const vector u)
{
    writeTime(filePtrs_[i]);

    filePtrs_[i]
        << tab << u.component(vector::X)
        << tab << u.component(vector::Y)
        << tab << u.component(vector::Z) << endl;
}


Foam::label Foam::turbulentSEMInletFvPatchVectorField::charWidth() const
{
    label addChars = 8;
    return IOstream::defaultPrecision() + addChars;
}

void Foam::turbulentSEMInletFvPatchVectorField::initStream(Ostream& os) const
{
    os.setf(ios_base::scientific, ios_base::floatfield);
    os.width(charWidth());
}

void Foam::turbulentSEMInletFvPatchVectorField::writeCommented
(
    Ostream& os,
    const string& str
) const
{
    os  << setw(1) << "#" << setw(1) << ' '
        << setf(ios_base::left) << setw(charWidth() - 2) << str.c_str();
}

void Foam::turbulentSEMInletFvPatchVectorField::writeTabbed
(
    Ostream& os,
    const string& str
) const
{
    os  << tab << setw(charWidth()) << str.c_str();
}

void Foam::turbulentSEMInletFvPatchVectorField::writeHeader
(
    Ostream& os,
    const string& str
) const
{
    os  << setw(1) << "#" << setw(1) << ' '
        << setf(ios_base::left) << setw(charWidth() - 2) << str.c_str() << nl;
}

template<class Type>
void Foam::turbulentSEMInletFvPatchVectorField::writeHeaderValue
(
    Ostream& os,
    const string& property,
    const Type& value
) const
{
    os  << setw(1) << '#' << setw(1) << ' '
        << setf(ios_base::left) << setw(charWidth() - 2) << property.c_str()
        << setw(1) << ':' << setw(1) << ' ' << value << nl;
}

void Foam::turbulentSEMInletFvPatchVectorField::writeTime(Ostream& os) const
{
    os  << setw(charWidth()) << db().time().timeName();
}

void Foam::turbulentSEMInletFvPatchVectorField::initialiseOutput()
{
    sort(outputFaceIndices_);

    if (!Pstream::parRun())
    {
        filePtrs_.setSize(nOutputFace_);
    }
    else
    {
        label start = 0;

        for (label i=0; i<Pstream::myProcNo(); i++)
        {
            start = start + patchSize_[i];
        }

        label currentNumb = 0;
        label currentStar = 0;

        forAll(outputFaceIndices_, i)
        {
            if (outputFaceIndices_[i]>=start&&outputFaceIndices_[i]<start+patchSize_[Pstream::myProcNo()])
            {
                if (currentNumb == 0)
                {
                    currentStar = i;
                }
                currentNumb += 1;
            }
        }

        labelList temp = outputFaceIndices_;

        if (currentNumb > 0)
        {
            nOutputFace_ = currentNumb;

            filePtrs_.setSize(nOutputFace_);
            outputFaceIndices_.setSize(nOutputFace_);

            forAll(outputFaceIndices_, i)
            {
                outputFaceIndices_[i] = temp[currentStar+i]-start;
            }
        }
        else
        {
            filePtrs_.setSize(0);
            outputFaceIndices_.setSize(0);
        }
    }

    createFiles();
}

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

    // Set number of patch faces for each processor
    patchSize_[Pstream::myProcNo()] = nf.size();

    Pstream::gatherList(patchSize_);
    Pstream::scatterList(patchSize_);

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


void Foam::turbulentSEMInletFvPatchVectorField::initialiseEddyBox()
{
    const scalarField& magSf = patch().magSf();

    // Maximum extent across all processors
    maxSigmaX_ = 2.0*cmptMax
    (
        vector
        (
            gMax(L_.component(tensor::XX)), 
            gMax(L_.component(tensor::XY)), 
            gMax(L_.component(tensor::XZ))
        )
    );

    // Eddy box volume
    v0_ = 2*gSum(magSf)*maxSigmaX_;

    if (debug)
    {
        Info<< "Patch: " << patch().patch().name() << " eddy box:" << nl
            << "    volume    : " << v0_ << nl
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
    // Initialise eddy properties
    scalar sumVolEddy = 0;
    scalar sumVolEddyAllProc = 0;

    if (nEddyLocal_ && !isCleanRestart_)
    {
        eddies_.setSize(nEddyLocal_);

        forAll(eddies_, k)
        {    
            label faceI = eddyLabel_[k];

            eddy e
            (
                eddyType_,
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
        DynamicList<eddy> eddies(size());

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
                    eddy e
                    (
                        eddyType_,
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
            const eddy& e = eddies_[k];
            eddyLabel_[k] = e.patchFaceI();
            eddyPosition_[k] = e.position0();
            eddyDistance_[k] = e.x();
            eddyScale_[k] = e.sigma();
            eddyIntensity_[k] = e.gamma();
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
        Info<< "Turbulent SEM patch: " << patch().name()
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


void Foam::turbulentSEMInletFvPatchVectorField::convectEddies
(
    const scalar deltaT
)
{
    // Note: all operations applied to local processor only

    label nRecycled = 0;

    forAll(eddies_, eddyI)
    {
        eddy& e = eddies_[eddyI];

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

               e = eddy
                    (
                        eddyType_,
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
                eddyIntensity_[eddyI] = e.gamma();

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


Foam::vectorField Foam::turbulentSEMInletFvPatchVectorField::uDashEddy
(
    const List<eddy>& eddies,
    const pointField& Cf
) const
{
    vectorField uDash(Cf.size(), vector::zero);

    forAll(eddies, k)
    {
        const eddy& e = eddies[k];
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

    if (!periodicInY_ && !periodicInZ_)
    {
        forAll(eddies_, i)
        {
            // Collect overlapping eddies
            const eddy& e = eddies_[i];

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
            const eddy& e = eddies_[i];

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
            const eddy& e = eddies_[i];

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
            const eddy& e = eddies_[i];

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


void Foam::turbulentSEMInletFvPatchVectorField::initialiseParameters()
{   
    const vectorField Cf(patch().Cf());
    
    forAll(U_, label)
    {
        tensor Lund(tensor::zero);
        
        if (U_[label] < 0)
        {
            Pout << "error: the patch-normal velocity magnitude at the point " << Cf[label]
                 << " is no larger than 0, please modify the input parameters for U" << endl;
        }
        
        if (R_[label].component(symmTensor::XX) < 0)
        {
            Pout << "error: the Reynolds stress component R_XX at the point " << Cf[label] 
                 << " is no larger than 0, please modify the input parameters for R" << endl;
        }
        else
        {
            Lund.component(tensor::XX) = sqrt(R_[label].component(symmTensor::XX));
            Lund.component(tensor::YX) = R_[label].component(symmTensor::XY)/Lund.component(tensor::XX);
            Lund.component(tensor::ZX) = R_[label].component(symmTensor::XZ)/Lund.component(tensor::XX);
            
            const scalar sqrLundYY = R_[label].component(symmTensor::YY)-sqr(Lund.component(tensor::YX));
            
            if (sqrLundYY < 0)
            {
                Pout << "error: the Reynolds stress component R_YY at the point " << Cf[label]
                     << " is no larger than the square of Lund_YX,"
                     << " Please modify the input parameters for R" << endl;
            }
            else
            {
                Lund.component(tensor::YY) = sqrt(sqrLundYY);
                Lund.component(tensor::ZY) = (R_[label].component(symmTensor::YZ)-Lund.component(tensor::YX)*Lund.component(tensor::ZX))/Lund.component(tensor::YY);
                
                const scalar sqrLundZZ = R_[label].component(symmTensor::ZZ)-sqr(Lund.component(tensor::ZX))-sqr(Lund.component(tensor::ZY));
                
                if (sqrLundZZ < 0)
                {
                    Pout << "error: the Reynolds stress component R_ZZ at the point " << Cf[label]
                         << " is no larger than sum of the squares of Lund_ZX and Lund_ZY,"
                         << " Please modify the input parameters for R" << endl;
                }
                else
                {
                    Lund.component(tensor::ZZ) = sqrt(sqrLundZZ);
                    
                    const tensor sqrLund = cmptMultiply(Lund, Lund);
                    const vector ex = sqrLund.x()/R_[label].xx();
                    const vector ey = sqrLund.y()/R_[label].yy();
                    const vector ez = sqrLund.z()/R_[label].zz();

                    const tensor E = tensor(ex, ey, ez);

                    const tensor L0 = inv(E)&L_[label];
        
                    forAll(L0, subLabel)
                    {
                        if (L0[subLabel] < 0)
                        {
                            Pout << "error: the " << subLabel+1 << "-th component of the converted length scales at the point " << Cf[label] 
                                 << " is no larger than 0, please modify the input parameters for L" << endl;
                        }
                    }
                }
            }
        }
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
    density_(Zero),

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
    patchSize_(Pstream::nProcs(), 0),
    triFace_(),
    triToFace_(),
    triCumulativeMagSf_(),
    sumTriMagSf_(Pstream::nProcs() + 1, Zero),
    periodicInY_(false),
    periodicInZ_(false),
    patchSpanY_(0),
    patchSpanZ_(0),

    eddies_(),
    eddyType_("gaussian"),
    nEddyGlobal_(Zero),
    nEddyLocal_(Zero),
    eddyLabel_(),
    eddyPosition_(),
    eddyDistance_(),
    eddyScale_(),
    eddyIntensity_(),

    patchNormal_(Zero),
    v0_(Zero),
    rndGen_((Pstream::myProcNo()+1)*time(NULL)),
    maxSigmaX_(Zero),
    curTimeIndex_(-1),
    patchBounds_(boundBox::invertedBox),
    singleProc_(false),
    isCleanRestart_(false),
    isRestart_(false),

    nOutputFace_(0),
    outputFaceIndices_(),
    filePtrs_()
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
    density_(ptf.density_),

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
    patchSize_(ptf.patchSize_),
    triFace_(ptf.triFace_),
    triToFace_(ptf.triToFace_),
    triCumulativeMagSf_(ptf.triCumulativeMagSf_),
    sumTriMagSf_(ptf.sumTriMagSf_),
    periodicInY_(ptf.periodicInY_),
    periodicInZ_(ptf.periodicInZ_),
    patchSpanY_(ptf.patchSpanY_),
    patchSpanZ_(ptf.patchSpanZ_),

    eddies_(),
    eddyType_(ptf.eddyType_),
    nEddyGlobal_(Zero),
    nEddyLocal_(Zero),
    eddyLabel_(),
    eddyPosition_(),
    eddyDistance_(),
    eddyScale_(),
    eddyIntensity_(),

    patchNormal_(ptf.patchNormal_),
    v0_(ptf.v0_),
    rndGen_(ptf.rndGen_),
    maxSigmaX_(ptf.maxSigmaX_),
    curTimeIndex_(ptf.curTimeIndex_),
    patchBounds_(ptf.patchBounds_),
    singleProc_(ptf.singleProc_),
    isCleanRestart_(ptf.isCleanRestart_),
    isRestart_(ptf.isRestart_),

    nOutputFace_(ptf.nOutputFace_),
    outputFaceIndices_(ptf.outputFaceIndices_),
    filePtrs_()
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
    density_(dict.lookupOrDefault<scalar>("density", 1.0)),

    perturb_(dict.lookupOrDefault<scalar>("perturb", 1e-5)),
    mapMethod_(dict.lookupOrDefault<word>("mapMethod", "nearestCell")),
    mapperPtr_(nullptr),
    interpolateR_(dict.lookupOrDefault<bool>("interpolateR", false)),
    interpolateL_(dict.lookupOrDefault<bool>("interpolateL", false)),
    interpolateU_(dict.lookupOrDefault<bool>("interpolateU", false)),
    R_(interpolateOrRead<symmTensor>("R", dict, interpolateR_)),
    L_(interpolateOrRead<tensor>("L", dict, interpolateL_)),
    U_(interpolateOrRead<scalar>("U", dict, interpolateU_)),
    UMean_(0),
    UMax_(0),

    patchArea_(-1),
    patchSize_(Pstream::nProcs(), 0),
    triFace_(),
    triToFace_(),
    triCumulativeMagSf_(),
    sumTriMagSf_(Pstream::nProcs() + 1, Zero),
    periodicInY_(dict.lookupOrDefault<bool>("periodicInY", false)),
    periodicInZ_(dict.lookupOrDefault<bool>("periodicInZ", false)),
    patchSpanY_(0),
    patchSpanZ_(0),

    eddies_(),
    eddyType_(dict.lookupOrDefault<word>("eddyType", "gaussian")),
    nEddyGlobal_(dict.lookupOrDefault<label>("nEddy", 0)),
    nEddyLocal_(dict.lookupOrDefault<label>("nEddyLocal", 0)),
    eddyLabel_(),
    eddyPosition_(),
    eddyDistance_(),
    eddyScale_(),
    eddyIntensity_(),

    patchNormal_(Zero),
    v0_(Zero),
    rndGen_((Pstream::myProcNo()+1)*time(NULL)),
    maxSigmaX_(Zero),
    curTimeIndex_(-1),
    patchBounds_(boundBox::invertedBox),
    singleProc_(false),
    isCleanRestart_(dict.lookupOrDefault<bool>("cleanRestart", false)),
    isRestart_(false),

    nOutputFace_(dict.lookupOrDefault<label>("nOutputFace", 0)),
    outputFaceIndices_(dict.lookupOrDefault<labelList>("outputFaceIndices", labelList(nOutputFace_, 0))),
    filePtrs_()
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
        is >> static_cast<List<tensor>&>(eddyScale_);

        is = dict.lookup("eddyIntensity");
        is >> static_cast<List<vector>&>(eddyIntensity_);
    }
}


Foam::turbulentSEMInletFvPatchVectorField::
turbulentSEMInletFvPatchVectorField
(
    const turbulentSEMInletFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    density_(ptf.density_),

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
    patchSize_(ptf.patchSize_),
    triFace_(ptf.triFace_),
    triToFace_(ptf.triToFace_),
    triCumulativeMagSf_(ptf.triCumulativeMagSf_),
    sumTriMagSf_(ptf.sumTriMagSf_),
    periodicInY_(ptf.periodicInY_),
    periodicInZ_(ptf.periodicInZ_),
    patchSpanY_(ptf.patchSpanY_),
    patchSpanZ_(ptf.patchSpanZ_),

    eddies_(),
    eddyType_(ptf.eddyType_),
    nEddyGlobal_(Zero),
    nEddyLocal_(Zero),
    eddyLabel_(),
    eddyPosition_(),
    eddyDistance_(),
    eddyScale_(),
    eddyIntensity_(),

    patchNormal_(ptf.patchNormal_),
    v0_(ptf.v0_),
    rndGen_(ptf.rndGen_),
    maxSigmaX_(ptf.maxSigmaX_),
    curTimeIndex_(ptf.curTimeIndex_),
    patchBounds_(ptf.patchBounds_),
    singleProc_(ptf.singleProc_),
    isCleanRestart_(ptf.isCleanRestart_),
    isRestart_(ptf.isRestart_),

    nOutputFace_(ptf.nOutputFace_),
    outputFaceIndices_(ptf.outputFaceIndices_),
    filePtrs_()
{}


Foam::turbulentSEMInletFvPatchVectorField::
turbulentSEMInletFvPatchVectorField
(
    const turbulentSEMInletFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    density_(ptf.density_),

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
    patchSize_(ptf.patchSize_),
    triFace_(ptf.triFace_),
    triToFace_(ptf.triToFace_),
    triCumulativeMagSf_(ptf.triCumulativeMagSf_),
    sumTriMagSf_(ptf.sumTriMagSf_),
    periodicInY_(ptf.periodicInY_),
    periodicInZ_(ptf.periodicInZ_),
    patchSpanY_(ptf.patchSpanY_),
    patchSpanZ_(ptf.patchSpanZ_),

    eddies_(),
    eddyType_(ptf.eddyType_),
    nEddyGlobal_(Zero),
    nEddyLocal_(Zero),
    eddyLabel_(),
    eddyPosition_(),
    eddyDistance_(),
    eddyScale_(),
    eddyIntensity_(),

    patchNormal_(ptf.patchNormal_),
    v0_(ptf.v0_),
    rndGen_(ptf.rndGen_),
    maxSigmaX_(ptf.maxSigmaX_),
    curTimeIndex_(ptf.curTimeIndex_),
    patchBounds_(ptf.patchBounds_),
    singleProc_(ptf.singleProc_),
    isCleanRestart_(ptf.isCleanRestart_),
    isRestart_(ptf.isRestart_),

    nOutputFace_(ptf.nOutputFace_),
    outputFaceIndices_(ptf.outputFaceIndices_),
    filePtrs_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::turbulentSEMInletFvPatchVectorField::autoMap
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


void Foam::turbulentSEMInletFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<vector>::rmap(ptf, addr);

    const turbulentSEMInletFvPatchVectorField& dfsemptf =
        refCast<const turbulentSEMInletFvPatchVectorField>(ptf);

    R_.rmap(dfsemptf.R_, addr);
    L_.rmap(dfsemptf.L_, addr);
    U_.rmap(dfsemptf.U_, addr);

    // Clear interpolator
    mapperPtr_.clear();
}


void Foam::turbulentSEMInletFvPatchVectorField::updateCoeffs()
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
        initialiseOutput();
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
            List<List<eddy>> overlappingEddies(Pstream::nProcs());
            calcOverlappingProcEddies(overlappingEddies);

            forAll(overlappingEddies, procI)
            {
                const List<eddy>& eddies = overlappingEddies[procI];

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

        forAll(outputFaceIndices_, i)
        {
            writeValues(i, U[outputFaceIndices_[i]]);
        }

        if (debug)
        {
            Info<< "Patch:" << patch().patch().name()
                << " min/max(U):" << gMin(U) << ", " << gMax(U) << endl;
        }

        curTimeIndex_ = db().time().timeIndex();
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::turbulentSEMInletFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    writeEntry(os, "value", *this);

    writeEntryIfDifferent<bool>(os, "periodicInY", false, periodicInY_);
    writeEntryIfDifferent<bool>(os, "periodicInZ", false, periodicInZ_);

    writeEntry(os, "U", U_);
    writeEntry(os, "R", R_);
    writeEntry(os, "L", L_);

    writeEntryIfDifferent<scalar>(os, "density", 1.0, density_);
    writeEntryIfDifferent<scalar>(os, "perturb", 1e-5, perturb_);

    writeEntryIfDifferent<label>(os, "nEddy", 0, nEddyGlobal_);
    writeEntryIfDifferent<label>(os, "nEddyLocal", 0, nEddyLocal_);
    writeEntryIfDifferent<word>(os, "eddyType", "gaussian", eddyType_);

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

    if (nOutputFace_ > 0)
    {
        os.writeKeyword("nOutputFace") << nOutputFace_ << token::END_STATEMENT << nl;
        writeEntry(os, "outputFaceIndices", outputFaceIndices_);
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
