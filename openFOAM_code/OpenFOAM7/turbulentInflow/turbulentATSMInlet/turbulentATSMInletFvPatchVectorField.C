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

#include "turbulentATSMInletFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "momentOfInertia.H"
#include "IFstream.H"
#include "OFstream.H"
#include "IOmanip.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::label Foam::turbulentATSMInletFvPatchVectorField::seedIterMax_ = 1000;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::turbulentATSMInletFvPatchVectorField::createFiles()
{
    if (nOutputFace_ > 0)
    {
        fileName baseDir = db().time().path();

        word outputPrefix = "postProcessing";
        word prefix = "turbulentATSMInlet";

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

void Foam::turbulentATSMInletFvPatchVectorField::writeFileHeader(const label i)
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
    writeHeaderValue(filePtrs_[i], "synthetic method", "turbulent spot");
    writeHeaderValue(filePtrs_[i], "vorton type", vortonType_);

    writeCommented(filePtrs_[i], "Time");

    writeTabbed(filePtrs_[i], "ux");
    writeTabbed(filePtrs_[i], "uy");
    writeTabbed(filePtrs_[i], "uz");

    filePtrs_[i] << endl;
}

void Foam::turbulentATSMInletFvPatchVectorField::writeValues(const label i, const vector u)
{
    writeTime(filePtrs_[i]);

    filePtrs_[i]
        << tab << u.component(vector::X)
        << tab << u.component(vector::Y)
        << tab << u.component(vector::Z) << endl;
}


Foam::label Foam::turbulentATSMInletFvPatchVectorField::charWidth() const
{
    label addChars = 8;
    return IOstream::defaultPrecision() + addChars;
}

void Foam::turbulentATSMInletFvPatchVectorField::initStream(Ostream& os) const
{
    os.setf(ios_base::scientific, ios_base::floatfield);
    os.width(charWidth());
}

void Foam::turbulentATSMInletFvPatchVectorField::writeCommented
(
    Ostream& os,
    const string& str
) const
{
    os  << setw(1) << "#" << setw(1) << ' '
        << setf(ios_base::left) << setw(charWidth() - 2) << str.c_str();
}

void Foam::turbulentATSMInletFvPatchVectorField::writeTabbed
(
    Ostream& os,
    const string& str
) const
{
    os  << tab << setw(charWidth()) << str.c_str();
}

void Foam::turbulentATSMInletFvPatchVectorField::writeHeader
(
    Ostream& os,
    const string& str
) const
{
    os  << setw(1) << "#" << setw(1) << ' '
        << setf(ios_base::left) << setw(charWidth() - 2) << str.c_str() << nl;
}

template<class Type>
void Foam::turbulentATSMInletFvPatchVectorField::writeHeaderValue
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

void Foam::turbulentATSMInletFvPatchVectorField::writeTime(Ostream& os) const
{
    os  << setw(charWidth()) << db().time().timeName();
}

void Foam::turbulentATSMInletFvPatchVectorField::initialiseOutput()
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
Foam::turbulentATSMInletFvPatchVectorField::patchMapper() const
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


void Foam::turbulentATSMInletFvPatchVectorField::initialisePatch()
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

    // Determine if all vortons spawned from a single processor
    singleProc_ = patch.size() == returnReduce(patch.size(), sumOp<label>());
    reduce(singleProc_, orOp<bool>());
}


void Foam::turbulentATSMInletFvPatchVectorField::initialiseVortonBox()
{
    const scalarField& magSf = patch().magSf();

    // Maximum extent across all processors
    maxSigmaX_ = 3.0/sqrt(constant::mathematical::pi)*cmptMax(gMax(L_));

    // Vorton box volume
    v0_ = 2*gSum(magSf)*maxSigmaX_;

    {
        Info<< "Patch: " << patch().patch().name() << " vorton box:" << nl
            << "    volume    : " << v0_ << nl
            << "    maxSigmaX : " << maxSigmaX_ << nl
            << endl;
    }
}


Foam::pointIndexHit Foam::turbulentATSMInletFvPatchVectorField::setNewPosition
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


void Foam::turbulentATSMInletFvPatchVectorField::initialiseVortons()
{
    // Initialise vorton properties
    scalar sumVolVorton = 0;
    scalar sumVolVortonAllProc = 0;

    if (nVortonLocal_ && !isCleanRestart_)
    {
        vortons_.setSize(nVortonLocal_);

        forAll(vortons_, k)
        {    
            label faceI = vortonLabel_[k];

            vorton v
            (
                vortonType_,
                faceI,
                vortonPosition_[k],
                vortonDistance_[k],
                R_[faceI],
                vortonScale_[k],
                vortonIntensity_[k]
            );

            vortons_[k] = v;

            sumVolVorton += v.volume();
        }

        sumVolVortonAllProc = returnReduce(sumVolVorton, sumOp<scalar>());
    }
    else
    {
        DynamicList<vorton> vortons(size());

        while (sumVolVortonAllProc/v0_ < density_)
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
                    vorton v
                    (
                        vortonType_,
                        faceI,
                        pos.hitPoint(),
                        rndGen_.scalarAB(-maxSigmaX_, maxSigmaX_),
                        L_[faceI],
                        R_[faceI],
                        rndGen_
                    );

                    // If vorton valid, patchFaceI is non-zero
                    if (v.patchFaceI() != -1)
                    {
                        vortons.append(v);
                        sumVolVorton += v.volume();
                        search = false;
                    }
                }
                // else vorton on remote processor

                reduce(search, andOp<bool>());
            }

            sumVolVortonAllProc = returnReduce(sumVolVorton, sumOp<scalar>());
        }

        vortons_.transfer(vortons);

        nVortonLocal_ = vortons_.size();

        vortonLabel_.setSize(nVortonLocal_);
        vortonPosition_.setSize(nVortonLocal_);
        vortonDistance_.setSize(nVortonLocal_);
        vortonScale_.setSize(nVortonLocal_);
        vortonIntensity_.setSize(nVortonLocal_);

        forAll(vortons_, k)
        {
            const vorton& v = vortons_[k];
            vortonLabel_[k] = v.patchFaceI();
            vortonPosition_[k] = v.position0();
            vortonDistance_[k] = v.x();
            vortonScale_[k] = v.sigma();
            vortonIntensity_[k] = v.gamma();
        }
    }

    nVortonGlobal_ = nVortonLocal_;
    reduce(nVortonGlobal_, sumOp<label>());

    if (debug)
    {
        Pout<< "Patch:" << patch().patch().name();

        if (Pstream::parRun())
        {
            Pout<< " processor:" << Pstream::myProcNo();
        }

        Pout<< " seeded:" << nVortonGlobal_ << " vortons" << endl;
    }

    if (nVortonGlobal_ > 0)
    {
        Info<< "Turbulent ATSM patch: " << patch().name()
            << " seeded " << nVortonGlobal_ << " vortons with total volume "
            << sumVolVortonAllProc
            << endl;
    }
    else
    {
        WarningInFunction
            << "Patch: " << patch().patch().name()
            << " on field " << internalField().name()
            << ": No vortons seeded - please check your set-up" << endl;
    }
}


void Foam::turbulentATSMInletFvPatchVectorField::convectVortons
(
    const scalar deltaT
)
{
    // Note: all operations applied to local processor only

    label nRecycled = 0;

    forAll(vortons_, vortonI)
    {
        vorton& v = vortons_[vortonI];

        v.move(deltaT*U_[v.patchFaceI()]);

        vortonDistance_[vortonI] = v.x();

        const scalar position0 = v.x();

        // Check to see if vorton has exited downstream box plane
        if (position0 > maxSigmaX_)
        {
            bool search = true;
            label iter = 0;

            while (search && iter++ < seedIterMax_)
            {
               // Spawn new vorton with new random properties (intensity etc)
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

               v = vorton
                    (
                        vortonType_,
                        faceI,
                        pos.hitPoint(),
                        -maxSigmaX_ + rndGen_.scalar01()*deltaT*U_[faceI],
                        L_[faceI],
                        R_[faceI],
                        rndGen_
                    );

                vortonLabel_[vortonI] = v.patchFaceI();
                vortonPosition_[vortonI] = v.position0();
                vortonDistance_[vortonI] = v.x();
                vortonScale_[vortonI] = v.sigma();
                vortonIntensity_[vortonI] = v.gamma();

                if (v.patchFaceI() != -1)
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
            << nRecycled << " vortons" << endl;
    }
}


Foam::vectorField Foam::turbulentATSMInletFvPatchVectorField::uDashVorton
(
    const List<vorton>& vortons,
    const pointField& Cf
) const
{
    vectorField uDash(Cf.size(), vector::zero);

    forAll(vortons, k)
    {
        const vorton& v = vortons[k];
        uDash += v.uDash(Cf, patchNormal_);

        if (periodicInY_)
        {
            const vector yOffSet = vector(0, patchSpanY_, 0);
            uDash += v.uDash(Cf+yOffSet, patchNormal_);
            uDash += v.uDash(Cf-yOffSet, patchNormal_);
        }

        if (periodicInZ_)
        {
            const vector zOffSet = vector(0, 0, patchSpanZ_);
            uDash += v.uDash(Cf+zOffSet, patchNormal_);
            uDash += v.uDash(Cf-zOffSet, patchNormal_);
        }

        if (periodicInY_&&periodicInZ_)
        {
            const vector yOffSet = vector(0, patchSpanY_, 0);
            const vector zOffSet = vector(0, 0, patchSpanZ_);

            uDash += v.uDash(Cf+yOffSet+zOffSet, patchNormal_);
            uDash += v.uDash(Cf+yOffSet-zOffSet, patchNormal_);
            uDash += v.uDash(Cf-yOffSet+zOffSet, patchNormal_);
            uDash += v.uDash(Cf-yOffSet-zOffSet, patchNormal_);
        }
    }

    return uDash;
}


void Foam::turbulentATSMInletFvPatchVectorField::calcOverlappingProcVortons
(
    List<List<vorton>>& overlappingVortons
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
        forAll(vortons_, i)
        {
            // Collect overlapping vortons
            const vorton& v = vortons_[i];

            // Vorton bounds
            point x = v.position(patchNormal_);

            boundBox vbb = v.bounds();
            vbb.min() += x;
            vbb.max() += x;

            forAll(patchBBs, procI)
            {
                // Not including intersection with local patch
                if (procI != Pstream::myProcNo())
                {
                    if (vbb.overlaps(patchBBs[procI]))
                    {
                        dynSendMap[procI].append(i);
                    }
                }      
            }
        }
    }
    else if (periodicInY_ && !periodicInZ_)
    {
        forAll(vortons_, i)
        {
            // Collect overlapping vortons
            const vorton& v = vortons_[i];

            // Vorton bounds
            point x = v.position(patchNormal_);

            List<boundBox> lvbb(3, v.bounds());

            lvbb[0].min() += x;
            lvbb[0].max() += x;

            lvbb[1].min() += x+vector(0, patchSpanY_, 0);
            lvbb[1].max() += x+vector(0, patchSpanY_, 0);

            lvbb[2].min() += x-vector(0, patchSpanY_, 0);
            lvbb[2].max() += x-vector(0, patchSpanY_, 0);

            forAll(patchBBs, procI)
            {
                // Not including intersection with local patch
                if (procI != Pstream::myProcNo())
                {
                    forAll(lvbb, indi)
                    {
                        if (lvbb[indi].overlaps(patchBBs[procI]))
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
        forAll(vortons_, i)
        {
            // Collect overlapping vortons
            const vorton& v = vortons_[i];

            // Vorton bounds
            point x = v.position(patchNormal_);

            List<boundBox> lvbb(3, v.bounds());

            lvbb[0].min() += x;
            lvbb[0].max() += x;

            lvbb[1].min() += x+vector(0, 0, patchSpanZ_);
            lvbb[1].max() += x+vector(0, 0, patchSpanZ_);

            lvbb[2].min() += x-vector(0, 0, patchSpanZ_);
            lvbb[2].max() += x-vector(0, 0, patchSpanZ_);

            forAll(patchBBs, procI)
            {
                // Not including intersection with local patch
                if (procI != Pstream::myProcNo())
                {
                    forAll(lvbb, indi)
                    {
                        if (lvbb[indi].overlaps(patchBBs[procI]))
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
        forAll(vortons_, i)
        {
            // Collect overlapping vortons
            const vorton& v = vortons_[i];

            // Vorton bounds
            point x = v.position(patchNormal_);

            List<boundBox> lvbb(9, v.bounds());

            lvbb[0].min() += x;
            lvbb[0].max() += x;

            lvbb[1].min() += x+vector(0, patchSpanY_, 0);
            lvbb[1].max() += x+vector(0, patchSpanY_, 0);

            lvbb[2].min() += x-vector(0, patchSpanY_, 0);
            lvbb[2].max() += x-vector(0, patchSpanY_, 0);

            lvbb[3].min() += x+vector(0, 0, patchSpanZ_);
            lvbb[3].max() += x+vector(0, 0, patchSpanZ_);

            lvbb[4].min() += x-vector(0, 0, patchSpanZ_);
            lvbb[4].max() += x-vector(0, 0, patchSpanZ_);

            lvbb[5].min() += x+vector(0, patchSpanY_, patchSpanZ_);
            lvbb[5].max() += x+vector(0, patchSpanY_, patchSpanZ_);

            lvbb[6].min() += x-vector(0, patchSpanY_, patchSpanZ_);
            lvbb[6].max() += x-vector(0, patchSpanY_, patchSpanZ_);

            lvbb[7].min() += x+vector(0, patchSpanY_,-patchSpanZ_);
            lvbb[7].max() += x+vector(0, patchSpanY_,-patchSpanZ_);

            lvbb[8].min() += x-vector(0, patchSpanY_,-patchSpanZ_);
            lvbb[8].max() += x-vector(0, patchSpanY_,-patchSpanZ_);

            forAll(patchBBs, procI)
            {
                // Not including intersection with local patch
                if (procI != Pstream::myProcNo())
                {
                    forAll(lvbb, indi)
                    {
                        if (lvbb[indi].overlaps(patchBBs[procI]))
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

    // Send the number of vortons for local processors to receive
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
            List<vorton> subVortons(UIndirectList<vorton>(vortons_, sendElems));

            UOPstream toDomain(domain, pBufs);

            toDomain<< subVortons;
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
                str >> overlappingVortons[domain];
            }
        }
    }

    // Restore tag
    UPstream::msgType() = oldTag;
}


void Foam::turbulentATSMInletFvPatchVectorField::initialiseParameters()
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
                    
                    // Principal stresses - eigenvalues returned in ascending order
                    vector Rp;

                    if (R_[label].component(symmTensor::XY)==0&&R_[label].component(symmTensor::XZ)==0&&R_[label].component(symmTensor::YZ)==0)
                    {
                        Rp = vector(R_[label].component(symmTensor::XX),R_[label].component(symmTensor::YY),R_[label].component(symmTensor::ZZ));
                    }
                    else
                    {
                        Rp = eigenValues(R_[label]);
                    }

                    // Eddy rotation from principal-to-global axes
                    // - given by the 3 eigenvectors of the Reynold stress tensor as rows in
                    //   the result tensor (transposed transformation tensor)
                    // - returned in ascending eigenvalue order
                    tensor Rpg = eigenVectors(R_[label], Rp);

                    if (debug)
                    {
                        Pout<< "Rpg & R & Rpg.T(): " << (Rpg & R_[label] & Rpg.T()) << endl;
                    }

                    const tensor sqrRpg = cmptMultiply(Rpg.T(), Rpg.T());

                    const vector ex = cmptMultiply(sqrRpg.x(), Rp)/R_[label].xx();
                    const vector ey = cmptMultiply(sqrRpg.y(), Rp)/R_[label].yy();
                    const vector ez = cmptMultiply(sqrRpg.z(), Rp)/R_[label].zz();

                    const tensor E = tensor(ex, ey, ez);

                    vector L0 = inv(E)&L_[label];

                    if (vortonType_ == "typeR" || L0.x() < 0.0)
                    {
                        L0.x() = (L0.y()*L0.z()*sqrt(Rp.x()))/(L0.y()*sqrt(Rp.z())+L0.z()*sqrt(Rp.y()));
                    }
                    
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
        
        if (!isPositiveDefinite)
        {
            Pout << "error: the Reynolds stress tensor at the point " << Cf[label]
                 << " is not positive definite, please modify the input parameters for R" << endl;
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbulentATSMInletFvPatchVectorField::
turbulentATSMInletFvPatchVectorField
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

    vortons_(),
    vortonType_("typeR"),
    nVortonGlobal_(Zero),
    nVortonLocal_(Zero),
    vortonLabel_(),
    vortonPosition_(),
    vortonDistance_(),
    vortonScale_(),
    vortonIntensity_(),

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


Foam::turbulentATSMInletFvPatchVectorField::
turbulentATSMInletFvPatchVectorField
(
    const turbulentATSMInletFvPatchVectorField& ptf,
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

    vortons_(),
    vortonType_(ptf.vortonType_),
    nVortonGlobal_(Zero),
    nVortonLocal_(Zero),
    vortonLabel_(),
    vortonPosition_(),
    vortonDistance_(),
    vortonScale_(),
    vortonIntensity_(),

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


Foam::turbulentATSMInletFvPatchVectorField::
turbulentATSMInletFvPatchVectorField
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
    L_(interpolateOrRead<vector>("L", dict, interpolateL_)),
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

    vortons_(),
    vortonType_(dict.lookupOrDefault<word>("vortonType", "typeR")),
    nVortonGlobal_(dict.lookupOrDefault<label>("nVorton", 0)),
    nVortonLocal_(dict.lookupOrDefault<label>("nVortonLocal", 0)),
    vortonLabel_(),
    vortonPosition_(),
    vortonDistance_(),
    vortonScale_(),
    vortonIntensity_(),

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

    if (nVortonLocal_ && !isCleanRestart_)
    {
        isRestart_ = true;

        ITstream& is = dict.lookup("vortonLabel");
        is >> static_cast<List<label>&>(vortonLabel_);

        is = dict.lookup("vortonPosition");
        is >> static_cast<List<vector>&>(vortonPosition_);

        is = dict.lookup("vortonDistance");
        is >> static_cast<List<scalar>&>(vortonDistance_);

        is = dict.lookup("vortonScale");
        is >> static_cast<List<vector>&>(vortonScale_);

        is = dict.lookup("vortonIntensity");
        is >> static_cast<List<vector>&>(vortonIntensity_);
    }
}


Foam::turbulentATSMInletFvPatchVectorField::
turbulentATSMInletFvPatchVectorField
(
    const turbulentATSMInletFvPatchVectorField& ptf
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

    vortons_(),
    vortonType_(ptf.vortonType_),
    nVortonGlobal_(Zero),
    nVortonLocal_(Zero),
    vortonLabel_(),
    vortonPosition_(),
    vortonDistance_(),
    vortonScale_(),
    vortonIntensity_(),

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


Foam::turbulentATSMInletFvPatchVectorField::
turbulentATSMInletFvPatchVectorField
(
    const turbulentATSMInletFvPatchVectorField& ptf,
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

    vortons_(),
    vortonType_(ptf.vortonType_),
    nVortonGlobal_(Zero),
    nVortonLocal_(Zero),
    vortonLabel_(),
    vortonPosition_(),
    vortonDistance_(),
    vortonScale_(),
    vortonIntensity_(),

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

void Foam::turbulentATSMInletFvPatchVectorField::autoMap
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


void Foam::turbulentATSMInletFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<vector>::rmap(ptf, addr);

    const turbulentATSMInletFvPatchVectorField& ATSMptf =
        refCast<const turbulentATSMInletFvPatchVectorField>(ptf);

    R_.rmap(ATSMptf.R_, addr);
    L_.rmap(ATSMptf.L_, addr);
    U_.rmap(ATSMptf.U_, addr);

    // Clear interpolator
    mapperPtr_.clear();
}


void Foam::turbulentATSMInletFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (curTimeIndex_ == -1)
    {
        initialisePatch();
        initialiseParameters();
        initialiseVortonBox();
        initialiseVortons();
        initialiseOutput();
    }

    if (curTimeIndex_ != db().time().timeIndex())
    {
        if (debug)
        {
            label n = vortons_.size();
            Info<< "Number of vortons: " << returnReduce(n, sumOp<label>())
                << endl;
        }

        const scalar deltaT = db().time().deltaTValue();

        // Move vortons using mean velocity
        convectVortons(deltaT);

        // Set velocity
        vectorField& U = *this;
        //U = UMean_;
        U = U_*patchNormal_;

        const pointField& Cf = patch().Cf();

        // Apply second part of normalisation coefficient
        const scalar c = Foam::sqrt(v0_)/Foam::sqrt(scalar(nVortonGlobal_));

        // In parallel, need to collect all vortons that will interact with
        // local faces

        if (singleProc_ || !Pstream::parRun())
        {
            U += c*uDashVorton(vortons_, Cf);
        }
        else
        {
            // Process local vorton contributions
            U += c*uDashVorton(vortons_, Cf);

            // Add contributions from overlapping vortons
            List<List<vorton>> overlappingVortons(Pstream::nProcs());
            calcOverlappingProcVortons(overlappingVortons);

            forAll(overlappingVortons, procI)
            {
                const List<vorton>& vortons = overlappingVortons[procI];

                if (vortons.size())
                {
                    //Pout<< "Applying " << vortons.size()
                    //    << " vortons from processor " << procI << endl;

                    U += c*uDashVorton(vortons, Cf);
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


void Foam::turbulentATSMInletFvPatchVectorField::write(Ostream& os) const
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

    writeEntryIfDifferent<label>(os, "nVorton", 0, nVortonGlobal_);
    writeEntryIfDifferent<label>(os, "nVortonLocal", 0, nVortonLocal_);
    writeEntryIfDifferent<word>(os, "vortonType", "typeR", vortonType_);

    if (nVortonLocal_)
    {
        writeEntry(os, "vortonLabel", vortonLabel_);
        writeEntry(os, "vortonPosition", vortonPosition_);
        writeEntry(os, "vortonDistance", vortonDistance_);
        writeEntry(os, "vortonScale", vortonScale_);
        writeEntry(os, "vortonIntensity", vortonIntensity_);
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
       turbulentATSMInletFvPatchVectorField
   );
}


// ************************************************************************* //
