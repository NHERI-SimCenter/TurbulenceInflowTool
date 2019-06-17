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

#include "turbulentDFInletFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::pointToPointPlanarInterpolation&
Foam::turbulentDFInletFvPatchVectorField::patchMapper() const
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

void Foam::turbulentDFInletFvPatchVectorField::initialisePatch()
{
    const vectorField nf(patch().nf());

    // Patch normal points into domain
    vector patchNormal = -gAverage(nf);

    // Check that patch is planar
    scalar error = max(magSqr(patchNormal + nf));

    if (error > SMALL)
    {
        WarningInFunction << "Patch " << patch().name() << " is not planar" << endl;
    }
     
    // Set number of patch faces for each processor
    patchSize_[Pstream::myProcNo()] = patch().Cf().size();

    Pstream::gatherList(patchSize_);
    Pstream::scatterList(patchSize_);
}

void Foam::turbulentDFInletFvPatchVectorField::initialiseVirtualGrid()
{
    // Set the origin and spacing of the virtual grid points
    boundBox patchBounds(patch().Cf());

    origin_ = patchBounds.min();

    delta_ = Foam::sqrt(gMin(patch().magSf()))*gridFactor_;

    My_ = patchBounds.span().component(vector::Y)/delta_ + 1;
    Mz_ = patchBounds.span().component(vector::Z)/delta_ + 1;

    if (debug)
    {
        Pout<< origin_ << " " << delta_ << " " << My_ << " " << Mz_ << endl;
    }

    // Get the patch face centres from all processors
    vectorField fCentres = gatherProc(patch().Cf());

    // Sort the face centres and set the patch faces to virtual girds mapping
    facesToIndices_.setSize(fCentres.size(), -1);
    sortfaceCentres(fCentres, facesToIndices_);

    vectorField fCentresSorted(fCentres, facesToIndices_);

    // Set the virtual girds to patch faces mapping
    labelList inverseMapping(facesToIndices_);

    sortedOrder(facesToIndices_, inverseMapping);

    label start = 0;

    if (Pstream::myProcNo() > 0)
    {
        for (label i = 0; i < Pstream::myProcNo(); i++)
        {
            start = start + patchSize_[i];
        }
    }

    if (patchSize_[Pstream::myProcNo()] > 0)
    {
        indicesToFaces_ = SubList<label>(inverseMapping, patchSize_[Pstream::myProcNo()], start);
    }

    vectorField ff(fCentresSorted, indicesToFaces_);

    // Set the virtual grid points size
    virualGridPoints_.setSize(fCentres.size(), vector::zero);

    // Set the coordinate of virtual grid points
    forAll(fCentresSorted, faceI)
    {
        scalar yCord = 0;
        scalar zCord = 0;

        scalar distY = fCentresSorted[faceI].component(vector::Y)-origin_.component(vector::Y);
        scalar distZ = fCentresSorted[faceI].component(vector::Z)-origin_.component(vector::Z);

        scalar remY = distY - floor(distY/delta_) * delta_;
        scalar remZ = distZ - floor(distZ/delta_) * delta_;

        if (remY > delta_/2) 
        {
            yCord = origin_.component(vector::Y) + distY - remY + delta_;
        }
        else 
        {
            yCord = origin_.component(vector::Y) + distY - remY;
        }

        if (remZ > delta_/2) 
        {
            zCord = origin_.component(vector::Z) + distZ - remZ + delta_;
        }
        else 
        {
            zCord = origin_.component(vector::Z) + distZ - remZ;
        }

        virualGridPoints_[faceI].component(vector::X) = origin_.component(vector::X);
        virualGridPoints_[faceI].component(vector::Y) = yCord;
        virualGridPoints_[faceI].component(vector::Z) = zCord;
    }

    yindices_.setSize(virualGridPoints_.size());
    zindices_.setSize(virualGridPoints_.size());

    forAll(virualGridPoints_, I)
    {
        yindices_[I] = round((virualGridPoints_[I].component(vector::Y)-origin_.component(vector::Y))/delta_);
        zindices_[I] = round((virualGridPoints_[I].component(vector::Z)-origin_.component(vector::Z))/delta_);
    }

    // Set the number of grids for each processsor
    indicesPerProc_ = 0;
    rest_ = 0;

    if (Pstream::master())
    {
        indicesPerProc_ = floor(scalar(virualGridPoints_.size())/(Pstream::nProcs()));
        rest_ = virualGridPoints_.size()-(indicesPerProc_*(Pstream::nProcs()));

        Info << nl << "Generating Inflow for " << virualGridPoints_.size() << " faces" << endl;
        Info << "Distributing Indices per Proc: " << indicesPerProc_ << endl;

        if (rest_ > 0)
        {
            Info << "First " << rest_ << " Procs will do +1 Indices" << nl << endl;
        }
    }

    Pstream::scatter(indicesPerProc_);
    Pstream::scatter(rest_);
}

void Foam::turbulentDFInletFvPatchVectorField::initialiseDigitalFilter()
{
    // Gather the length scale field from each processors
    vectorField Lu(gatherProc(Lu_), facesToIndices_);
    vectorField Lv(gatherProc(Lv_), facesToIndices_);
    vectorField Lw(gatherProc(Lw_), facesToIndices_);

    ny_.setSize(Lu.size());
    nz_.setSize(Lu.size());

    forAll (Lu, label)
    {
        ny_[label] = vector(round(Lu[label].component(1)/delta_),round(Lv[label].component(1)/delta_),round(Lw[label].component(1)/delta_));
        nz_[label] = vector(round(Lu[label].component(2)/delta_),round(Lv[label].component(2)/delta_),round(Lw[label].component(2)/delta_));
    }

    nyMax_ = gMax(ny_);
    nzMax_ = gMax(nz_);

    // Set size of virtual grid
    rndSize_.component(0) = (My_+nfK_*nyMax_.component(0))*(Mz_+nfK_*nzMax_.component(0));
    rndSize_.component(1) = (My_+nfK_*nyMax_.component(1))*(Mz_+nfK_*nzMax_.component(1));
    rndSize_.component(2) = (My_+nfK_*nyMax_.component(2))*(Mz_+nfK_*nzMax_.component(2));

    //For current processor, get start and end index in virtual list.
    label start;
    label size;

    if (Pstream::myProcNo() < rest_)
    {
        start = Pstream::myProcNo()*indicesPerProc_ + Pstream::myProcNo();
        size = indicesPerProc_+1;
    }
    else
    {
        start = Pstream::myProcNo()*indicesPerProc_ + rest_;
        size = indicesPerProc_;
    }

    //only fill element of current proc
    filterCoeffProc_u.setSize(size);
    filterCoeffProc_v.setSize(size);
    filterCoeffProc_w.setSize(size);

    //loop through all indices in one proc, and set filter kernel size and fill with 0.0
    for (label subI = 0; subI < size; subI++)
    {
        //index in full array
        label I = subI+start;

        filterCoeffProc_u[subI].setSize((nfK_*ny_[I].component(0)+1)*(nfK_*nz_[I].component(0)+1), 0.0);
        filterCoeffProc_v[subI].setSize((nfK_*ny_[I].component(1)+1)*(nfK_*nz_[I].component(1)+1), 0.0);
        filterCoeffProc_w[subI].setSize((nfK_*ny_[I].component(2)+1)*(nfK_*nz_[I].component(2)+1), 0.0);

        get2DFilterCoeff(filterCoeffProc_u[subI], ny_[I].component(0), nz_[I].component(0));
        get2DFilterCoeff(filterCoeffProc_v[subI], ny_[I].component(1), nz_[I].component(1));
        get2DFilterCoeff(filterCoeffProc_w[subI], ny_[I].component(2), nz_[I].component(2));
    }
}

void Foam::turbulentDFInletFvPatchVectorField::sortfaceCentres
(
    const vectorField& fCentres_,
    labelList& order_
)
{
    List<scalar> yCoordinate;
    yCoordinate.setSize(fCentres_.size());

    forAll(yCoordinate, i)
    {
        yCoordinate[i] = fCentres_[i].component(vector::Y);
    }

    List<label> yOrder;
    sortedOrder(yCoordinate, yOrder);

    sort(yCoordinate);

    List<label> unique;
    uniqueOrder(yCoordinate, unique);

    for (int i = 0; i < unique.size(); i++)
    {
        label size, start;

        if (i == 0)
        {
            size = unique[i] + 1;
            start = 0;
        }
        else
        {
            size = unique[i]-unique[i-1];
            start = unique[i-1] + 1;
        }

        List<scalar> zCoordinate(size, 0.0);

        for (int j = 0; j < size; j++)
        {
            zCoordinate[j] = fCentres_[yOrder[start+j]].component(vector::Z);
        }

        List<label> zOrder;
        sortedOrder(zCoordinate, zOrder);

        for (int j = 0; j < size; j++)
        {
            order_[start+j] = yOrder[start+zOrder[j]];
        }
    }
}


inline Foam::label 
Foam::turbulentDFInletFvPatchVectorField::get1DIndex(label x, label y, label yDim)
{
    return x * yDim + y;
}


void Foam::turbulentDFInletFvPatchVectorField::get1DFilterCoeff(scalarList& b, label n)
{
    const scalar pi = constant::mathematical::pi;

    if (n == 0)
    {
        b.setSize(1, 0.0);
    }
    else
    {
        scalar sum = 0.0;

        label N = nfK_*n/2;

        b.setSize(2*N+1);

        for (label k = 0; k < 2*N+1; k++)
        {
            if (filterShape_ == "exponential")
            {
                b[k] = Foam::exp(-fabs(1.0*(k-N)/n));
            }
            else if (filterShape_ == "gaussian")
            {
                b[k] = Foam::exp(-pi*Foam::sqr(scalar(k-N))/(2.0*Foam::sqr(scalar(n))));
            }
            else
            {
                Info << "filter coefficient function" << filterShape_ << " does not exist (ERROR)" << endl;
            }

            sum += b[k]*b[k];
        }

        sum = Foam::sqrt(sum);
        b = b/sum;
    }
}

void Foam::turbulentDFInletFvPatchVectorField::get2DFilterCoeff(scalarList& filter,label ny, label nz)
{
    scalarList by;
    scalarList bz;

    get1DFilterCoeff(by, ny);
    get1DFilterCoeff(bz, nz);

    for (label i = 0; i < nfK_*ny+1; i++)
    {
        for (label j = 0; j < nfK_*nz+1; j++)
        {
            filter[get1DIndex(i, j, nfK_*nz+1)] = by[i]*bz[j];
        }
    }
}

Foam::scalarField
Foam::turbulentDFInletFvPatchVectorField::getRandomField(label totalSize)
{
    List<scalarField> virtualRandomFieldProc;

    virtualRandomFieldProc.setSize(Pstream::nProcs());

    label size = floor(totalSize/Pstream::nProcs());
    label rest = totalSize-size*Pstream::nProcs();

    if (Pstream::myProcNo() < rest)
    {
        virtualRandomFieldProc[Pstream::myProcNo()].setSize(size+1);
    }
    else
    {
        virtualRandomFieldProc[Pstream::myProcNo()].setSize(size);
    }

    forAll (virtualRandomFieldProc[Pstream::myProcNo()], i)
    {
        virtualRandomFieldProc[Pstream::myProcNo()][i] = rndGen_.scalarNormal();
    }

    Pstream::gatherList(virtualRandomFieldProc);
    Pstream::scatterList(virtualRandomFieldProc);

    return ListListOps::combine<scalarField>(virtualRandomFieldProc, accessOp<scalarField>());
}

void Foam::turbulentDFInletFvPatchVectorField::spatialCorr()
{
    scalarField virtualRandomField_u = getRandomField(rndSize_.component(0));
    scalarField virtualRandomField_v = getRandomField(rndSize_.component(1));
    scalarField virtualRandomField_w = getRandomField(rndSize_.component(2));

    label start, size;

    if (Pstream::myProcNo() < rest_)
    {
        start = Pstream::myProcNo()*indicesPerProc_ + Pstream::myProcNo();
        size = indicesPerProc_ + 1;
    }
    else
    {
        start = Pstream::myProcNo()*indicesPerProc_ + rest_;
        size = indicesPerProc_;
    }
    
    List<vectorField> virtualFilteredFieldProc;

    virtualFilteredFieldProc.setSize(Pstream::nProcs());
    virtualFilteredFieldProc[Pstream::myProcNo()].setSize(size, vector::zero);

    labelVector yOffset = nfK_*nyMax_/2;
    labelVector zOffset = nfK_*nzMax_/2;

    forAll (virtualFilteredFieldProc[Pstream::myProcNo()], subI)
    {
        label I = subI+start;

        label i = yindices_[I]; // i = yindices on virtual Grid
        label j = zindices_[I]; // j = zindices on virtual Grid

        vector u = vector::zero;

        for (label ii = 0; ii < nfK_*ny_[I].component(0)+1; ii++)
        {
            label start_rnd = get1DIndex
            (
                i+yOffset.component(0)-nfK_/2*ny_[I].component(0)+ii,
                j+zOffset.component(0)-nfK_/2*nz_[I].component(0),
                Mz_+nfK_*nzMax_.component(0)
            );

            label size_rnd = nfK_*nz_[I].component(0)+1;
            label start_filt = get1DIndex(ii, 0, nfK_*nz_[I].component(0)+1);

            scalarField rnd = SubField<scalar>(virtualRandomField_u, size_rnd, start_rnd);
            scalarField filt = SubField<scalar>(filterCoeffProc_u[subI], size_rnd, start_filt);

            u.component(0) += sumProd(rnd,filt);
        }

        for (label ii = 0; ii < nfK_*ny_[I].component(1)+1; ii++)
        {
            label start_rnd = get1DIndex
            (
                i+yOffset.component(1)-nfK_/2*ny_[I].component(1)+ii,
                j+zOffset.component(1)-nfK_/2*nz_[I].component(1),
                Mz_+nfK_*nzMax_.component(1)
            );

            label size_rnd = nfK_*nz_[I].component(1)+1;
            label start_filt = get1DIndex(ii, 0, nfK_*nz_[I].component(1)+1);

            scalarField rnd = SubField<scalar>(virtualRandomField_v, size_rnd, start_rnd);
            scalarField filt = SubField<scalar>(filterCoeffProc_v[subI], size_rnd, start_filt);

            u.component(1) += sumProd(rnd,filt);
        }

        for (label ii = 0; ii < nfK_*ny_[I].component(2)+1; ii++)
        {
            label start_rnd = get1DIndex
            (
                i+yOffset.component(2)-nfK_/2*ny_[I].component(2)+ii,
                j+zOffset.component(2)-nfK_/2*nz_[I].component(2),
                Mz_+nfK_*nzMax_.component(2)
            );

            label size_rnd = nfK_*nz_[I].component(2)+1;
            label start_filt = get1DIndex(ii, 0, nfK_*nz_[I].component(2)+1);

            scalarField rnd = SubField<scalar>(virtualRandomField_w, size_rnd, start_rnd);
            scalarField filt = SubField<scalar>(filterCoeffProc_w[subI], size_rnd, start_filt);

            u.component(2) += sumProd(rnd,filt);
        }

        virtualFilteredFieldProc[Pstream::myProcNo()][subI] = u;
    }

    Pstream::gatherList(virtualFilteredFieldProc);
    Pstream::scatterList(virtualFilteredFieldProc);

    vectorField virtualFilteredField = ListListOps::combine<vectorField>(virtualFilteredFieldProc, accessOp<vectorField>());

    uFluctFiltered_ = vectorField(virtualFilteredField, indicesToFaces_);
}

void Foam::turbulentDFInletFvPatchVectorField::temporalCorr()
{
    scalar dt = db().time().deltaT().value();

    forAll(uFluctTemporal_, faceI)
    {
        vector T = vector(Lu_[faceI].component(0),Lv_[faceI].component(0),Lw_[faceI].component(0))/U_[faceI].component(0);

        uFluctTemporal_[faceI].component(0) = uFluctTemporalOld_[faceI].component(0)*Foam::exp(-dt/T.component(0))
                                            + uFluctFiltered_[faceI].component(0)*Foam::sqrt(1.0-Foam::exp(-2.0*dt/T.component(0)));

        uFluctTemporal_[faceI].component(1) = uFluctTemporalOld_[faceI].component(1)*Foam::exp(-dt/T.component(1))
                                            + uFluctFiltered_[faceI].component(1)*Foam::sqrt(1.0-Foam::exp(-2.0*dt/T.component(1)));

        uFluctTemporal_[faceI].component(2) = uFluctTemporalOld_[faceI].component(2)*Foam::exp(-dt/T.component(2))
                                            + uFluctFiltered_[faceI].component(2)*Foam::sqrt(1.0-Foam::exp(-2.0*dt/T.component(2)));
    }

    uFluctTemporalOld_ = uFluctTemporal_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbulentDFInletFvPatchVectorField::
turbulentDFInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),

    curTimeIndex_(-1),
    patchSize_(Pstream::nProcs(), 0),
    uFluctFiltered_(p.size(), pTraits<vector>::zero),
    uFluctTemporalOld_(p.size(), pTraits<vector>::zero),
    uFluctTemporal_(p.size(), pTraits<vector>::zero),

    perturb_(1e-5),
    mapMethod_("planarInterpolation"),
    mapperPtr_(nullptr),
    interpolateR_(false),
    R_(p.size(), pTraits<symmTensor>::zero),
    interpolateL_(false),
    Lu_(p.size(), pTraits<vector>::zero),
    Lv_(p.size(), pTraits<vector>::zero),
    Lw_(p.size(), pTraits<vector>::zero),
    interpolateU_(false),
    U_(p.size(), pTraits<vector>::zero),
    Lund_(p.size(), pTraits<tensor>::zero),

    isInitialized_(false),
    gridFactor_(1.0),
    origin_(vector::zero),
    My_(0),
    Mz_(0),
    delta_(0),
    ny_(),
    nz_(),
    nyMax_(vector::zero),
    nzMax_(vector::zero),
    nfK_(4),
    yindices_(),
    zindices_(),

    indicesPerProc_(0),
    rest_(0),
    facesToIndices_(),
    indicesToFaces_(),

    rndGen_((Pstream::myProcNo()+1)*time(NULL)),
    filterShape_("exponential"),
    virualGridPoints_(),
    rndSize_(vector::zero),
    filterCoeffProc_u(),
    filterCoeffProc_v(),
    filterCoeffProc_w()
{}


Foam::turbulentDFInletFvPatchVectorField::
turbulentDFInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),

    curTimeIndex_(-1),
    patchSize_(Pstream::nProcs(), 0),
    uFluctFiltered_(p.size(), pTraits<vector>::zero),
    uFluctTemporalOld_(p.size(), pTraits<vector>::zero),
    uFluctTemporal_(p.size(), pTraits<vector>::zero),

    perturb_(dict.lookupOrDefault<scalar>("perturb", 1e-5)),
    mapMethod_(dict.lookup("mapMethod")),
    mapperPtr_(nullptr),
    interpolateR_(false),
    R_(interpolateOrRead<symmTensor>("R", dict, interpolateR_)),
    interpolateL_(false),
    Lu_(interpolateOrRead<vector>("Lu", dict, interpolateL_)),
    Lv_(interpolateOrRead<vector>("Lv", dict, interpolateL_)),
    Lw_(interpolateOrRead<vector>("Lw", dict, interpolateL_)),
    interpolateU_(false),
    U_(interpolateOrRead<vector>("U", dict, interpolateU_)),
    Lund_(p.size(), pTraits<tensor>::zero),

    isInitialized_(false),
    gridFactor_(dict.lookupOrDefault("gridFactor", 1.0)),
    origin_(vector::zero),
    My_(0),
    Mz_(0),
    delta_(0),
    ny_(),
    nz_(),
    nyMax_(vector::zero),
    nzMax_(vector::zero),
    nfK_(dict.lookupOrDefault<label>("filterFactor", 4)),
    yindices_(),
    zindices_(),

    indicesPerProc_(0),
    rest_(0),
    facesToIndices_(),
    indicesToFaces_(),

    rndGen_((Pstream::myProcNo()+1)*time(NULL)),
    filterShape_(dict.lookupOrDefault<word>("filterShape", "exponential")),
    virualGridPoints_(),
    rndSize_(vector::zero),
    filterCoeffProc_u(),
    filterCoeffProc_v(),
    filterCoeffProc_w()
{
    Lund_.replace(tensor::XX, sqrt(R_.component(symmTensor::XX)));
    Lund_.replace(tensor::YX, R_.component(symmTensor::XY)/Lund_.component(tensor::XX));
    Lund_.replace(tensor::ZX, R_.component(symmTensor::XZ)/Lund_.component(tensor::XX));
    Lund_.replace(tensor::YY, sqrt(R_.component(symmTensor::YY)-sqr(Lund_.component(tensor::YX))));
    Lund_.replace(tensor::ZY, (R_.component(symmTensor::YZ) - Lund_.component(tensor::YX)*Lund_.component(tensor::ZX) )/Lund_.component(tensor::YY));
    Lund_.replace(tensor::ZZ, sqrt(R_.component(symmTensor::ZZ) - sqr(Lund_.component(tensor::ZX))-sqr(Lund_.component(tensor::ZY))));
}


Foam::turbulentDFInletFvPatchVectorField::
turbulentDFInletFvPatchVectorField
(
    const turbulentDFInletFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),

    curTimeIndex_(ptf.curTimeIndex_),
    patchSize_(ptf.patchSize_),
    uFluctFiltered_(ptf.uFluctFiltered_, mapper),
    uFluctTemporalOld_(ptf.uFluctTemporalOld_, mapper),
    uFluctTemporal_(ptf.uFluctTemporal_, mapper),

    perturb_(ptf.perturb_),
    mapMethod_(ptf.mapMethod_),
    mapperPtr_(nullptr),
    interpolateR_(ptf.interpolateR_),
    R_(ptf.R_, mapper),
    interpolateL_(ptf.interpolateL_),
    Lu_(ptf.Lu_, mapper),
    Lv_(ptf.Lv_, mapper),
    Lw_(ptf.Lw_, mapper),
    interpolateU_(ptf.interpolateU_),
    U_(ptf.U_, mapper),
    Lund_(ptf.Lund_, mapper),

    isInitialized_(ptf.isInitialized_),
    gridFactor_(ptf.gridFactor_),
    origin_(ptf.origin_),
    My_(ptf.My_),
    Mz_(ptf.Mz_),
    delta_(ptf.delta_),
    ny_(ptf.ny_),
    nz_(ptf.nz_),
    nyMax_(ptf.nyMax_),
    nzMax_(ptf.nzMax_),
    nfK_(ptf.nfK_),
    yindices_(ptf.yindices_),
    zindices_(ptf.zindices_),

    indicesPerProc_(ptf.indicesPerProc_),
    rest_(ptf.rest_),
    facesToIndices_(ptf.facesToIndices_),
    indicesToFaces_(ptf.indicesToFaces_),

    rndGen_(ptf.rndGen_),
    filterShape_(ptf.filterShape_),
    virualGridPoints_(),
    rndSize_(ptf.rndSize_),
    filterCoeffProc_u(),
    filterCoeffProc_v(),
    filterCoeffProc_w()
{}


Foam::turbulentDFInletFvPatchVectorField::
turbulentDFInletFvPatchVectorField
(
    const turbulentDFInletFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),

    curTimeIndex_(ptf.curTimeIndex_),
    patchSize_(ptf.patchSize_),
    uFluctFiltered_(ptf.uFluctFiltered_),
    uFluctTemporalOld_(ptf.uFluctTemporalOld_),
    uFluctTemporal_(ptf.uFluctTemporal_),

    perturb_(ptf.perturb_),
    mapMethod_(ptf.mapMethod_),
    mapperPtr_(nullptr),
    interpolateR_(ptf.interpolateR_),
    R_(ptf.R_),
    interpolateL_(ptf.interpolateL_),
    Lu_(ptf.Lu_),
    Lv_(ptf.Lv_),
    Lw_(ptf.Lw_),
    interpolateU_(ptf.interpolateU_),
    U_(ptf.U_),
    Lund_(ptf.Lund_),

    isInitialized_(ptf.isInitialized_),
    gridFactor_(ptf.gridFactor_),
    origin_(ptf.origin_),
    My_(ptf.My_),
    Mz_(ptf.Mz_),
    delta_(ptf.delta_),
    ny_(ptf.ny_),
    nz_(ptf.nz_),
    nyMax_(ptf.nyMax_),
    nzMax_(ptf.nzMax_),
    nfK_(ptf.nfK_),
    yindices_(ptf.yindices_),
    zindices_(ptf.zindices_),

    indicesPerProc_(ptf.indicesPerProc_),
    rest_(ptf.rest_),
    facesToIndices_(ptf.facesToIndices_),
    indicesToFaces_(ptf.indicesToFaces_),

    rndGen_(ptf.rndGen_),
    filterShape_(ptf.filterShape_),
    virualGridPoints_(),
    rndSize_(ptf.rndSize_),
    filterCoeffProc_u(),
    filterCoeffProc_v(),
    filterCoeffProc_w()
{}


Foam::turbulentDFInletFvPatchVectorField::
turbulentDFInletFvPatchVectorField
(
    const turbulentDFInletFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),

    curTimeIndex_(ptf.curTimeIndex_),
    patchSize_(ptf.patchSize_),
    uFluctFiltered_(ptf.uFluctFiltered_),
    uFluctTemporalOld_(ptf.uFluctTemporalOld_),
    uFluctTemporal_(ptf.uFluctTemporal_),

    perturb_(ptf.perturb_),
    mapMethod_(ptf.mapMethod_),
    mapperPtr_(nullptr),
    interpolateR_(ptf.interpolateR_),
    R_(ptf.R_),
    interpolateL_(ptf.interpolateL_),
    Lu_(ptf.Lu_),
    Lv_(ptf.Lv_),
    Lw_(ptf.Lw_),
    interpolateU_(ptf.interpolateU_),
    U_(ptf.U_),
    Lund_(ptf.Lund_),

    isInitialized_(ptf.isInitialized_),
    gridFactor_(ptf.gridFactor_),
    origin_(ptf.origin_),
    My_(ptf.My_),
    Mz_(ptf.Mz_),
    delta_(ptf.delta_),
    ny_(ptf.ny_),
    nz_(ptf.nz_),
    nyMax_(ptf.nyMax_),
    nzMax_(ptf.nzMax_),
    nfK_(ptf.nfK_),
    yindices_(ptf.yindices_),
    zindices_(ptf.zindices_),

    indicesPerProc_(ptf.indicesPerProc_),
    rest_(ptf.rest_),
    facesToIndices_(ptf.facesToIndices_),
    indicesToFaces_(ptf.indicesToFaces_),

    rndGen_(ptf.rndGen_),
    filterShape_(ptf.filterShape_),
    virualGridPoints_(),
    rndSize_(ptf.rndSize_),
    filterCoeffProc_u(),
    filterCoeffProc_v(),
    filterCoeffProc_w()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::turbulentDFInletFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    //initialize virtual grid points
    if(!isInitialized_)
    {
        initialisePatch();
        initialiseVirtualGrid();
        initialiseDigitalFilter();

        spatialCorr();
        uFluctTemporalOld_ = uFluctFiltered_;

        isInitialized_=true;
    }



    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        curTimeIndex_ = db().time().timeIndex();

        // Filter random field
        spatialCorr();

        // create new temporally correlated slice
        temporalCorr();

        // Set final velocity field
        vectorField& U = *this;
        U = U_ + (Lund_ & uFluctTemporal_);

        // Re-scale to ensure correct flow rate
        scalar fCorr = gSum(U_ & patch().Sf())/gSum(U & patch().Sf());

        if (Pstream::master())
        {
            Info<< "mass flow correction coefficient: " << fCorr << endl;
        }

        U *= fCorr;

    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::turbulentDFInletFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry("value", os);

    os.writeKeyword("perturb") << perturb_ << token::END_STATEMENT << nl;

    os.writeKeyword("gridFactor") << gridFactor_ << token::END_STATEMENT << nl;
    os.writeKeyword("filterShape") << filterShape_ << token::END_STATEMENT << nl;
    os.writeKeyword("filterFactor") << nfK_ << token::END_STATEMENT << nl;

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
        turbulentDFInletFvPatchVectorField
    );
}

// ************************************************************************* //
