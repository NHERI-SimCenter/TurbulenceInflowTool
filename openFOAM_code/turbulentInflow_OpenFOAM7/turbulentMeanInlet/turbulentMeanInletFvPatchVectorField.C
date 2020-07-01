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

#include "turbulentMeanInletFvPatchVectorField.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "IFstream.H"
#include "OFstream.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::turbulentMeanInletFvPatchVectorField::initialise()
{
    const vectorField nf(patch().nf());

    // Patch normal points into domain
    patchNormal_ = -gAverage(nf);

    // Check that patch is planar
    scalar error = max(magSqr(patchNormal_+nf));

    if (error > SMALL)
    {
        WarningInFunction << "Patch " << patch().name() << " is not planar" << endl;
    }

    forAll (patchNormal_, label)
    {
        if (mag(patchNormal_[label]) <= SMALL)
        {
            patchNormal_[label] = 0;
        }
    }
}

const Foam::pointToPointPlanarInterpolation&
Foam::turbulentMeanInletFvPatchVectorField::patchMapper() const
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
                patch().Cf(),
                perturb_,
                nearestOnly
            )
        );
    }

    return *mapperPtr_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbulentMeanInletFvPatchVectorField::
turbulentMeanInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),

    isInitialized_(false),
    patchNormal_(vector::zero),

    perturb_(1e-5),
    mapMethod_("planarInterpolation"),
    mapperPtr_(nullptr),

    interpolateU_(false),
    U_(p.size(), 0.0)
{}


Foam::turbulentMeanInletFvPatchVectorField::
turbulentMeanInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),

    isInitialized_(false),
    patchNormal_(vector::zero),

    perturb_(dict.lookupOrDefault<scalar>("perturb", 1e-5)),
    mapMethod_(dict.lookupOrDefault<word>("mapMethod", "nearestCell")),
    mapperPtr_(nullptr),

    interpolateU_(false),
    U_(interpolateOrRead<scalar>("U", dict, interpolateU_))
{}


Foam::turbulentMeanInletFvPatchVectorField::
turbulentMeanInletFvPatchVectorField
(
    const turbulentMeanInletFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),

    isInitialized_(ptf.isInitialized_),
    patchNormal_(ptf.patchNormal_),

    perturb_(ptf.perturb_),
    mapMethod_(ptf.mapMethod_),
    mapperPtr_(nullptr),

    interpolateU_(ptf.interpolateU_),
    U_(mapper(ptf.U_))
{}


Foam::turbulentMeanInletFvPatchVectorField::
turbulentMeanInletFvPatchVectorField
(
    const turbulentMeanInletFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),

    isInitialized_(ptf.isInitialized_),
    patchNormal_(ptf.patchNormal_),

    perturb_(ptf.perturb_),
    mapMethod_(ptf.mapMethod_),
    mapperPtr_(nullptr),

    interpolateU_(ptf.interpolateU_),
    U_(ptf.U_)
{}


Foam::turbulentMeanInletFvPatchVectorField::
turbulentMeanInletFvPatchVectorField
(
    const turbulentMeanInletFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),

    isInitialized_(ptf.isInitialized_),
    patchNormal_(ptf.patchNormal_),
    
    perturb_(ptf.perturb_),
    mapMethod_(ptf.mapMethod_),
    mapperPtr_(nullptr),

    interpolateU_(ptf.interpolateU_),
    U_(ptf.U_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::turbulentMeanInletFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if(!isInitialized_)
    {
        initialise();
        isInitialized_ = true;
    }

    vectorField& U = *this;
    U = (U_*patchNormal_);
    
    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::turbulentMeanInletFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "value", *this);
    writeEntry(os, "U", U_);
}

void Foam::turbulentMeanInletFvPatchVectorField::autoMap(const fvPatchFieldMapper& mapper)
{
    fixedValueFvPatchVectorField::autoMap(mapper);
    mapper(U_, U_);
}

void Foam::turbulentMeanInletFvPatchVectorField::rmap(const fvPatchField<vector>& ptf, const labelList& addr)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    const turbulentMeanInletFvPatchVectorField& tiptf = refCast<const turbulentMeanInletFvPatchVectorField>(ptf);

    U_.rmap(tiptf.U_, addr);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        turbulentMeanInletFvPatchVectorField
    );
}

// ************************************************************************* //
