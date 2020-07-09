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

#include "turbulentKineticEnergyFvPatchScalarField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::pointToPointPlanarInterpolation&
Foam::turbulentKineticEnergyFvPatchScalarField::patchMapper() const
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbulentKineticEnergyFvPatchScalarField::
turbulentKineticEnergyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    curTimeIndex_(-1),
    mapperPtr_(nullptr),
    mapMethod_("nearestCell"),
    perturb_(1e-5),
    interpolateR_(false),
    R_(p.size(), symmTensor::zero),
    k_(p.size(), 0)
{}


Foam::turbulentKineticEnergyFvPatchScalarField::
turbulentKineticEnergyFvPatchScalarField
(
    const turbulentKineticEnergyFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    curTimeIndex_(-1),
    mapperPtr_(nullptr),
    mapMethod_(ptf.mapMethod_),
    perturb_(ptf.perturb_),
    interpolateR_(false),
    R_(mapper(ptf.R_)),
    k_(mapper(ptf.k_))
{}


Foam::turbulentKineticEnergyFvPatchScalarField::
turbulentKineticEnergyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    curTimeIndex_(-1),
    mapperPtr_(nullptr),
    mapMethod_(dict.lookupOrDefault<word>("mapMethod", "nearestCell")),
    perturb_(dict.lookupOrDefault<scalar>("perturb", 1e-5)),
    interpolateR_(false),
    R_(interpolateOrRead<symmTensor>("R", dict, interpolateR_)),
    k_(0.5*(R_.component(symmTensor::XX)+R_.component(symmTensor::YY)+R_.component(symmTensor::ZZ)))
{}


Foam::turbulentKineticEnergyFvPatchScalarField::
turbulentKineticEnergyFvPatchScalarField
(
    const turbulentKineticEnergyFvPatchScalarField& ptf
)
:
    fixedValueFvPatchField<scalar>(ptf),
    curTimeIndex_(-1),
    mapperPtr_(nullptr),
    mapMethod_(ptf.mapMethod_),
    perturb_(ptf.perturb_),
    interpolateR_(false),
    R_(ptf.R_),
    k_(ptf.k_)
{}


Foam::turbulentKineticEnergyFvPatchScalarField::
turbulentKineticEnergyFvPatchScalarField
(
    const turbulentKineticEnergyFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ptf, iF),
    curTimeIndex_(-1),
    mapperPtr_(nullptr),
    mapMethod_(ptf.mapMethod_),
    perturb_(ptf.perturb_),
    interpolateR_(false),
    R_(ptf.R_),
    k_(ptf.k_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::turbulentKineticEnergyFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (curTimeIndex_ != db().time().timeIndex())
    {
        scalarField& k = *this;

        k = k_;

        curTimeIndex_ = db().time().timeIndex();
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::turbulentKineticEnergyFvPatchScalarField::write
(
    Ostream& os
) const
{
    fixedValueFvPatchField<scalar>::write(os);

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
    
    writeEntry(os, "R", R_);
}

void Foam::turbulentKineticEnergyFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    m(R_,R_);
}


void Foam::turbulentKineticEnergyFvPatchScalarField::rmap
(
    const fvPatchField<scalar>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const turbulentKineticEnergyFvPatchScalarField& tkptf =
        refCast<const turbulentKineticEnergyFvPatchScalarField>(ptf);

    R_.rmap(tkptf.R_, addr);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchScalarField,
       turbulentKineticEnergyFvPatchScalarField
   );
}


// ************************************************************************* //
