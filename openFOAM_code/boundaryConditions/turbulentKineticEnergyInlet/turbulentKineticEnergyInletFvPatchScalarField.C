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

#include "turbulentKineticEnergyInletFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbulentKineticEnergyInletFvPatchScalarField::
turbulentKineticEnergyInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    UName_("U")
{}

Foam::turbulentKineticEnergyInletFvPatchScalarField::
turbulentKineticEnergyInletFvPatchScalarField
(
    const turbulentKineticEnergyInletFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_)
{}

Foam::turbulentKineticEnergyInletFvPatchScalarField::
turbulentKineticEnergyInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    UName_(dict.lookupOrDefault<word>("U", "U"))
{
    const fvMesh& mesh = internalField().mesh();

    IOdictionary inflowDict
    (
        IOobject
        (
            "inflowProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    );

    bool calculateU_(false);
    bool calculateI_(false);

    vectorField U_(calculateOrRead<vector>("U", inflowDict, calculateU_));
    symmTensorField I_(calculateOrRead<symmTensor>("I", inflowDict, calculateI_));

    scalarField value(p.size(), 0.0);

    forAll (value, label)
    {
        value[label] = sqr(I_[label].component(symmTensor::XX))+sqr(I_[label].component(symmTensor::YY))+sqr(I_[label].component(symmTensor::ZZ));
    }  

    forAll (value, label)
    {
        value[label] = 0.5*sqr(mag(U_[label])*value[label]);
    }

    fvPatchScalarField::operator=(value);
}

Foam::turbulentKineticEnergyInletFvPatchScalarField::
turbulentKineticEnergyInletFvPatchScalarField
(
    const turbulentKineticEnergyInletFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    UName_(ptf.UName_)
{}


Foam::turbulentKineticEnergyInletFvPatchScalarField::
turbulentKineticEnergyInletFvPatchScalarField
(
    const turbulentKineticEnergyInletFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    UName_(ptf.UName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::turbulentKineticEnergyInletFvPatchScalarField::
updateCoeffs()
{
    if (updated())
    {
        return;
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::turbulentKineticEnergyInletFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        turbulentKineticEnergyInletFvPatchScalarField
    );
}

// ************************************************************************* //
