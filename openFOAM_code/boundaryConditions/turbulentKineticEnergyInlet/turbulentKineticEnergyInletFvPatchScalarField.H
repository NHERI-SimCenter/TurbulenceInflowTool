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

Class
    Foam::turbulentKineticEnergyInletFvPatchScalarField

Description
    This boundary condition provides a turbulent kinetic energy condition,
    based on user-supplied turbulence intensity, defined as a fraction of the
    mean velocity:

        \f[
            k_p = 1.5 (I |U|)^2
        \f]

    where

    \vartable
        k_p     | kinetic energy at the patch
        I       | turbulence intensity
        U       | velocity field
    \endvartable

    In the event of reverse flow, a zero-gradient condition is applied.

Usage
    \table
        Property     | Description             | Required    | Default value
        intensity    | fraction of mean field [0-1] | yes    |
        U            | velocity field name     | no          | U
        phi          | flux field name         | no          | phi
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type        turbulentKineticEnergyInlet;
        intensity   0.05;           // 5% turbulence
        value       uniform 1;      // placeholder
    }
    \endverbatim

See also
    Foam::inletOutletFvPatchField

SourceFiles
    turbulentKineticEnergyInletFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentKineticEnergyInletFvPatchScalarField_H
#define turbulentKineticEnergyInletFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
   Class turbulentKineticEnergyInletFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class turbulentKineticEnergyInletFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Name of the velocity field
        word UName_;

    // Private Member Functions

        //- Helper function to calculate values from the boundary data or
        //  read from dictionary
        template<class Type>
        tmp<Field<Type>> calculateOrRead
        (
            const word& fieldName,
            const dictionary& dict,
            bool& calculateField
        ) const;

        //- Helper function to calculate values from the boundary data
        template<class Type>
        tmp<Field<Type>> calculateBoundaryData
        (
            const word& fieldName,
            const dictionary& dict
        ) const;

public:

    //- Runtime type information
    TypeName("turbulentKineticEnergyInlet");


    // Constructors

        //- Construct from patch and internal field
        turbulentKineticEnergyInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentKineticEnergyInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  turbulentKineticEnergyInletFvPatchScalarField
        //  onto a new patch
        turbulentKineticEnergyInletFvPatchScalarField
        (
            const turbulentKineticEnergyInletFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        turbulentKineticEnergyInletFvPatchScalarField
        (
            const turbulentKineticEnergyInletFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new turbulentKineticEnergyInletFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        turbulentKineticEnergyInletFvPatchScalarField
        (
            const turbulentKineticEnergyInletFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new turbulentKineticEnergyInletFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "turbulentKineticEnergyInletFvPatchScalarFieldTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
