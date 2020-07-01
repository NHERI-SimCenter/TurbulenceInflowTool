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

#ifndef turbulentKineticEnergyFvPatchScalarField_H
#define turbulentKineticEnergyFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "Random.H"
#include <functional>
#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class pointToPointPlanarInterpolation;

/*---------------------------------------------------------------------------*\
        Class turbulentKineticEnergyFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class turbulentKineticEnergyFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private Data

        //- Current time index
        label curTimeIndex_;

        //- 2D interpolation (for 'planarInterpolation' mapMethod)
        mutable autoPtr<pointToPointPlanarInterpolation> mapperPtr_;

        //- Method for interpolation between a patch and turbulence plane
        const word mapMethod_;

        //- Fraction of perturbation (fraction of bounding box) to add (for
        const scalar perturb_;

        //- Flag: interpolate R field (default=false)
        bool interpolateR_;

        //- Reynold stress tensor
        symmTensorField R_;

        //- Kinetic energy
        scalarField k_;


    // Private Member Functions

        //- Return a reference to the patch mapper object
        const pointToPointPlanarInterpolation& patchMapper() const;

        //- Helper function to interpolate values from the boundary data or
        //- read from dictionary
        template<class Type>
        tmp<Field<Type>> interpolateOrRead
        (
            const word& fieldName,
            const dictionary& dict,
            bool& interpolateField
        ) const;

        //- Helper function to interpolate values from the boundary data
        template<class Type>
        tmp<Field<Type>> interpolateBoundaryData
        (
            const word& fieldName
        ) const;

        template<class Type>
        tmp<Field<Type>> calculateBoundaryData
        (
            const word& fieldName,
            const dictionary& dict
        ) const;


public:

   //- Runtime type information
   TypeName("turbulentKineticEnergy");

    // Constructors

        //- Construct from patch and internal field
        turbulentKineticEnergyFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentKineticEnergyFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- turbulentKineticEnergyFvPatchScalarField onto a new patch
        turbulentKineticEnergyFvPatchScalarField
        (
            const turbulentKineticEnergyFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        turbulentKineticEnergyFvPatchScalarField
        (
            const turbulentKineticEnergyFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new turbulentKineticEnergyFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        turbulentKineticEnergyFvPatchScalarField
        (
            const turbulentKineticEnergyFvPatchScalarField&,
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
                new turbulentKineticEnergyFvPatchScalarField(*this, iF)
            );
        }


    //- Destructor
    virtual ~turbulentKineticEnergyFvPatchScalarField() = default;


    // Member Functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


       //- Write
       virtual void write(Ostream&) const;


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            //  Used to update fields following mesh topology change
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            //  Used to reconstruct fields
            virtual void rmap(const fvPatchField<scalar>&, const labelList&);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "turbulentKineticEnergyFvPatchScalarFieldTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
