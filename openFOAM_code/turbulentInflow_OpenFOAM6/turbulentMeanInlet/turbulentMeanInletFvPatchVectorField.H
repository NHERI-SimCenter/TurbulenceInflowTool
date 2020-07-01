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
    Foam::turbulentMeanInletFvPatchVectorField

Description
    Turbulent velocity boundary condition based on the Digital Filter Method.

    Reference: 
        M. Klein, A. Sadiki, J. Janicka,
        "A ditial filter based generation of inflow data for spatially
        developing direct numerical or large eddy simulations",
        Journal of Computational Physics, (2003) 186:652:665

        Zheng-Tong Xie, Ian P. Castro,
        "Effcient generation of inflow conditions for large eddy simulation
        of street-scale flows",
        Flow Turbulence Combust, (2008) 81:449-470

    Reynolds stress, mean velocity and turbulence length scale values can 
    either be sepcified directly or computed using a given distribution function        

See also
    Foam::fixedValueFvPatchVectorField

SourceFiles
    turbulentMeanInletFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentMeanInletFvPatchVectorField_H
#define turbulentMeanInletFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "Random.H"
#include "labelVector.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class pointToPointPlanarInterpolation;

/*---------------------------------------------------------------------------*\
             Class turbulentMeanInletFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class turbulentMeanInletFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{

    // Private data

        // Patch information

            //- Flag to identify if initialized
            bool isInitialized_;

            //- local x-axis
            vector patchNormal_;

        // Table reading for patch inlet flow properties

            //- Fraction of perturbation (fraction of bounding box) to add
            scalar perturb_;

            //- Interpolation scheme to use
            word mapMethod_;

            //- 2D interpolation (for 'planarInterpolation' mapMethod)
            mutable autoPtr<pointToPointPlanarInterpolation> mapperPtr_;

            //- Flag to identify to calculate the U field
            bool interpolateU_;

            //- Inlet mean velocity
            scalarField U_;


    // Private Member Functions

        //- Initialise
        void initialise();

        //- Helper function to calculate values from the boundary data
        //- Return a reference to the patch mapper object
        const pointToPointPlanarInterpolation& patchMapper() const;

        //- Helper function to calculate values from the boundary data or
        //  read from dictionary
        template<class Type>
        tmp<Field<Type>> interpolateOrRead
        (
            const word& fieldName,
            const dictionary& dict,
            bool& calculateField
        ) const;

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
    TypeName("turbulentMeanInlet");


    // Constructors

        //- Construct from patch and internal field
        turbulentMeanInletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentMeanInletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given turbulentMeanInletFvPatchVectorField
        //  onto a new patch
        turbulentMeanInletFvPatchVectorField
        (
            const turbulentMeanInletFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        turbulentMeanInletFvPatchVectorField
        (
            const turbulentMeanInletFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new turbulentMeanInletFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        turbulentMeanInletFvPatchVectorField
        (
            const turbulentMeanInletFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new turbulentMeanInletFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<vector>&,
                const labelList&
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "turbulentMeanInletFvPatchVectorFieldTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
