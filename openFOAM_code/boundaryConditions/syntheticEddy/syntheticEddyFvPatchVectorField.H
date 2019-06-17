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
    Foam::syntheticEddyFvPatchVectorField

Description
    Turbulent velocity boundary condition based on the Synthetic-Eddy Method.

    Reference: 
        N. Jarrin, S. Benhamadouche, D. Laurence, R. Prosser,
        "A synthetic-eddy-method for generating inflow conditions for large
        eddy simulations",
        International Journal of Heat and Fluid Flow (2006) 27:585-593

        N. Jarrin, R. Prosser, J.-C. Uribe, S. Benhamadouche, D. Laurence,
        "Reconstruction of turbulent fluctuations for hubrid RANS/LES simulations
        using a Synthetic-Eddy Method",
        International Journal of Heat and Fluid Flow (2009) 30:435-442

    Reynolds stress, mean velocity and turbulence length scale values can 
    either be sepcified directly or computed using a given distribution function


Usage
    Property    | Description                             | Required | Default value
    U           | Mean velocity field                     | yes      |
    I           | Turbulent intensity tensor field        | yes      |
    Lu          | Turbulence length scale for u component | yes      |
    Lv          | Turbulence length scale for v component | yes      |
    Lw          | Turbulence length scale for w component | yes      |
    shape       | Velocity distribution function          | no       | tent
    eddyDensity | Eddy density                            | no       | 1

SeeAlso
    timeVaryingMappedFixedValueFvPatchField

SourceFiles
    syntheticEddyFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef syntheticEddyFvPatchVectorField_H
#define syntheticEddyFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "Random.H"
#include "eddy.H"
#include "pointIndexHit.H"
#include "instantList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class pointToPointPlanarInterpolation;

/*---------------------------------------------------------------------------*\
           Class syntheticEddyFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class syntheticEddyFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Current time index (used for updating)
        label curTimeIndex_;

        // Table reading for patch inlet flow properties

            //- Dynamic mesh dictionary
            dictionary inflowProperties_;

            //- Flag to identify to calculate the U field
            bool calculateU_;

            //- Flag to identify to calculate the I field
            bool calculateI_;

            //- Flag to identify to calculate the L field
            bool calculateL_;

            //- Inlet mean velocity
            vectorField U_;

            //- turbulent intensity tensor
            symmTensorField I_;

            //- Reynold stress tensor
            symmTensorField R_;

            //- Length scale for u component
            vectorField Lu_;

            //- Length scale for v component
            vectorField Lv_;

            //- Length scale for w component
            vectorField Lw_;


        // Patch information

            //- Patch bounds (local processor)
            boundBox patchBounds_;

            //- Patch area - total across all processors
            scalar patchArea_;

            //- Patch normal into the domain
            vector patchNormal_;

            //- Decomposed patch faces as a list of triangles
            faceList triFace_;

            //- Addressing from per triangle to patch face
            labelList triToFace_;

            //- Cumulative triangle area per triangle face
            scalarList triCumulativeMagSf_;

            //- Cumulative area fractions per processor
            scalarList sumTriMagSf_;

            //- Eddy box volume
            scalar V0_;


        // Eddy information

            //- Ratio of sum of eddy volumes to eddy box volume; default = 1
            const scalar eddyDensity_;

            //- Global number of eddies
            label nEddy_;

            //- List of eddies
            List<eddy> eddies_;

            //- velocity distribution shape
            const word velocityShape_;

            //- Random number generator
            Random rndGen_;

            //- Maximum length scale (across all processors)
            scalar maxSigmaX_;

            //- Single processor contains all eddies (flag)
            bool singleProc_;

            //- Maximum number of attempts when seeding eddies
            static label seedIterMax_;


    // Private Member Functions

        //- Initialise info for patch point search
        void initialisePatch();

        //- Initialise the eddy box
        void initialiseEddyBox();

        //- Set a new eddy position
        pointIndexHit setNewPosition(const bool global);

        //- Initialise eddies
        void initialiseEddies();

        //- Convect the eddies
        void convectEddies(const scalar deltaT);

        //- Calculate the velocity fluctuation at a point
        vector uDashEddy(const List<eddy>& eddies, const point& globalX) const;

        //- Return sigma to length scale ratio
        scalar sigmaToLengthscaleRatio() const;

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

        //- Return a reference to the patch mapper object
        const pointToPointPlanarInterpolation& patchMapper() const;

        //- Return eddies from remote processors that interact with local
        //  processor
        void calcOverlappingProcEddies
        (
            List<List<eddy>>& overlappingEddies
        ) const;


public:

   //- Runtime type information
   TypeName("syntheticEddy");


   // Constructors

        //- Construct from patch and internal field
        syntheticEddyFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        syntheticEddyFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given syntheticEddyFvPatchVectorField
        //  onto a new patch
        syntheticEddyFvPatchVectorField
        (
            const syntheticEddyFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        syntheticEddyFvPatchVectorField
        (
            const syntheticEddyFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new syntheticEddyFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        syntheticEddyFvPatchVectorField
        (
            const syntheticEddyFvPatchVectorField&,
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
                new syntheticEddyFvPatchVectorField(*this, iF)
            );
        }

    virtual ~syntheticEddyFvPatchVectorField();


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "syntheticEddyFvPatchVectorFieldTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
