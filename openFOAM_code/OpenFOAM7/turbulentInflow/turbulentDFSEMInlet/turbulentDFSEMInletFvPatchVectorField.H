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
    Foam::turbulentDFSEMInletFvPatchVectorField

SeeAlso
    fixedValueFvPatchField

SourceFiles
    turbulentDFSEMInletFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentDFSEMInletFvPatchVectorField_H
#define turbulentDFSEMInletFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "Random.H"
#include "dfeddy.H"
#include "pointIndexHit.H"
#include "instantList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class pointToPointPlanarInterpolation;

/*---------------------------------------------------------------------------*\
           Class turbulentDFSEMInletFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class turbulentDFSEMInletFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private Data

        //- Maximum number of attempts when seeding eddies
        static label seedIterMax_;

        //- Characteristic length scale, e.g. half channel height
        const scalar delta_;

        //- Ratio of sum of eddy volumes to eddy box volume; default = 1
        const scalar density_;

        //- Von Karman constant
        const scalar kappa_;


        // Table reading for patch inlet flow properties

            //- Fraction of perturbation (fraction of bounding box) to add
            scalar perturb_;

            //- Interpolation scheme to use (nearestCell | planarInterpolation)
            word mapMethod_;

            //- 2D interpolation (for 'planarInterpolation' mapMethod)
            mutable autoPtr<pointToPointPlanarInterpolation> mapperPtr_;

            //- Flag to identify to interpolate the R field
            bool interpolateR_;

            //- Flag to identify to interpolate the L field
            bool interpolateL_;

            //- Flag to identify to interpolate the U field
            bool interpolateU_;

            //- Reynolds stress tensor
            symmTensorField R_;

            //- Length scale
            scalarField L_;

            //- Inlet velocity
            scalarField U_;

            //- Mean inlet velocity
            scalar UMean_;

            //- Max normal velocity
            scalar UMax_;


        // Patch information

            //- Patch area - total across all processors
            scalar patchArea_;

            //- Decomposed patch faces as a list of triangles
            faceList triFace_;

            //- Addressing from per triangle to patch face
            labelList triToFace_;

            //- Cumulative triangle area per triangle face
            scalarList triCumulativeMagSf_;

            //- Cumulative area fractions per processor
            scalarList sumTriMagSf_;

            bool periodicInY_;

            bool periodicInZ_;

            scalar patchSpanY_;

            scalar patchSpanZ_;


        //- List of eddies
        List<dfeddy> eddies_;

        //- Global number of eddies
        label nEddyGlobal_;

        //- Local number of eddies
        label nEddyLocal_;

        List<label> eddyLabel_;

        List<vector> eddyPosition_;

        List<scalar> eddyDistance_;

        List<vector> eddyScale_;

        List<vector> eddyIntensity_;

        //- Minimum number of cells required to resolve an eddy
        label nCellPerEddy_;

        //- Patch normal into the domain
        vector patchNormal_;

        //- Eddy box volume
        scalar v0_;

        //- Random number generator
        Random rndGen_;

        //- Length scale per patch face
        scalarField sigmax_;

        //- Maximum length scale (across all processors)
        scalar maxSigmaX_;

        //- Current time index (used for updating)
        label curTimeIndex_;

        //- Patch bounds (local processor)
        boundBox patchBounds_;

        //- Single processor contains all eddies (flag)
        bool singleProc_;

        //- Flag to identify if clean restart
        bool isCleanRestart_;

        //- Flag to identify if restart
        bool isRestart_;


    // Private Member Functions

        //- Initialise info for patch point search
        void initialisePatch();
        
        //- Initialise and check turbulence parameters
        void initialiseParameters();

        //- Initialise the eddy box
        void initialiseEddyBox();

        //- Set a new eddy position
        pointIndexHit setNewPosition(const bool global);

        //- Initialise eddies
        void initialiseEddies();

        //- Convect the eddies
        void convectEddies(const scalar deltaT);

        //- Calculate the velocity fluctuation at a point
        vectorField uDashEddy(const List<dfeddy>&, const pointField&) const;

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

        //- Return a reference to the patch mapper object
        const pointToPointPlanarInterpolation& patchMapper() const;

        //- Return eddies from remote processors that interact with local
        //- processor
        void calcOverlappingProcEddies
        (
            List<List<dfeddy>>& overlappingEddies
        ) const;


public:

   //- Runtime type information
   TypeName("turbulentDFSEMInlet");


   // Constructors

        //- Construct from patch and internal field
        turbulentDFSEMInletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentDFSEMInletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given turbulentDFSEMInletFvPatchVectorField
        //- onto a new patch
        turbulentDFSEMInletFvPatchVectorField
        (
            const turbulentDFSEMInletFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        turbulentDFSEMInletFvPatchVectorField
        (
            const turbulentDFSEMInletFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new turbulentDFSEMInletFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        turbulentDFSEMInletFvPatchVectorField
        (
            const turbulentDFSEMInletFvPatchVectorField&,
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
                new turbulentDFSEMInletFvPatchVectorField(*this, iF)
            );
        }


    //- Destructor
    virtual ~turbulentDFSEMInletFvPatchVectorField() = default;


    // Member Functions

        //- Helper function to check that Reynold stresses are valid
        static bool checkStresses(const symmTensorField& Rf);


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const fvPatchFieldMapper& m);

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchVectorField& ptf,
                const labelList& addr
            );


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
    #include "turbulentDFSEMInletFvPatchVectorFieldTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
