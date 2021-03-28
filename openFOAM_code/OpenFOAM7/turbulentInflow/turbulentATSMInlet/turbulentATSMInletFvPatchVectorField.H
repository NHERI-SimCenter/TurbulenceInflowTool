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
    Foam::turbulentATSMInletFvPatchVectorField

SeeAlso
    fixedValueFvPatchField

SourceFiles
    turbulentATSMInletFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentATSMInletFvPatchVectorField_H
#define turbulentATSMInletFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "Random.H"
#include "vorton.H"
#include "pointIndexHit.H"
#include "instantList.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class pointToPointPlanarInterpolation;

/*---------------------------------------------------------------------------*\
           Class turbulentATSMInletFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class turbulentATSMInletFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private Data

        //- Maximum number of attempts when seeding vortons
        static label seedIterMax_;

        //- Ratio of sum of eddy volumes to eddy box volume; default = 1
        const scalar density_;

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
            vectorField L_;

            //- Inlet velocity
            scalarField U_;

            //- Mean inlet velocity
            scalar UMean_;

            //- Max normal velocity
            scalar UMax_;


        // Patch information

            //- Patch area - total across all processors
            scalar patchArea_;

            //- patch size
            labelList patchSize_;

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


        //- List of vortons
        List<vorton> vortons_;

        word vortonType_;

        //- Global number of vortons
        label nVortonGlobal_;

        //- Local number of vortons
        label nVortonLocal_;

        List<label> vortonLabel_;

        List<vector> vortonPosition_;

        List<scalar> vortonDistance_;

        List<vector> vortonScale_;

        List<vector> vortonIntensity_;

        //- Patch normal into the domain
        vector patchNormal_;

        //- Vorton box volume
        scalar v0_;

        //- Random number generator
        Random rndGen_;

        //- Maximum length scale (across all processors)
        scalar maxSigmaX_;

        //- Current time index (used for updating)
        label curTimeIndex_;

        //- Patch bounds (local processor)
        boundBox patchBounds_;

        //- Single processor contains all vortons (flag)
        bool singleProc_;

        //- Flag to identify if clean restart
        bool isCleanRestart_;

        //- Flag to identify if restart
        bool isRestart_;

        // Output information

            //- number of output faces
            label nOutputFace_;

            //- output face indices
            labelList outputFaceIndices_;

            //- File pointer
            PtrList<OFstream> filePtrs_;


    // Private Member Functions

        //- create files
        void createFiles();

        void writeFileHeader(const label);

        void writeValues(const label, const vector);

        //- Return width of character stream output
        label charWidth() const;

        //- Initialise the output stream for writing
        void initStream(Ostream& os) const;

        //- Write a commented string to stream
        void writeCommented(Ostream& os, const string& str) const;

        //- Write a tabbed string to stream
        void writeTabbed(Ostream& os, const string& str) const;

        //- Write a commented header to stream
        void writeHeader(Ostream& os, const string& str) const;

        //- Write a (commented) header property and value pair
        template<class Type>
        void writeHeaderValue
        (
            Ostream& os,
            const string& property,
            const Type& value
        ) const;

        //- Write the current time to stream
        void writeTime(Ostream& os) const;

        //- Initialise info for output velocities
        void initialiseOutput();

        //- Initialise info for patch point search
        void initialisePatch();
        
        //- Initialise and check turbulence parameters
        void initialiseParameters();

        //- Initialise the eddy box
        void initialiseVortonBox();

        //- Set a new eddy position
        pointIndexHit setNewPosition(const bool global);

        //- Initialise vortons
        void initialiseVortons();

        //- Convect the vortons
        void convectVortons(const scalar deltaT);

        //- Calculate the velocity fluctuation at a point
        vectorField uDashVorton(const List<vorton>&, const pointField&) const;

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

        //- Return vortons from remote processors that interact with local
        //- processor
        void calcOverlappingProcVortons
        (
            List<List<vorton>>& overlappingVortons
        ) const;


public:

   //- Runtime type information
   TypeName("turbulentATSMInlet");


   // Constructors

        //- Construct from patch and internal field
        turbulentATSMInletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentATSMInletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given turbulentATSMInletFvPatchVectorField
        //- onto a new patch
        turbulentATSMInletFvPatchVectorField
        (
            const turbulentATSMInletFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        turbulentATSMInletFvPatchVectorField
        (
            const turbulentATSMInletFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new turbulentATSMInletFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        turbulentATSMInletFvPatchVectorField
        (
            const turbulentATSMInletFvPatchVectorField&,
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
                new turbulentATSMInletFvPatchVectorField(*this, iF)
            );
        }


    //- Destructor
    virtual ~turbulentATSMInletFvPatchVectorField() = default;


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
    #include "turbulentATSMInletFvPatchVectorFieldTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
