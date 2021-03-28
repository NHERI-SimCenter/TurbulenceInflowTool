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
    Foam::turbulentDFMInletFvPatchVectorField

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
    turbulentDFMInletFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentDFMInletFvPatchVectorField_H
#define turbulentDFMInletFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "Random.H"
#include "labelVector.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class pointToPointPlanarInterpolation;

/*---------------------------------------------------------------------------*\
             Class turbulentDFMInletFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class turbulentDFMInletFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{

    // Private data

        // Patch information

            //- Current time index (used for updating)
            label curTimeIndex_;

            //- patch size
            labelList patchSize_;

            //- local x-axis
            vector patchNormal_;

            bool periodicInY_;

            bool periodicInZ_;

            //- Intermediate spatially correlated fluctuation
            vectorField uFluctFiltered_;

            //- Old temporal correlated fluctuation
            vectorField uFluctTemporalOld_;

            //- New temporal correlated fluctuation
            vectorField uFluctTemporal_;


        // Table reading for patch inlet flow properties

            //- Fraction of perturbation (fraction of bounding box) to add
            scalar perturb_;

            //- Interpolation scheme to use
            word mapMethod_;

            //- 2D interpolation (for 'planarInterpolation' mapMethod)
            mutable autoPtr<pointToPointPlanarInterpolation> mapperPtr_;

            //- Flag to identify to calculate the U field
            bool interpolateU_;

            //- Flag to identify to calculate the I field
            bool interpolateR_;

            //- Flag to identify to calculate the L field
            bool interpolateL_;

            //- Inlet mean velocity
            scalarField U_;

            //- turbulent intensity tensor
            symmTensorField R_;

            //- Length scales
            tensorField L_;
            
            //- Converted length scales
            tensorField L0_;

            //- Lund coefficients
            tensorField Lund_;


        // Virtual grid points infomation

            //- Flag to identify if initialized
            bool isInitialized_;

            //- Flag to identify if clean restart
            bool isCleanRestart_;

            //- Flag to identify if restart
            bool isRestart_;

            //- Ratio of virtual gird spacing to mesh size
            scalar gridFactor_;

            //- Origian of virtual grid points
            vector origin_;

            //- Number of girds along y axis
            label My_;

            //- Number of girds along z axis
            label Mz_;

            //- Grid spacing
            scalar delta_;

            //- y direction length scale to grid spacing ratio
            Field<labelVector> ny_;

            //- z direction length scale to grid spacing ratio
            Field<labelVector> nz_;

            //- Filter width to length scale ratio
            label nfK_;

            //- y indices
            labelList yindices_;

            //- z indices
            labelList zindices_;


    // Parallel processing information

            //- Distributing virtual grids per Proc
            label indicesPerProc_;

            //- Procs will do +1 Indices
            label rest_;

    // Random generation

            //- Random number generator
            Random rndGen_;

            //- filter function shape
            const word filterType_;

            //- Random field
            labelVector rndSize_;

            //- Filter coefficients for u component
            scalarListList filterCoeffProcx;

            //- Filter coefficients for v component
            scalarListList filterCoeffProcy;

            //- Filter coefficients for v component
            scalarListList filterCoeffProcz;


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

        //- Initialise
        void initialise();

        //- Initialise info for output velocities
        void initialiseOutput();

        //- Initialise virtual grid points
        void initialisePatch();

        //- Initialise and check turbulence parameters
        void initialiseParameters();

        //- Initialise virtual grid points
        void initialiseVirtualGrid();

        //- Initialise digital filters
        void initialiseFilterCoeff();

        //- Helper function to calculate values from the boundary data or
        //  read from dictionary
        template<class Type>
        tmp<Field<Type>> interpolateOrRead
        (
            const word& fieldName,
            const dictionary& dict,
            bool& calculateField
        ) const;

        //- Helper function to calculate values from the boundary data
        //- Return a reference to the patch mapper object
        const pointToPointPlanarInterpolation& patchMapper() const;

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

        //- Helper function to gather values from all processors
        template<class Type>
        Field<Type> gatherProc
        (
            const Field<Type>& valsProc
        ) const;

        //- Returns index of a 2D array in list format
        inline label get1DIndex(label x, label y, label yDim);

        //- Get one dimensional filter coefficents
        void get1DFilterCoeff(scalarList& b, label n);

        //- Get one dimensional filter coefficents
        void get2DFilterCoeff(scalarList& b, label ny, label nz);

        //- Generate normal distributed random field
        scalarField getRandomField(label);

        //- Create spatial correlated random field
        void spatialCorr();

        //- Create temporal correlation
        void temporalCorr();

        //- Modified Bessel function of the second kind
        scalar bessi0(const scalar);

        //- Modified Bessel function of the second kind
        scalar bessk0(const scalar);

public:

    //- Runtime type information
    TypeName("turbulentDFMInlet");


    // Constructors

        //- Construct from patch and internal field
        turbulentDFMInletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentDFMInletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given turbulentDFMInletFvPatchVectorField
        //  onto a new patch
        turbulentDFMInletFvPatchVectorField
        (
            const turbulentDFMInletFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        turbulentDFMInletFvPatchVectorField
        (
            const turbulentDFMInletFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new turbulentDFMInletFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        turbulentDFMInletFvPatchVectorField
        (
            const turbulentDFMInletFvPatchVectorField&,
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
                new turbulentDFMInletFvPatchVectorField(*this, iF)
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
    #include "turbulentDFMInletFvPatchVectorFieldTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
