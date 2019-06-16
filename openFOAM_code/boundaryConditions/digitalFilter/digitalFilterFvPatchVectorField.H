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
    Foam::digitalFilterFvPatchVectorField

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

Usage
    Property    | Description                             | Required | Default value
    U           | Mean velocity field                     | yes      |
    I           | Turbulent intensity tensor field        | yes      |
    Lu          | Turbulence length scale for u component | yes      |
    Lv          | Turbulence length scale for v component | yes      |
    Lw          | Turbulence length scale for w component | yes      |
    shape       | Digital filter function                 | no       | exponential
    girdFactor  | Gird spacing to mesh size ratio         | no       | 1
    nfK         | Filter length to length scale ratio     | no       | 4          

See also
    Foam::fixedValueFvPatchVectorField

SourceFiles
    digitalFilterFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef digitalFilterFvPatchVectorField_H
#define digitalFilterFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "Random.H"
#include "labelVector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class pointToPointPlanarInterpolation;

/*---------------------------------------------------------------------------*\
             Class digitalFilterFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class digitalFilterFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{

    // Private data

        // Patch information

            //- Current time index (used for updating)
            label curTimeIndex_;

            //- patch size
            labelList patchSize_;

            //- Intermediate spatially correlated fluctuation
            vectorField uFluctFiltered_;

            //- Old temporal correlated fluctuation
            vectorField uFluctTemporalOld_;

            //- New temporal correlated fluctuation
            vectorField uFluctTemporal_;


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

            //- Length scale for u component
            vectorField Lu_;

            //- Length scale for v component
            vectorField Lv_;

            //- Length scale for w component
            vectorField Lw_;

            //- Lund coefficients
            tensorField Lund_;


        // Virtual grid points infomation

            //- Flag to identify if initialized
            bool isInitialized_;

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

            //- Maximum y direction length scale to grid spacing ratio
            labelVector nyMax_;

            //- Maximum z direction length scale to grid spacing ratio
            labelVector nzMax_;

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

            //- Patch to virtual grid mapping (global to gloabl)
            labelList facesToIndices_;

            //- Virtual grid to patch mapping (global to local)
            labelList indicesToFaces_;

    // Random generation

            //- Random number generator
            Random rndGen_;

            //- filter function shape
            const word filterShape_;

            //- Virtual grid points
            pointField virualGridPoints_;

            //- Random field
            labelVector rndSize_;

            //- Filter coefficients for u component
            scalarListList filterCoeffProc_u;

            //- Filter coefficients for v component
            scalarListList filterCoeffProc_v;

            //- Filter coefficients for v component
            scalarListList filterCoeffProc_w;


    // Private Member Functions

        //- Initialise virtual grid points
        void initialisePatch();

        //- Initialise virtual grid points
        void initialiseVirtualGrid();

        //- Initialise digital filters
        void initialiseDigitalFilter();

        void sortfaceCentres(const vectorField& fCentres_, labelList& order_);

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

public:

    //- Runtime type information
    TypeName("digitalFilter");


    // Constructors

        //- Construct from patch and internal field
        digitalFilterFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        digitalFilterFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given digitalFilterFvPatchVectorField
        //  onto a new patch
        digitalFilterFvPatchVectorField
        (
            const digitalFilterFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        digitalFilterFvPatchVectorField
        (
            const digitalFilterFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new digitalFilterFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        digitalFilterFvPatchVectorField
        (
            const digitalFilterFvPatchVectorField&,
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
                new digitalFilterFvPatchVectorField(*this, iF)
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
    #include "digitalFilterFvPatchVectorFieldTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
