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
    Foam::turbulentDFInletFvPatchVectorField

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
    either be sepcified directly, or mapped. If mapping, the values should
    be entered in the same form as the timeVaryingMappedFixedValue condition,
    except that no interpolation in time is supported. These should be located
    in the directory:

        \$FOAM_CASE/constant/boundaryData/<patchName>/points
        \$FOAM_CASE/constant/boundaryData/<patchName>/0/U|R|Lu|Lv|Lw

Usage
    Property    | Description                             | Required | Default value
    U           | Mean velocity field                     | yes      |
    R           | Renoyld Stress tensor field             | yes      |
    Lu          | Turbulence length scale for u component | yes      |
    Lv          | Turbulence length scale for v component | yes      |
    Lw          | Turbulence length scale for w component | yes      |
    shape       | Digital filter function                 | no       | exponential
    girdFactor  | Gird spacing to mesh size ratio         | no       | 1
    nfK         | Filter length to length scale ratio     | no       | 4          
    mapMethod   | Method to map reference values          | no       | nearestCell
    perturb     | Point perturbation for interpolation    | no       | 1e-5

See also
    Foam::fixedValueFvPatchVectorField

SourceFiles
    turbulentDFInletFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentDFInletFvPatchVectorField_H
#define turbulentDFInletFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "Random.H"
#include "labelVector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class pointToPointPlanarInterpolation;

/*---------------------------------------------------------------------------*\
             Class turbulentDFInletFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class turbulentDFInletFvPatchVectorField
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

            //- Fraction of perturbation (fraction of bounding box) to add
            scalar perturb_;

            //- Interpolation scheme to use
            word mapMethod_;

            //- 2D interpolation (for 'planarInterpolation' mapMethod)
            mutable autoPtr<pointToPointPlanarInterpolation> mapperPtr_;

            //- Flag to identify to interpolate the R field
            bool interpolateR_;

            //- Reynolds stress tensor
            symmTensorField R_;

            //- Flag to identify to interpolate the L field
            bool interpolateL_;

            //- Length scale for u component
            vectorField Lu_;

            //- Length scale for v component
            vectorField Lv_;

            //- Length scale for w component
            vectorField Lw_;

            //- Flag to identify to interpolate the U field
            bool interpolateU_;

            //- Inlet mean velocity
            vectorField U_;

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

        //- Return a reference to the patch mapper object
        const pointToPointPlanarInterpolation& patchMapper() const;

        //- Helper function to interpolate values from the boundary data or
        //  read from dictionary
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
    TypeName("turbulentDFInlet");


    // Constructors

        //- Construct from patch and internal field
        turbulentDFInletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentDFInletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given turbulentDFInletFvPatchVectorField
        //  onto a new patch
        turbulentDFInletFvPatchVectorField
        (
            const turbulentDFInletFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        turbulentDFInletFvPatchVectorField
        (
            const turbulentDFInletFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new turbulentDFInletFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        turbulentDFInletFvPatchVectorField
        (
            const turbulentDFInletFvPatchVectorField&,
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
                new turbulentDFInletFvPatchVectorField(*this, iF)
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
    #include "turbulentDFInletFvPatchVectorFieldTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
