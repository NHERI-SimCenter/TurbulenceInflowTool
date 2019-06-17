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

#include "pointToPointPlanarInterpolation.H"
#include "Time.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::turbulentSEMInletFvPatchVectorField::interpolateOrRead
(
    const word& fieldName,
    const dictionary& dict,
    bool& interpolateField
) const
{
    if (dict.found(fieldName))
    {
        tmp<Field<Type>> tFld
        (
            new Field<Type>
            (
                fieldName,
                dict,
                this->patch().size()
            )
        );

        interpolateField = false;
        return tFld;
    }
    else
    {
        interpolateField = true;
        return interpolateBoundaryData<Type>(fieldName);
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::turbulentSEMInletFvPatchVectorField::interpolateBoundaryData
(
    const word& fieldName
) const
{
    const word& patchName = this->patch().name();

    fileName valsFile
    (
        fileHandler().filePath
        (
            fileName
            (
                this->db().time().path()
               /this->db().time().caseConstant()
               /"boundaryData"
               /patchName
               /"0"
               /fieldName
            )
        )
    );

    autoPtr<ISstream> isPtr
    (
        fileHandler().NewIFstream
        (
            valsFile
        )
    );

    Field<Type> vals(isPtr());

    Info<< "Turbulent SEM patch " << patchName
        << ": interpolating field " << fieldName
        << " from " << valsFile << endl;

    return patchMapper().interpolate(vals);
}


// ************************************************************************* //
