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

#include "Time.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::turbulentKineticEnergyInletFvPatchScalarField::calculateOrRead
(
    const word& fieldName,
    const dictionary& dict,
    bool& calculateField
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

        calculateField = false;
        return tFld;
    }
    else
    {
        calculateField = true;
        return calculateBoundaryData<Type>(fieldName, dict);
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::turbulentKineticEnergyInletFvPatchScalarField::calculateBoundaryData
(
    const word& fieldName,
    const dictionary& dict
) const
{
    word subDictName(fieldName+"_Profile");

    Field<Type> value(this->patch().size(), pTraits<Type>::zero);

    if (dict.found(subDictName))
    {
        const dictionary& subDict = dict.subDict(subDictName);

        word profile(subDict.lookupOrDefault<word>("profile", "uniform"));
        Type referenceValue(subDict.lookupOrDefault<Type>("referenceValue", pTraits<Type>::zero));

        if (profile == "uniform")
        {
            value = referenceValue;
        }
        else if (profile == "exponential")
        {
            vector dir(subDict.lookupOrDefault<vector>("direction", vector(0,1,0)));
            vector origin(subDict.lookupOrDefault<vector>("origin", vector::zero));
            vector refPoint(subDict.lookupOrDefault<vector>("referencePoint", vector(0,1,0)));

            Type alpha(subDict.lookupOrDefault<Type>("alpha", pTraits<Type>::zero));
            Type constant(subDict.lookupOrDefault<Type>("constant", pTraits<Type>::zero));

            const scalarField dirCmpt((dir&(this->patch().Cf()-origin))/(dir&(refPoint-origin)));

            forAll(value, label)
            {
                value[label] = cmptMultiply
                (
                    referenceValue,
                    cmptPow(pTraits<Type>::one*fabs(dirCmpt[label]), alpha)
                );
           }

            value = value + constant;
        }
        else if (profile == "parabolic")
        {
            vector dir(subDict.lookupOrDefault<vector>("direction", vector(0,1,0)));
            vector origin(subDict.lookupOrDefault<vector>("origin", vector::zero));
            vector refPoint(subDict.lookupOrDefault<vector>("referencePoint", vector(0,1,0)));

            Type constant(subDict.lookupOrDefault<Type>("constant", pTraits<Type>::zero));

            const scalarField dirCmpt((dir&(this->patch().Cf()-origin))/(dir&(refPoint-origin)));

            forAll(value, label)
            {
                value[label] = referenceValue*(1-sqr(dirCmpt[label]-1));
            }
            value = value + constant;
        }
        else
        {
            Info << "profile " << profile << " does not exist (ERROR)" << endl;
        }
    }
    else
    {
        Info << "parameters for " << fieldName << " does not exist (ERROR)" << endl;
    }

    tmp<Field<Type>> tFld(new Field<Type>(value));

    return tFld;
}


// ************************************************************************* //
