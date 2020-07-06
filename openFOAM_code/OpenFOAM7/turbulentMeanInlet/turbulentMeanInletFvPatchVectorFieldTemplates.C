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
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::turbulentMeanInletFvPatchVectorField::interpolateOrRead
(
    const word& fieldName,
    const dictionary& dict,
    bool& interpolateField
) const
{
    const word calculateName("calculate"+fieldName);
    const bool calculateFlag(dict.lookupOrDefault<bool>(calculateName, false));

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
        if (calculateFlag)
        {

            IOdictionary inflowProperties
            (
                IOobject
                (
                    "inflowProperties",
                    this->db().time().constant(),
                    this->internalField().mesh(),
                    IOobject::MUST_READ_IF_MODIFIED,
                    IOobject::NO_WRITE,
                    false
                )
            );

            return calculateBoundaryData<Type>(fieldName, inflowProperties);
        }
        else
        {
            interpolateField = true;
            return interpolateBoundaryData<Type>(fieldName);
        }
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::turbulentMeanInletFvPatchVectorField::interpolateBoundaryData
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

    Info<< "Turbulent Mean patch " << patchName
        << ": interpolating field " << fieldName
        << " from " << valsFile << endl;

    return patchMapper().interpolate(vals);
}

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::turbulentMeanInletFvPatchVectorField::calculateBoundaryData
(
    const word& fieldName,
    const dictionary& dict
) const
{
    word subDictName(fieldName+"Dict");
    Field<Type> value(this->patch().size(), pTraits<Type>::zero);

    const vectorField nf(patch().nf());
    const vector e1 = -gAverage(nf);

    vector N(dict.lookupOrDefault<vector>("Naxis", vector(0,0,0)));
    vector e3(vector::zero);

    if (mag(N) > 0)
    {
        N /= mag(N);
        e3 = e1^N;

        if (mag(e3) == 0)
        {
            scalar beta(dict.lookupOrDefault<scalar>("beta", 0));
            e3 = vector
            (
                ::sin(beta*constant::mathematical::twoPi/360)*N[1],
               -::sin(beta*constant::mathematical::twoPi/360)*N[0],
                ::cos(beta*constant::mathematical::twoPi/360)
            );
        }
    }
    else
    {
        e3 = vector(0,0,1);
    }

    coordinateSystem patchCoord
    (
        "patchCoord",
        vector::zero,
        e3,
        e1
    );

    const polyPatch& polyPatch = this->patch().patch();
    const pointField localPoints = patchCoord.localPosition(polyPatch.points());

    boundBox patchBounds(localPoints);

    const vector offset(dict.lookupOrDefault<vector>("offset", vector::zero));
    const vector origin = patchBounds.min()+offset;

    patchCoord.origin() = patchCoord.globalPosition(origin);

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
            vectorField fCentres(patchCoord.localPosition(this->patch().Cf()));

            const scalar refAngl(subDict.lookupOrDefault<scalar>("referenceAngl", 0.0));
            const scalar refDist(subDict.lookupOrDefault<scalar>("referenceDist", 1.0));

            vector refDire(vector::zero);
            refDire.component(vector::Y) = ::sin(refAngl*constant::mathematical::twoPi/360);
            refDire.component(vector::Z) = ::cos(refAngl*constant::mathematical::twoPi/360);

            if (refDist <= 0)
            {
                Info<<"reference distance of the " << fieldName << " field is no larger than zero (ERROR)" << endl;
            }

            const scalarField dirCmpt((refDire&fCentres)/refDist);

            Type alpha(subDict.lookupOrDefault<Type>("alpha", pTraits<Type>::zero));

            forAll(value, label)
            {
                value[label] = cmptMultiply
                (
                    referenceValue,
                    cmptPow(pTraits<Type>::one*fabs(dirCmpt[label]), alpha)
                );
            }
        }
        else if (profile == "linear")
        {
            vectorField fCentres(patchCoord.localPosition(this->patch().Cf()));

            const scalar refAngl(subDict.lookupOrDefault<scalar>("referenceAngl", 0.0));
            const scalar refDist(subDict.lookupOrDefault<scalar>("referenceDist", 1.0));

            vector refDire(vector::zero);
            refDire.component(vector::Y) = ::sin(refAngl*constant::mathematical::twoPi/360);
            refDire.component(vector::Z) = ::cos(refAngl*constant::mathematical::twoPi/360);

            if (refDist <= 0)
            {
                Info<<"reference distance of the " << fieldName << " field is no larger than zero (ERROR)" << endl;
            }

            const scalarField dirCmpt((refDire&fCentres)/refDist);

            Type alpha(subDict.lookupOrDefault<Type>("alpha", pTraits<Type>::zero));

            forAll(value, label)
            {
                value[label] = cmptMultiply
                (
                    alpha*(dirCmpt[label]-1.0)+pTraits<Type>::one,
                    referenceValue
                );
            }
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

    Info<< "Turbulent Mean patch " << this->patch().name()
        << ": calculating field " << fieldName
        << " from " << dict.name() << endl;

    return tFld;
}

template<>
Foam::tmp<Foam::symmTensorField>
Foam::turbulentMeanInletFvPatchVectorField::calculateBoundaryData<Foam::symmTensor>
(
    const word& fieldName,
    const dictionary& dict
) const
{
    word subDictName(fieldName+"Dict");
    symmTensorField value(this->patch().size(), symmTensor::zero);

    const vectorField nf(patch().nf());
    const vector e1 = -gAverage(nf);

    vector N(dict.lookupOrDefault<vector>("Naxis", vector(0,0,0)));
    vector e3(vector::zero);

    if (mag(N) > 0)
    {
        N /= mag(N);
        e3 = e1^N;

        if (mag(e3) == 0)
        {
            scalar beta(dict.lookupOrDefault<scalar>("beta", 0));
            e3 = vector
            (
                ::sin(beta*constant::mathematical::twoPi/360)*N[1],
               -::sin(beta*constant::mathematical::twoPi/360)*N[0],
                ::cos(beta*constant::mathematical::twoPi/360)
            );
        }
    }
    else
    {
        e3 = vector(0,0,1);
    }

    coordinateSystem patchCoord
    (
        "patchCoord",
        vector::zero,
        e3,
        e1
    );

    const polyPatch& polyPatch = this->patch().patch();
    const pointField localPoints = patchCoord.localPosition(polyPatch.points());

    boundBox patchBounds(localPoints);

    const vector offset(dict.lookupOrDefault<vector>("offset", vector::zero));
    const vector origin = patchBounds.min()+offset;

    patchCoord.origin() = patchCoord.globalPosition(origin);

    if (dict.found(subDictName))
    {
        const dictionary& subDict = dict.subDict(subDictName);

        word profile(subDict.lookupOrDefault<word>("profile", "uniform"));
        symmTensor referenceValue(subDict.lookupOrDefault<symmTensor>("referenceValue", symmTensor::zero));

        if (profile == "uniform")
        {
            value = referenceValue;
        }
        else if (profile == "exponential")
        {
            vectorField fCentres(patchCoord.localPosition(this->patch().Cf()));

            const scalar refAngl(subDict.lookupOrDefault<scalar>("referenceAngl", 0.0));
            const scalar refDist(subDict.lookupOrDefault<scalar>("referenceDist", 1.0));

            vector refDire(vector::zero);
            refDire.component(vector::Y) = ::sin(refAngl*constant::mathematical::twoPi/360);
            refDire.component(vector::Z) = ::cos(refAngl*constant::mathematical::twoPi/360);

            if (refDist <= 0)
            {
                Info<<"reference distance of the " << fieldName << " field is no larger than zero (ERROR)" << endl;
            }

            const scalarField dirCmpt((refDire&fCentres)/refDist);

            vector eigenVal(vector::zero);
            tensor eigenVec(tensor::zero);

            if (referenceValue.component(symmTensor::XY)==0&&referenceValue.component(symmTensor::XZ)==0&&referenceValue.component(symmTensor::YZ)==0)
            {
                eigenVal = vector(referenceValue.component(symmTensor::XX),referenceValue.component(symmTensor::YY),referenceValue.component(symmTensor::ZZ));
                eigenVec.component(tensor::XX) = 1.0;
                eigenVec.component(tensor::YY) = 1.0;
                eigenVec.component(tensor::ZZ) = 1.0;
            }
            else
            {
                eigenVal = eigenValues(referenceValue);
                eigenVec = eigenVectors(referenceValue);
            }

            vector a1(eigenVec.component(tensor::XX),eigenVec.component(tensor::XY),eigenVec.component(tensor::XZ));
            vector a2(eigenVec.component(tensor::YX),eigenVec.component(tensor::YY),eigenVec.component(tensor::YZ));
            vector a3(eigenVec.component(tensor::ZX),eigenVec.component(tensor::ZY),eigenVec.component(tensor::ZZ));

            symmTensor M1(a1[0]*a1[0],a1[0]*a1[1],a1[0]*a1[2],a1[1]*a1[1],a1[1]*a1[2],a1[2]*a1[2]);
            symmTensor M2(a2[0]*a2[0],a2[0]*a2[1],a2[0]*a2[2],a2[1]*a2[1],a2[1]*a2[2],a2[2]*a2[2]);
            symmTensor M3(a3[0]*a3[0],a3[0]*a3[1],a3[0]*a3[2],a3[1]*a3[1],a3[1]*a3[2],a3[2]*a3[2]);

            M1 = M1*eigenVal[0];
            M2 = M2*eigenVal[1];
            M3 = M3*eigenVal[2];

            vector alpha(subDict.lookupOrDefault<vector>("alpha", vector::zero));

            forAll(value, label)
            {
                vector powCoef(cmptPow(vector::one*fabs(dirCmpt[label]), alpha));
                value[label] = powCoef[0]*M1+powCoef[1]*M2+powCoef[2]*M3;
            }
        }
        else if (profile == "linear")
        {
            vectorField fCentres(patchCoord.localPosition(this->patch().Cf()));

            const scalar refAngl(subDict.lookupOrDefault<scalar>("referenceAngl", 0.0));
            const scalar refDist(subDict.lookupOrDefault<scalar>("referenceDist", 1.0));

            vector refDire(vector::zero);
            refDire.component(vector::Y) = ::sin(refAngl*constant::mathematical::twoPi/360);
            refDire.component(vector::Z) = ::cos(refAngl*constant::mathematical::twoPi/360);

            if (refDist <= 0)
            {
                Info<<"reference distance of the " << fieldName << " field is no larger than zero (ERROR)" << endl;
            }

            const scalarField dirCmpt((refDire&fCentres)/refDist);

            vector eigenVal(vector::zero);
            tensor eigenVec(tensor::zero);

            if (referenceValue.component(symmTensor::XY)==0&&referenceValue.component(symmTensor::XZ)==0&&referenceValue.component(symmTensor::YZ)==0)
            {
                eigenVal = vector(referenceValue.component(symmTensor::XX),referenceValue.component(symmTensor::YY),referenceValue.component(symmTensor::ZZ));
                eigenVec.component(tensor::XX) = 1.0;
                eigenVec.component(tensor::YY) = 1.0;
                eigenVec.component(tensor::ZZ) = 1.0;
            }
            else
            {
                eigenVal = eigenValues(referenceValue);
                eigenVec = eigenVectors(referenceValue);
            }

            vector a1(eigenVec.component(tensor::XX),eigenVec.component(tensor::XY),eigenVec.component(tensor::XZ));
            vector a2(eigenVec.component(tensor::YX),eigenVec.component(tensor::YY),eigenVec.component(tensor::YZ));
            vector a3(eigenVec.component(tensor::ZX),eigenVec.component(tensor::ZY),eigenVec.component(tensor::ZZ));

            symmTensor M1(a1[0]*a1[0],a1[0]*a1[1],a1[0]*a1[2],a1[1]*a1[1],a1[1]*a1[2],a1[2]*a1[2]);
            symmTensor M2(a2[0]*a2[0],a2[0]*a2[1],a2[0]*a2[2],a2[1]*a2[1],a2[1]*a2[2],a2[2]*a2[2]);
            symmTensor M3(a3[0]*a3[0],a3[0]*a3[1],a3[0]*a3[2],a3[1]*a3[1],a3[1]*a3[2],a3[2]*a3[2]);

            M1 = M1*eigenVal[0];
            M2 = M2*eigenVal[1];
            M3 = M3*eigenVal[2];

            vector alpha(subDict.lookupOrDefault<vector>("alpha", vector::zero));

            forAll(value, label)
            {
                vector powCoef(alpha*(dirCmpt[label]-1.0)+vector::one);
                value[label] = powCoef[0]*M1+powCoef[1]*M2+powCoef[2]*M3;
            }
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

    tmp<symmTensorField> tFld(new symmTensorField(value));

    Info<< "Turbulent Mean patch " << this->patch().name()
        << ": calculating field " << fieldName
        << " from " << dict.name() << endl;

    return tFld;
}

// ************************************************************************* //
