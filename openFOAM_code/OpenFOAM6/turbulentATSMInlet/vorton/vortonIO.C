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

#include "vorton.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vorton::vorton(Istream& is)
:
    type_(is),
    patchFaceI_(readLabel(is)),
    position0_(is),
    x_(readScalar(is)),
    sigma_(is),
    gamma_(is),
    Rpg_(is)
{
    is.check(FUNCTION_NAME);
}


// * * * * * * * * * * * * * * * * * Operators * * * * * * * * * * * * * * * //

void Foam::vorton::operator=(const vorton& v)
{
    type_ = v.type_;
    patchFaceI_ = v.patchFaceI_;
    position0_ = v.position0_;
    x_ = v.x_;
    sigma_ = v.sigma_;
    gamma_ = v.gamma_;
    Rpg_ = v.Rpg_;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, vorton& v)
{
    is.check(FUNCTION_NAME);

    is  >> v.type_
        >> v.patchFaceI_
        >> v.position0_
        >> v.x_
        >> v.sigma_
        >> v.gamma_
        >> v.Rpg_;

    is.check(FUNCTION_NAME);
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const vorton& v)
{
    os.check(FUNCTION_NAME);

    os  << v.type_ << token::SPACE
        << v.patchFaceI_ << token::SPACE
        << v.position0_ << token::SPACE
        << v.x_ << token::SPACE
        << v.sigma_ << token::SPACE
        << v.gamma_ << token::SPACE
        << v.Rpg_;

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
