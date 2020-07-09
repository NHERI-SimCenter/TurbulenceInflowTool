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

#include "eddy.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eddy::eddy(Istream& is)
:
    type_(is),
    patchFaceI_(readLabel(is)),
    position0_(is),
    x_(readScalar(is)),
    sigma_(is),
    gamma_(is),
    Lund_(is)
{
    is.check(FUNCTION_NAME);
}


// * * * * * * * * * * * * * * * * * Operators * * * * * * * * * * * * * * * //

void Foam::eddy::operator=(const eddy& e)
{
    type_ = e.type_;
    patchFaceI_ = e.patchFaceI_;
    position0_ = e.position0_;
    x_ = e.x_;
    sigma_ = e.sigma_;
    gamma_ = e.gamma_;
    Lund_ = e.Lund_;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, eddy& e)
{
    is.check(FUNCTION_NAME);

    is  >> e.type_
        >> e.patchFaceI_
        >> e.position0_
        >> e.x_
        >> e.sigma_
        >> e.gamma_
        >> e.Lund_;

    is.check(FUNCTION_NAME);
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const eddy& e)
{
    os.check(FUNCTION_NAME);

    os  << e.type_ << token::SPACE
        << e.patchFaceI_ << token::SPACE
        << e.position0_ << token::SPACE
        << e.x_ << token::SPACE
        << e.sigma_ << token::SPACE
        << e.gamma_ << token::SPACE
        << e.Lund_;

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
