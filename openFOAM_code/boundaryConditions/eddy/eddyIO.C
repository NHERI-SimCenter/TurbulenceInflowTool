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
    shape_(is),
    patchFaceI_(readLabel(is)),
    position0_(is),
    x_(readScalar(is)),
    convectVelocity_(readScalar(is)),
    sigma_u_(is),
    sigma_v_(is),
    sigma_w_(is),
    epsilon_(is),
    Lund_(is)
{
    is.check(FUNCTION_NAME);
}


// * * * * * * * * * * * * * * * * * Operators * * * * * * * * * * * * * * * //

void Foam::eddy::operator=(const eddy& e)
{
    shape_ = e.shape_;
    patchFaceI_ = e.patchFaceI_;
    position0_ = e.position0_;
    x_ = e.x_;
    convectVelocity_ = e.convectVelocity_;
    sigma_u_ = e.sigma_u_;
    sigma_v_ = e.sigma_v_;
    sigma_w_ = e.sigma_w_;
    epsilon_ = e.epsilon_;
    Lund_ = e.Lund_;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, eddy& e)
{
    is.check(FUNCTION_NAME);

    is  >> e.shape_
        >> e.patchFaceI_
        >> e.position0_
        >> e.x_
        >> e.convectVelocity_
        >> e.sigma_u_
        >> e.sigma_v_
        >> e.sigma_w_
        >> e.epsilon_
        >> e.Lund_;

    is.check(FUNCTION_NAME);
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const eddy& e)
{
    os.check(FUNCTION_NAME);

    os  << e.shape_ << token::SPACE
        << e.patchFaceI_ << token::SPACE
        << e.position0_ << token::SPACE
        << e.x_ << token::SPACE
        << e.convectVelocity_ << token::SPACE
        << e.sigma_u_ << token::SPACE
        << e.sigma_v_ << token::SPACE
        << e.sigma_w_ << token::SPACE
        << e.epsilon_ << token::SPACE
        << e.Lund_;

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
