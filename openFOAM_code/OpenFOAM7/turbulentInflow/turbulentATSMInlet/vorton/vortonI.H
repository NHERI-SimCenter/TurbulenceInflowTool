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

#include "mathematicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::scalar Foam::vorton::epsi(Random& rndGen) const
{
    // Random number with zero mean and unit variance
    return rndGen.sample01<scalar>() > 0.5 ? 1 : -1;
}


inline Foam::word Foam::vorton::type() const
{
    return type_;
}


inline Foam::label Foam::vorton::patchFaceI() const
{
    return patchFaceI_;
}


inline const Foam::point& Foam::vorton::position0() const
{
    return position0_;
}


inline Foam::scalar Foam::vorton::x() const
{
    return x_;
}


inline const Foam::vector& Foam::vorton::sigma() const
{
    return sigma_;
}


inline const Foam::vector& Foam::vorton::gamma() const
{
    return gamma_;
}


inline const Foam::tensor& Foam::vorton::Rpg() const
{
    return Rpg_;
}


inline Foam::point Foam::vorton::position(const vector& n) const
{
    return position0_ + n*x_;
}


inline Foam::scalar Foam::vorton::volume() const
{
    return cmptProduct(sigma());
}


void Foam::vorton::move(const scalar dx)
{
    x_ += dx;
}


Foam::boundBox Foam::vorton::bounds(const bool global) const
{
    boundBox bb;

    bb.min() = -3.0*sigma();
    bb.max() =  3.0*sigma();

    return bb;
}
// ************************************************************************* //
