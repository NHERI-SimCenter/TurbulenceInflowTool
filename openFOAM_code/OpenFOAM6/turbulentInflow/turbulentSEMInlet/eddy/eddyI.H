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

Foam::scalar Foam::eddy::epsi(Random& rndGen) const
{
    // Random number with zero mean and unit variance
    return rndGen.sample01<scalar>() > 0.5 ? 1 : -1;
}


inline Foam::word Foam::eddy::type() const
{
    return type_;
}


inline Foam::label Foam::eddy::patchFaceI() const
{
    return patchFaceI_;
}


inline const Foam::point& Foam::eddy::position0() const
{
    return position0_;
}


inline Foam::scalar Foam::eddy::x() const
{
    return x_;
}


inline const Foam::tensor& Foam::eddy::sigma() const
{
    return sigma_;
}


inline const Foam::vector& Foam::eddy::gamma() const
{
    return gamma_;
}


inline const Foam::tensor& Foam::eddy::Lund() const
{
    return Lund_;
}


inline Foam::point Foam::eddy::position(const vector& n) const
{
    return position0_ + n*x_;
}


Foam::vector Foam::eddy::epsilon(Random& rndGen) const
{
    return vector(epsi(rndGen), epsi(rndGen), epsi(rndGen));
}

inline Foam::scalar Foam::eddy::volume() const
{
    const vector sigmaMin = vector
    (
        cmptMin(sigma().x()), 
        cmptMin(sigma().y()), 
        cmptMin(sigma().z())
    );

    return cmptProduct(sigmaMin);
}


void Foam::eddy::move(const scalar dx)
{
    x_ += dx;
}


Foam::boundBox Foam::eddy::bounds(const bool global) const
{
    boundBox bb;

    const vector sigmaMax = vector
    (
        cmptMax(sigma().x()), 
        cmptMax(sigma().y()), 
        cmptMax(sigma().z())
    );

    bb.min() = -sigmaMax;
    bb.max() =  sigmaMax;

    return bb;
}
// ************************************************************************* //
