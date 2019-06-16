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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::scalar Foam::eddy::epsi(Random& rndGen) const
{
    // Random number with zero mean and unit variance
    return rndGen.scalar01() > 0.5 ? 1 : -1;
}


inline const Foam::word Foam::eddy::shape() const
{
    return shape_;
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


inline const Foam::scalar& Foam::eddy::convectVelocity() const
{
    return convectVelocity_;
}


inline const Foam::vector& Foam::eddy::sigma_u() const
{
    return sigma_u_;
}


inline const Foam::vector& Foam::eddy::sigma_v() const
{
    return sigma_v_;
}


inline const Foam::vector& Foam::eddy::sigma_w() const
{
    return sigma_w_;
}


inline const Foam::vector& Foam::eddy::epsilon() const
{
    return epsilon_;
}


inline const Foam::tensor& Foam::eddy::Lund() const
{
    return Lund_;
}


inline Foam::point Foam::eddy::position(const vector& n) const
{
    return position0_ + n*x_;
}


inline Foam::vector Foam::eddy::epsilon(Random& rndGen) const
{
    return vector(epsi(rndGen), epsi(rndGen), epsi(rndGen));
}


inline Foam::scalar Foam::eddy::volume() const
{
    vector sigma = (sigma_u_+sigma_v_+sigma_w_)/3.0;

    return scalar(sigma.component(0)*sigma.component(1)*sigma.component(2));
}

void Foam::eddy::move(const scalar dt)
{
    x_ += convectVelocity_*dt;
}

Foam::boundBox Foam::eddy::bounds() const
{
    boundBox bb;

    vector sigma_x = vector(sigma_u_.component(0),sigma_v_.component(0),sigma_w_.component(0));
    vector sigma_y = vector(sigma_u_.component(1),sigma_v_.component(1),sigma_w_.component(1));
    vector sigma_z = vector(sigma_u_.component(2),sigma_v_.component(2),sigma_w_.component(2));

    vector sigma = vector(cmptMax(sigma_x),cmptMax(sigma_y),cmptMax(sigma_z));

    bb.min() = -sigma;
    bb.max() = sigma;

    return bb;
}

// ************************************************************************* //
