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

// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * //

Foam::eddy::eddy()
:
    shape_("tent"),
    patchFaceI_(-1),
    position0_(vector::zero),
    x_(0),
    convectVelocity_(0),
    sigma_u_(vector::zero),
    sigma_v_(vector::zero),
    sigma_w_(vector::zero),
    epsilon_(vector::zero),
    Lund_(tensor::zero)
{}


Foam::eddy::eddy
(
    const word shape,
    const label patchFaceI,
    const point& position0,
    const scalar x,
    const scalar convectVelocity,
    const vector sigma_u,
    const vector sigma_v,
    const vector sigma_w,
    const symmTensor& R,
    Random& rndGen
)
:
    shape_(shape),
    patchFaceI_(patchFaceI),
    position0_(position0),
    x_(x),
    convectVelocity_(convectVelocity),
    sigma_u_(sigma_u),
    sigma_v_(sigma_v),
    sigma_w_(sigma_w),
    epsilon_(epsilon(rndGen)),
    Lund_(tensor::zero)
{
    Lund_.replace(tensor::XX, sqrt(R.component(symmTensor::XX)));
    Lund_.replace(tensor::YX, R.component(symmTensor::XY)/Lund_.component(tensor::XX));
    Lund_.replace(tensor::ZX, R.component(symmTensor::XZ)/Lund_.component(tensor::XX));
    Lund_.replace(tensor::YY, sqrt(R.component(symmTensor::YY)-sqr(Lund_.component(tensor::YX))));
    Lund_.replace(tensor::ZY, (R.component(symmTensor::YZ) - Lund_.component(tensor::YX)*Lund_.component(tensor::ZX) )/Lund_.component(tensor::YY));
    Lund_.replace(tensor::ZZ, sqrt(R.component(symmTensor::ZZ) - sqr(Lund_.component(tensor::ZX))-sqr(Lund_.component(tensor::ZY))));
}


Foam::eddy::eddy(const eddy& e)
:
    shape_(e.shape_),
    patchFaceI_(e.patchFaceI_),
    position0_(e.position0_),
    x_(e.x_),
    convectVelocity_(e.convectVelocity_),
    sigma_u_(e.sigma_u_),
    sigma_v_(e.sigma_v_),
    sigma_w_(e.sigma_w_),
    epsilon_(e.epsilon_),
    Lund_(e.Lund_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vector Foam::eddy::uDash(const point& xp, const vector& n) const
{
    // Relative position inside eddy (global system)
    const vector r_u = cmptDivide(xp-position(n), sigma_u_);
    const vector r_v = cmptDivide(xp-position(n), sigma_v_);
    const vector r_w = cmptDivide(xp-position(n), sigma_w_);

    vector uDash = vector::zero;

    uDash.component(0) = epsilon_.component(0)*shape(r_u.component(0))*shape(r_u.component(1))*shape(r_u.component(2))
                       / sqrt(cmptProduct(sigma_u_));

    uDash.component(1) = epsilon_.component(1)*shape(r_v.component(0))*shape(r_v.component(1))*shape(r_v.component(2))
                       / sqrt(cmptProduct(sigma_v_));

    uDash.component(2) = epsilon_.component(2)*shape(r_w.component(0))*shape(r_w.component(1))*shape(r_w.component(2))
                       / sqrt(cmptProduct(sigma_w_));

    uDash = Lund_&uDash;

    return uDash;
}


Foam::scalar Foam::eddy::shape(const scalar& x) const
{
    if (fabs(x) >= 1)
    {
        return 0;
    }
    else
    {
        if (shape_ == "tent")
        {
            return sqrt(1.5)*(1.0-fabs(x));
        }
        else if (shape_ == "step")
        {
            return sqrt(0.5);
        }
        else if (shape_ == "gaussian")
        {
            return 1.3010*exp(-4.5*x*x);
        }
        else
        {
            Info << "shape" << shape_ 
                 << "does not exist (ERROR)" << endl;
            return 0;
        }
    }
}

// ************************************************************************* //
