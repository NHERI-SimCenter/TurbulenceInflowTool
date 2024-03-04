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
#include "mathematicalConstants.H"
#include "UList.H"

// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * //

Foam::eddy::eddy()
:
    type_("gaussian"),
    patchFaceI_(-1),
    position0_(vector::zero),
    x_(0),
    sigma_(tensor::zero),
    gamma_(vector::zero),
    Lund_(tensor::I)
{}


Foam::eddy::eddy
(
    const word type,
    const label patchFaceI,
    const point& position0,
    const scalar x,
    const tensor L,
    const symmTensor& R,
    Random& rndGen
)
:
    type_(type),
    patchFaceI_(patchFaceI),
    position0_(position0),
    x_(x),
    sigma_(tensor::zero),
    gamma_(epsilon(rndGen)),
    Lund_(tensor::I)
{
    Lund_.xx() = sqrt(R.xx());
    Lund_.yx() = R.xy()/Lund_.xx();
    Lund_.zx() = R.xz()/Lund_.xx();
    Lund_.yy() = sqrt(R.yy()-sqr(Lund_.yx()));
    Lund_.zy() = (R.yz() - Lund_.yx()*Lund_.zx())/Lund_.yy();
    Lund_.zz() = sqrt(R.zz() - sqr(Lund_.zx()) - sqr(Lund_.zy()));

    const tensor Lund2 = cmptMultiply(Lund_, Lund_);

    const vector ex = Lund2.x()/R.xx();
    const vector ey = Lund2.y()/R.yy();
    const vector ez = Lund2.z()/R.zz();

    const tensor E = tensor(ex, ey, ez);

    sigma_ = 2.0*(inv(E)&L);
}

Foam::eddy::eddy
(
    const word type,
    const label patchFaceI,
    const point& position0,
    const scalar x,
    const symmTensor& R,
    const tensor sigma,
    const vector gamma
)
:
    type_(type),
    patchFaceI_(patchFaceI),
    position0_(position0),
    x_(x),
    sigma_(sigma),
    gamma_(gamma),
    Lund_(tensor::I)
{
    Lund_.xx() = sqrt(R.xx());
    Lund_.yx() = R.xy()/Lund_.xx();
    Lund_.zx() = R.xz()/Lund_.xx();
    Lund_.yy() = sqrt(R.yy()-sqr(Lund_.yx()));
    Lund_.zy() = (R.yz() - Lund_.yx()*Lund_.zx())/Lund_.yy();
    Lund_.zz() = sqrt(R.zz() - sqr(Lund_.zx()) - sqr(Lund_.zy()));
}

Foam::eddy::eddy(const eddy& e)
:
    type_(e.type_),
    patchFaceI_(e.patchFaceI_),
    position0_(e.position0_),
    x_(e.x_),
    sigma_(e.sigma_),
    gamma_(e.gamma_),
    Lund_(e.Lund_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vectorField Foam::eddy::uDash(const pointField& xp, const vector& n) const
{
    const scalar pi = constant::mathematical::pi;

    vectorField uDash(xp.size(), vector::zero);

    if (type() == "gaussian")
    {
        const scalar cx = gamma().x()*sqrt(8.0)/sqrt(cmptProduct(sigma().x()));
        const scalar cy = gamma().y()*sqrt(8.0)/sqrt(cmptProduct(sigma().y()));
        const scalar cz = gamma().z()*sqrt(8.0)/sqrt(cmptProduct(sigma().z()));

        forAll(uDash, label)
        {
            const vector rx = cmptDivide(xp[label]-position(n), sigma().x());
            const vector ry = cmptDivide(xp[label]-position(n), sigma().y());
            const vector rz = cmptDivide(xp[label]-position(n), sigma().z());

            const scalar rx2 = cmptSum(cmptMultiply(rx, rx));
            const scalar ry2 = cmptSum(cmptMultiply(ry, ry));
            const scalar rz2 = cmptSum(cmptMultiply(rz, rz));

            if (rx2 <= 1.0)
            {
                uDash[label].x() = cx*exp(-2.0*pi*rx2);
            }

            if (ry2 <= 1.0)
            {
                uDash[label].y() = cy*exp(-2.0*pi*ry2);
            }

            if (rz2 <= 1.0)
            {
                uDash[label].z() = cz*exp(-2.0*pi*rz2);
            }
        }
    }
    else if (type() == "tent")
    {
        const scalar cx = gamma().x()*sqrt(3.375)/sqrt(cmptProduct(sigma().x()));
        const scalar cy = gamma().y()*sqrt(3.375)/sqrt(cmptProduct(sigma().y()));
        const scalar cz = gamma().z()*sqrt(3.375)/sqrt(cmptProduct(sigma().z()));

        forAll(uDash, label)
        {
            const vector rx = cmptDivide(cmptMag(xp[label]-position(n)), sigma().x());
            const vector ry = cmptDivide(cmptMag(xp[label]-position(n)), sigma().y());
            const vector rz = cmptDivide(cmptMag(xp[label]-position(n)), sigma().z());

            if (cmptMax(rx) <= 1.0)
            {
                uDash[label].x() = cx*cmptProduct(vector::one-rx);
            }

            if (cmptMax(ry) <= 1.0)
            {
                uDash[label].y() = cy*cmptProduct(vector::one-ry);
            }

            if (cmptMax(rz) <= 1.0)
            {
                uDash[label].z() = cz*cmptProduct(vector::one-rz);
            }
        }
    }
    else if (type() == "step")
    {
        const scalar cx = gamma().x()*sqrt(0.125)/sqrt(cmptProduct(sigma().x()));
        const scalar cy = gamma().y()*sqrt(0.125)/sqrt(cmptProduct(sigma().y()));
        const scalar cz = gamma().z()*sqrt(0.125)/sqrt(cmptProduct(sigma().z()));

        forAll(uDash, label)
        {
            const vector rx = cmptDivide(cmptMag(xp[label]-position(n)), sigma().x());
            const vector ry = cmptDivide(cmptMag(xp[label]-position(n)), sigma().y());
            const vector rz = cmptDivide(cmptMag(xp[label]-position(n)), sigma().z());

            if (cmptMax(rx) <= 1.0)
            {
                uDash[label].x() = cx;
            }

            if (cmptMax(ry) <= 1.0)
            {
                uDash[label].y() = cy;
            }

            if (cmptMax(rz) <= 1.0)
            {
                uDash[label].z() = cz;
            }
        }
    }
    else
    {
        Info << "eddy type: " << type_
             << "does not exist (ERROR)" << endl;
    }

    return (Lund()&uDash);
}

// ************************************************************************* //
