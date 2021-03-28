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
#include "mathematicalConstants.H"
#include "UList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::vorton::debug = 0;

// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * //

Foam::vorton::vorton()
:
    type_("typeR"),
    patchFaceI_(-1),
    position0_(vector::zero),
    x_(0),
    sigma_(vector::zero),
    gamma_(vector::zero),
    Rpg_(tensor::I)
{}


Foam::vorton::vorton
(
    const word type,
    const label patchFaceI,
    const point& position0,
    const scalar x,
    const vector L0,
    const symmTensor& R0,
    Random& rndGen
)
:
    type_(type),
    patchFaceI_(patchFaceI),
    position0_(position0),
    x_(x),
    sigma_(vector::zero),
    gamma_(vector::zero),
    Rpg_(tensor::I)
{
    const scalar pi = constant::mathematical::pi;

    // Principal stresses - eigenvalues returned in ascending order
    vector R;

    if (R0.component(symmTensor::XY)==0&&R0.component(symmTensor::XZ)==0&&R0.component(symmTensor::YZ)==0)
    {
        R = vector(R0.component(symmTensor::XX),R0.component(symmTensor::YY),R0.component(symmTensor::ZZ));
    }
    else
    {
        R = eigenValues(R0);
    }

    // Eddy rotation from principal-to-global axes
    // - given by the 3 eigenvectors of the Reynold stress tensor as rows in
    //   the result tensor (transposed transformation tensor)
    // - returned in ascending eigenvalue order
    Rpg_ = eigenVectors(R0, R);

    if (debug)
    {
        Pout<< "Rpg & R & Rpg.T(): " << (Rpg_ & R0 & Rpg_.T()) << endl;
    }

    const tensor Rpg2 = cmptMultiply(Rpg_.T(), Rpg_.T());

    const vector ex = cmptMultiply(Rpg2.x(), R)/R0.xx();
    const vector ey = cmptMultiply(Rpg2.y(), R)/R0.yy();
    const vector ez = cmptMultiply(Rpg2.z(), R)/R0.zz();

    const tensor E = tensor(ex, ey, ez);

    vector L = inv(E)&L0;

    if (type == "typeR" || L.x() < 0.0)
    {
        L.x() = (L.y()*L.z()*sqrt(R.x()))/(L.y()*sqrt(R.z())+L.z()*sqrt(R.y()));
    }

    sigma_ = L/sqrt(pi);

    scalar f1 = 2.0*sqrt(R.x()*L.y()*L.z()/L.x())/pi;
    scalar f2 = 2.0*sqrt(R.y()*L.x()*L.z()/L.y())/pi;
    scalar f3 = 2.0*sqrt(R.z()*L.x()*L.y()/L.z())/pi;

    if (type == "typeR")
    {
        if (f1 + f3 != f2)
        {
            if (f1 < f3)
            {
                f1 = -f1;
            }
        }

        vector gamma1 = vector::zero;
        vector gamma2 = vector::zero;

        gamma1.z() = 1.0;
        gamma2.z() = 1.0;

        gamma1.x() = (sqr(sigma_.z())+f2)/sqr(sigma_.x());
        gamma2.x() = (sqr(sigma_.z())-f2)/sqr(sigma_.x());

        gamma1.y() = (sqr(sigma_.z())+f1)/sqr(sigma_.y());
        gamma2.y() = (sqr(sigma_.z())-f1)/sqr(sigma_.y());

        if (mag(gamma1-vector::one) < mag(gamma2-vector::one))
        {
            gamma_ = gamma1;
        }
        else
        {
            gamma_ = gamma2;
        }
    }
    else if (type == "typeL")
    {
        tensor M = tensor::zero;

        M.xy() = sqr(sigma_.y());
        M.xz() =-sqr(sigma_.z());

        M.yx() = sqr(sigma_.x());
        M.yz() =-sqr(sigma_.z());

        M.zx() = sqr(sigma_.x());
        M.zy() =-sqr(sigma_.y());

        tensor pinvM = pinv(M);

        List<vector> sign(8, vector::zero);

        label indi = 0;

        for (label i=-1; i<=1; i=i+2)
        {
            for (label j=-1; j<=1; j=j+2)
            {
                for (label k=-1; k<=1; k=k+2)
                {
                    sign[indi] = vector(i,j,k);
                    indi = indi+1;
                }
            }
        }

        List<vector> Gamma(8, vector::zero);
        List<vector> RGamma(8, vector::zero);
        List<scalar> RGammaValue(8);

        forAll(Gamma, label)
        {
            vector f = cmptMultiply(sign[label], vector(f1,f2,f3));
            Gamma[label] = pinvM&f;

            RGamma[label].x() = 0.25*pi*sqrt(pi)*sigma_.x()*sqr(Gamma[label].y()*sqr(sigma_.y())-Gamma[label].z()*sqr(sigma_.z()))/sigma_.y()/sigma_.z();
            RGamma[label].y() = 0.25*pi*sqrt(pi)*sigma_.y()*sqr(Gamma[label].x()*sqr(sigma_.x())-Gamma[label].z()*sqr(sigma_.z()))/sigma_.x()/sigma_.z();
            RGamma[label].z() = 0.25*pi*sqrt(pi)*sigma_.z()*sqr(Gamma[label].x()*sqr(sigma_.x())-Gamma[label].y()*sqr(sigma_.y()))/sigma_.x()/sigma_.y();

            RGamma[label] = RGamma[label]-R;
            RGammaValue[label] = mag(RGamma[label]);
        }

        List<label> order;
        sortedOrder(RGammaValue, order);

        gamma_ = Gamma[order[0]];     
    }
    else
    {
        Info << "vorton type " << type << " does not exist (ERROR)" << endl;
    }

    gamma_ *= epsi(rndGen);
}

Foam::vorton::vorton
(
    const word type,
    const label patchFaceI,
    const point& position0,
    const scalar x,
    const symmTensor& R0,
    const vector sigma,
    const vector gamma
)
:
    type_(type),
    patchFaceI_(patchFaceI),
    position0_(position0),
    x_(x),
    sigma_(sigma),
    gamma_(gamma),
    Rpg_(tensor::I)
{
    // Principal stresses - eigenvalues returned in ascending order
    const vector R = eigenValues(R0);

    // Eddy rotation from principal-to-global axes
    // - given by the 3 eigenvectors of the Reynold stress tensor as rows in
    //   the result tensor (transposed transformation tensor)
    // - returned in ascending eigenvalue order
    Rpg_ = eigenVectors(R0, R);
}

Foam::vorton::vorton(const vorton& v)
:
    type_(v.type_),
    patchFaceI_(v.patchFaceI_),
    position0_(v.position0_),
    x_(v.x_),
    sigma_(v.sigma_),
    gamma_(v.gamma_),
    Rpg_(v.Rpg_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tensor Foam::vorton::pinv(const tensor& A) const
{
    const tensor ATA = A.T()&A;

    vector lambda = eigenValues(ATA);
    tensor V = eigenVectors(ATA, lambda);

    V = tensor(V.z(), V.y(), V.x()).T();

    tensor S = tensor::zero;
    tensor invS = tensor::zero;

    S.xx() = sqrt(lambda.z());
    S.yy() = sqrt(lambda.y());

    invS.xx() = 1.0/S.xx();
    invS.yy() = 1.0/S.yy();

    tensor U = (A&V)&invS;

    return (V&invS)&U.T();
}


Foam::vectorField Foam::vorton::uDash(const pointField& xp, const vector& n) const
{
    vectorField uDash(xp.size(), vector::zero);

    forAll(uDash, label)
    {
        const vector x = Rpg()&(xp[label]-position(n));
        const vector r = cmptDivide(x, sigma());
        const scalar c = exp(-0.5*cmptSum(cmptMultiply(r, r)));

        uDash[label].x() = c*(gamma().y()/sqr(sigma().z())-gamma().z()/sqr(sigma().y()))*x.y()*x.z();
        uDash[label].y() = c*(gamma().z()/sqr(sigma().x())-gamma().x()/sqr(sigma().z()))*x.x()*x.z();
        uDash[label].z() = c*(gamma().x()/sqr(sigma().y())-gamma().y()/sqr(sigma().x()))*x.x()*x.y();
    }

    return (Rpg().T()&uDash);
}

// ************************************************************************* //
