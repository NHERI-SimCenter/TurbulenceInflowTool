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

#include "dfeddy.H"
#include "mathematicalConstants.H"
#include "UList.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::label Foam::dfeddy::Gamma2Values[] = {1, 2, 3, 4, 5, 6, 7, 8};
Foam::UList<Foam::label> Foam::dfeddy::Gamma2(&Gamma2Values[0], 8);
int Foam::dfeddy::debug = 0;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::dfeddy::setScales
(
    const scalar sigmaX,
    const label gamma2,
    const vector& e,
    const vector& lambda,
    vector& sigma,
    vector& alpha
) const
{
    // Static array of gamma^2 vs c2 coefficient
    static const scalar gamma2VsC2[8] =
        {2, 1.875, 1.737, 1.75, 0.91, 0.825, 0.806, 1.5};

    scalar gamma = Foam::sqrt(scalar(gamma2));

    // c2 coefficient retrieved from array
    scalar c2 = gamma2VsC2[gamma2 - 1];

    // Length scale in largest eigenvalue direction
    label d1 = dir1_;
    label d2 = (d1 + 1) % 3;
    label d3 = (d1 + 2) % 3;

    sigma[d1] = sigmaX;

    // Note: sigma_average = 1/3*(sigma_x + sigma_y + sigma_z)
    // Substituting for sigma_y = sigma_x/gamma and sigma_z = sigma_y
    //sigma[d1] = 3*sigmaX/(1 + 2/gamma);
    // Other length scales equal, as function of major axis length and gamma
    sigma[d2] = sigma[d1]/gamma;
    sigma[d3] = sigma[d2];

    vector sigma2 = cmptMultiply(sigma, sigma);
    scalar slos2 = cmptSum(cmptDivide(lambda, sigma2));

    bool ok = true;

    for (label beta = 0; beta < 3; ++beta)
    {
        scalar x = slos2 - 2*lambda[beta]/sigma2[beta];

        if (x < 0)
        {
            alpha[beta] = 0;
            ok = false;
        }
        else
        {
            alpha[beta] = e[beta]*sqrt(x/(2*c2));
        }
    }

    if (debug > 1)
    {
        Pout<< "c2:" << c2
            << ", gamma2:" << gamma2
            << ", gamma:" << gamma
            << ", lambda:" << lambda
            << ", sigma2: " << sigma2
            << ", slos2: " << slos2
            << ", sigmaX:" << sigmaX
            << ", sigma:" << sigma
            << ", alpha:" << alpha
            << endl;
    }

    return ok;
}


// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * //

Foam::dfeddy::dfeddy()
:
    patchFaceI_(-1),
    position0_(vector::zero),
    x_(0),
    sigma_(vector::zero),
    alpha_(vector::zero),
    Rpg_(tensor::I),
    c1_(-1),
    dir1_(0)
{}


Foam::dfeddy::dfeddy
(
    const label patchFaceI,
    const point& position0,
    const scalar x,
    const scalar sigmaX,
    const symmTensor& R,
    Random& rndGen
)
:
    patchFaceI_(patchFaceI),
    position0_(position0),
    x_(x),
    sigma_(vector::zero),
    alpha_(vector::zero),
    Rpg_(tensor::I),
    c1_(-1),
    dir1_(0)
{
    // Principal stresses - eigenvalues returned in ascending order
    vector lambda;

    if (R.component(symmTensor::XY)==0&&R.component(symmTensor::XZ)==0&&R.component(symmTensor::YZ)==0)
    {
        lambda = vector(R.component(symmTensor::XX),R.component(symmTensor::YY),R.component(symmTensor::ZZ));
    }
    else
    {
        lambda = eigenValues(R);
    }

    // Eddy rotation from principal-to-global axes
    // - given by the 3 eigenvectors of the Reynold stress tensor as rows in
    //   the result tensor (transposed transformation tensor)
    // - returned in ascending eigenvalue order
    Rpg_ = eigenVectors(R, lambda).T();

    if (debug)
    {
        // Global->Principal transform = Rgp = Rpg.T()
        // Rgp & R & Rgp.T() should have eigenvalues on its diagonal and
        // zeros for all other components
        Pout<< "Rpg.T() & R & Rpg: " << (Rpg_.T() & R & Rpg_) << endl;
    }

    // Set the dfeddy orientation to position of max eigenvalue
    // (direction of dfeddy major axis, sigma_x in reference)
    dir1_ = 2;

    // Random vector of 1's and -1's
    const vector e(epsilon(rndGen));

    // Set intensities and length scales
    bool found = false;
    forAll(Gamma2, i)
    {
        // Random length scale ratio, gamma = sigmax/sigmay = sigmax/sigmaz
        // - using gamma^2 to ease lookup of c2 coefficient
        label g2 = Gamma2[i];

        if (setScales(sigmaX, g2, e, lambda, sigma_, alpha_))
        {
            found = true;
            break;
        }
    }

    // Normalisation coefficient (eq. 11)
    // Note: sqrt(10*V)/sqrt(nEddy) applied outside when computing uDash
    c1_ = cmptAv(sigma_)/cmptProduct(sigma_)*cmptMin(sigma_);

    if (found)
    {
        // Shuffle the gamma^2 values
        shuffle(Gamma2);
    }
    else
    {
        if (debug)
        {
            // If not found typically means that the stress has a repeated
            // eigenvalue/not covered by the selection of Gamma values, e.g.
            // as seen by range of applicability on Lumley diagram
            WarningInFunction
                << "Unable to set dfeddy intensity for dfeddy: " << *this
                << endl;
        }

        // Remove the influence of this dfeddy/indicate that its initialisation
        // failed
        patchFaceI_ = -1;
    }
}


Foam::dfeddy::dfeddy
(
    const label patchFaceI,
    const point& position0,
    const scalar x,
    const symmTensor& R,
    const vector sigma,
    const vector alpha
)
:
    patchFaceI_(patchFaceI),
    position0_(position0),
    x_(x),
    sigma_(sigma),
    alpha_(alpha),
    Rpg_(tensor::I),
    c1_(-1),
    dir1_(0)
{
    // Principal stresses - eigenvalues returned in ascending order
    vector lambda = eigenValues(R);

    // Eddy rotation from principal-to-global axes
    // - given by the 3 eigenvectors of the Reynold stress tensor as rows in
    //   the result tensor (transposed transformation tensor)
    // - returned in ascending eigenvalue order
    Rpg_ = eigenVectors(R, lambda).T();

    if (debug)
    {
        // Global->Principal transform = Rgp = Rpg.T()
        // Rgp & R & Rgp.T() should have eigenvalues on its diagonal and
        // zeros for all other components
        Pout<< "Rpg.T() & R & Rpg: " << (Rpg_.T() & R & Rpg_) << endl;
    }

    // Set the dfeddy orientation to position of max eigenvalue
    // (direction of dfeddy major axis, sigma_x in reference)
    dir1_ = 2;

    // Normalisation coefficient (eq. 11)
    // Note: sqrt(10*V)/sqrt(nEddy) applied outside when computing uDash
    c1_ = cmptAv(sigma_)/cmptProduct(sigma_)*cmptMin(sigma_);
}


Foam::dfeddy::dfeddy(const dfeddy& e)
:
    patchFaceI_(e.patchFaceI_),
    position0_(e.position0_),
    x_(e.x_),
    sigma_(e.sigma_),
    alpha_(e.alpha_),
    Rpg_(e.Rpg_),
    c1_(e.c1_),
    dir1_(e.dir1_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vectorField Foam::dfeddy::uDash(const pointField& xp, const vector& n) const
{
    vectorField uDash(xp.size(), vector::zero);

    forAll(uDash, label)
    {
        // Relative position inside dfeddy (global system)
        const vector r = cmptDivide(xp[label] - position(n), sigma_);

        if (mag(r) <= 1)
        {
            // Relative position inside dfeddy (dfeddy principal system)
            const vector rp = Rpg_.T() & r;

            // Shape function (dfeddy principal system)
            const vector q = cmptMultiply(sigma_, vector::one - cmptMultiply(rp, rp));

            // Fluctuating velocity (dfeddy principal system) (eq. 8)
            uDash[label] = cmptMultiply(q, rp^alpha_);
        }

    }

    // Convert into global system (eq. 10)
    return c1_*(Rpg_ & uDash);
}

// ************************************************************************* //
