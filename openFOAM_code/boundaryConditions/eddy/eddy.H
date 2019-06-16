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

Class
    Foam::eddy

Description
    Class to describe eddies for the turbulentSEMInletFvPatchVectorField
    boundary condition.

SourceFiles
    eddy.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentSEMInletFvPatchVectorField_eddy_H
#define turbulentSEMInletFvPatchVectorField_eddy_H

#include "vector.H"
#include "point.H"
#include "tensor.H"
#include "Random.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class eddy;
bool operator==(const eddy& a, const eddy& b);
bool operator!=(const eddy& a, const eddy& b);
Istream& operator>>(Istream& is, eddy& e);
Ostream& operator<<(Ostream& os, const eddy& e);

/*---------------------------------------------------------------------------*\
                            Class eddy Declaration
\*---------------------------------------------------------------------------*/

class eddy
{
    // Private data

        //- velocity distribution shape
        word shape_;

        //- Patch face index that spawned the eddy
        label patchFaceI_;

        //- Reference position
        point position0_;

        //- Distance from reference position in normal direction
        scalar x_;

        //- Convecting velocity
        scalar convectVelocity_;

        //- Length scales for u component
        vector sigma_u_;

        //- Length scales for v component
        vector sigma_v_;

        //- Length scales for w component
        vector sigma_w_;

        //- Time-averaged intensity
        vector epsilon_;

        //- Tensor
        tensor Lund_; 

    // Private Member Functions

        //- Return a number with zero mean and unit variance
        inline scalar epsi(Random& rndGen) const;


public:

    // Constructors

        //- Construct null
        eddy();

        //- Construct from Istream
        eddy(Istream& is);

        //- Construct from components
        eddy
        (
            const word shape,             // velocity distribution shape
            const label patchFaceI,       // patch face index
            const point& position0,       // reference position
            const scalar x,               // distance from reference position
            const scalar convectVelocity, // convecting velocity
            const vector sigma_u,         // length scale for u component
            const vector sigma_v,         // length scale for v component
            const vector sigma_w,         // length scale for w component
            const symmTensor& R,          // Stress tensor
            Random& rndGen
        );

        //- Construct copy
        eddy(const eddy& e);


    // Public Member Functions

        // Access

            //- Return the velocity distribution shape type
            inline const word shape() const;

            //- Return the patch face index that spawned the eddy
            inline label patchFaceI() const;

            //- Return the reference position
            inline const point& position0() const;

            //- Return the distance from the reference position
            inline scalar x() const;

            //- Return the convecting velocity
            inline const scalar& convectVelocity() const;

            //- Return the Length scales for u component
            inline const vector& sigma_u() const;

            //- Return the Length scales for v component
            inline const vector& sigma_v() const;

            //- Return the Length scales for w component
            inline const vector& sigma_w() const;

            //- Return the time-averaged intensity
            inline const vector& epsilon() const;

            //- Return the Lund tensor
            inline const tensor& Lund() const;

            //- Return the eddy position
            inline point position(const vector& n) const;

            //- Return random vector of -1 and 1's
            inline vector epsilon(Random& rndGen) const;


        // Helper functions

            //- Volume
            inline scalar volume() const;

            //- Move the eddy
            inline void move(const scalar dt);

            //- Eddy bounds
            inline boundBox bounds() const;

        // Evaluate

            //- Return the fluctuating velocity contribution at local point xp
            vector uDash(const point& xp, const vector& n) const;

            //- Return the value of velocity distirbution function
            scalar shape(const scalar& x) const;


    // Member Operators

        void operator=(const eddy& e);


    // Friend Operators

        friend bool operator==(const eddy& a, const eddy& b)
        {
            return
                a.shape_ == b.shape_
             && a.patchFaceI_ == b.patchFaceI_
             && a.position0_ == b.position0_
             && a.x_ == b.x_
             && a.convectVelocity_ == b.convectVelocity_
             && a.sigma_u_ == b.sigma_u_
             && a.sigma_v_ == b.sigma_v_
             && a.sigma_w_ == b.sigma_w_
             && a.epsilon_ == b.epsilon_
             && a.Lund_ == b.Lund_;
        }

        friend bool operator!=(const eddy& a, const eddy& b)
        {
            return !(a == b);
        }


    // IOstream Operators

        friend Istream& operator>>(Istream& is, eddy& e);
        friend Ostream& operator<<(Ostream& os, const eddy& e);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "eddyI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
