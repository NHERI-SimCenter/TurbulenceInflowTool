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
    Foam::vorton

Description
    Class to describe eddies for the turbulentDFSEMInletFvPatchVectorField
    boundary condition.

SourceFiles
    vorton.C
    vortonI.H
    vortonIO.C

\*---------------------------------------------------------------------------*/

#ifndef vorton_H
#define vorton_H

#include "vector.H"
#include "point.H"
#include "tensor.H"
#include "Random.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class vorton;
class Istream;
class Ostream;

bool operator==(const vorton& a, const vorton& b);
bool operator!=(const vorton& a, const vorton& b);
Istream& operator>>(Istream& is, vorton& v);
Ostream& operator<<(Ostream& os, const vorton& v);


/*---------------------------------------------------------------------------*\
                            Class vorton Declaration
\*---------------------------------------------------------------------------*/

class vorton
{
    // Private data
        //- Vorton type
        word type_;

        //- Patch face index that spawned the vorton
        label patchFaceI_;

        //- Reference position
        point position0_;

        //- Distance from reference position in normal direction
        scalar x_;

        //- Length scales in 3-D space
        vector sigma_;

        //- Time-averaged intensity
        vector gamma_;

        //- Coordinate system transformation from local to global axes
        //  X-direction aligned with max stress eigenvalue
        tensor Rpg_;


    // Private Member Functions

        //- Return a number with zero mean and unit variance
        inline scalar epsi(Random& rndGen) const;

        tensor pinv(const tensor& A) const;


public:

    // Constructors

        //- Construct null
        vorton();

        //- Construct from Istream
        vorton(Istream& is);

        //- Construct from components
        vorton
        (
            const word type,
            const label patchFaceI,     // patch face index
            const point& position0,     // reference position
            const scalar x,             // distance from reference position
            const vector L0,            // length scale
            const symmTensor& R0,       // Stress tensor
            Random& rndGen
        );

        //- Construct from components
        vorton
        (
            const word type,
            const label patchFaceI,     // patch face index
            const point& position0,     // reference position
            const scalar x,             // distance from reference position
            const symmTensor& R0,       // Stress tensor
            const vector sigma,         // length scale
            const vector gamma
        );

        //- Construct copy
        vorton(const vorton& v);

        static int debug;


    // Public Member Functions

        // Access

            //- Return the vorton type
            inline word type() const;

            //- Return the patch face index that spawned the vorton
            inline label patchFaceI() const;

            //- Return the reference position
            inline const point& position0() const;

            //- Return the distance from the reference position
            inline scalar x() const;

            //- Return the Length scales in 3-D space
            inline const vector& sigma() const;

            //- Return the time-averaged intensity
            inline const vector& gamma() const;

            //- Return the coordinate system transformation from local
            //  principal to global axes
            inline const tensor& Rpg() const;

            //- Return the vorton position
            inline point position(const vector& n) const;


        // Helper functions

            //- Volume
            inline scalar volume() const;

            //- Move the vorton
            inline void move(const scalar dx);

            //- Eddy bounds
            inline boundBox bounds(const bool global = true) const;


        // Evaluate

            //- Return the fluctuating velocity contribution at local point xp
            vectorField uDash(const pointField& xp, const vector& n) const;


    // Member Operators

        void operator=(const vorton& v);


    // Friend Operators

        friend bool operator==(const vorton& a, const vorton& b)
        {
            return
                a.type_ == b.type_
             && a.patchFaceI_ == b.patchFaceI_
             && a.position0_ == b.position0_
             && a.x_ == b.x_
             && a.sigma_ == b.sigma_
             && a.gamma_ == b.gamma_
             && a.Rpg_ == b.Rpg_;
        }

        friend bool operator!=(const vorton& a, const vorton& b)
        {
            return !(a == b);
        }


    // IOstream Operators

        friend Istream& operator>>(Istream& is, vorton& v);
        friend Ostream& operator<<(Ostream& os, const vorton& v);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "vortonI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
