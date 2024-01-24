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
    Foam::dfeddy

Description
    Class to describe eddies for the turbulentDFSEMInletFvPatchVectorField
    boundary condition.

SourceFiles
    dfeddy.C
    dfeddyI.H
    dfeddyIO.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentDFSEMInletFvPatchVectorField_dfeddy_H
#define turbulentDFSEMInletFvPatchVectorField_dfeddy_H

#include "vector.H"
#include "point.H"
#include "tensor.H"
#include "Random.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class dfeddy;
class Istream;
class Ostream;

bool operator==(const dfeddy& a, const dfeddy& b);
bool operator!=(const dfeddy& a, const dfeddy& b);
Istream& operator>>(Istream& is, dfeddy& e);
Ostream& operator<<(Ostream& os, const dfeddy& e);


/*---------------------------------------------------------------------------*\
                            Class dfeddy Declaration
\*---------------------------------------------------------------------------*/

class dfeddy
{
    // Private data

        static label Gamma2Values[8];
        static UList<label> Gamma2;

        //- Patch face index that spawned the dfeddy
        label patchFaceI_;

        //- Reference position
        point position0_;

        //- Distance from reference position in normal direction
        scalar x_;

        //- Length scales in 3-D space
        vector sigma_;

        //- Time-averaged intensity
        vector alpha_;

        //- Coordinate system transformation from local to global axes
        //  X-direction aligned with max stress eigenvalue
        tensor Rpg_;

        //- Model coefficient c1
        scalar c1_;

        //- Index of streamwise direction
        label dir1_;


    // Private Member Functions

        //- Set the dfeddy scales: length, intensity
        bool setScales
        (
            const scalar sigmaX,
            const label gamma2,
            const vector& e,
            const vector& lambda,
            vector& sigma,
            vector& alpha
        ) const;

        //- Return a number with zero mean and unit variance
        inline scalar epsi(Random& rndGen) const;


public:

    // Constructors

        //- Construct null
        dfeddy();

        //- Construct from Istream
        dfeddy(Istream& is);

        //- Construct from components
        dfeddy
        (
            const label patchFaceI,     // patch face index
            const point& position0,     // reference position
            const scalar x,             // distance from reference position
            const scalar sigmaX,        // length scale
            const symmTensor& R,        // Stress tensor
            Random& rndGen
        );

        //- Construct from components
        dfeddy
        (
            const label patchFaceI,     // patch face index
            const point& position0,     // reference position
            const scalar x,             // distance from reference position
            const symmTensor& R,        // Stress tensor
            const vector sigma,         // length scale
            const vector alpha
        );

        //- Construct copy
        dfeddy(const dfeddy& e);

        static int debug;


    // Public Member Functions

        // Access

            //- Return the patch face index that spawned the dfeddy
            inline label patchFaceI() const;

            //- Return the reference position
            inline const point& position0() const;

            //- Return the distance from the reference position
            inline scalar x() const;

            //- Return the lLength scales in 3-D space
            inline const vector& sigma() const;

            //- Return the time-averaged intensity
            inline const vector& alpha() const;

            //- Return the coordinate system transformation from local
            //  principal to global axes
            inline const tensor& Rpg() const;

            //- Return the model coefficient c1
            inline scalar c1() const;

            //- Return the dfeddy position
            inline point position(const vector& n) const;

            //- Return the index of the streamwise direction
            inline label dir1() const;

            //- Return random vector of -1 and 1's
            inline vector epsilon(Random& rndGen) const;


        // Helper functions

            //- Volume
            inline scalar volume() const;

            //- Move the dfeddy
            inline void move(const scalar dx);

            //- Eddy bounds
            inline boundBox bounds(const bool global = true) const;


        // Evaluate

            //- Return the fluctuating velocity contribution at local point xp
            vectorField uDash(const pointField& xp, const vector& n) const;


    // Member Operators

        void operator=(const dfeddy& e);


    // Friend Operators

        friend bool operator==(const dfeddy& a, const dfeddy& b)
        {
            return
                a.patchFaceI_ == b.patchFaceI_
             && a.position0_ == b.position0_
             && a.x_ == b.x_
             && a.sigma_ == b.sigma_
             && a.alpha_ == b.alpha_
             && a.Rpg_ == b.Rpg_
             && a.c1_ == b.c1_
             && a.dir1_ == b.dir1_;
        }

        friend bool operator!=(const dfeddy& a, const dfeddy& b)
        {
            return !(a == b);
        }


    // IOstream Operators

        friend Istream& operator>>(Istream& is, dfeddy& e);
        friend Ostream& operator<<(Ostream& os, const dfeddy& e);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dfeddyI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
