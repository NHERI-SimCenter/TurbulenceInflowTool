/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  6
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      inflowProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        // mean velocity
        U_Profile
        {
                        referenceValue    (10 0 0);
                        profile           exponential;
                        direction         (0 0 1);
                        origin            (0 0 0);
                        referencePoint    (0 0 0.364);
                        alpha             (0.326 0 0);
        }


        // turbulent intensity (symmTensorField)
        I_Profile
        {
                        profile           exponential;
                        referenceValue    (0.208 0 0 0.182 0 0.152);
                        direction         (0 0 1);
                        origin            (0 0 0);                 
                        referencePoint    (0 0 0.364);
                        alpha             (-0.191 0 0 -0.123 0 -0.005);
        }

        // turbulence length scale profile for u component
        Lu_Profile
        {
                        profile           exponential;
                        referenceValue    (0.302 0.302 0.302);
                        direction         (0 0 1);
                        origin            (0 0 0);                 
                        referencePoint    (0 0 0.254);      
                        alpha             (0.473 0.473 0.473);         
        }

        // turbulence length scale profile for v component
        Lv_Profile
        {
                        profile           exponential;
                        referenceValue    (0.0815 0.0815 0.0815);
                        direction         (0 0 1);
                        origin            (0 0 0);                 
                        referencePoint    (0 0 0.254);
                        alpha             (0.881 0.881 0.881);
        }

        // turbulence length scale profile for w component
        Lw_Profile
        {
                        profile           exponential;
                        referenceValue    (0.0326 0.0326 0.0326);
                        direction         (0 0 1);
                        origin            (0 0 0);                 
                        referencePoint    (0 0 0.254);
                        alpha             (1.539 1.539 1.539);
        }


// ************************************************************************* //
