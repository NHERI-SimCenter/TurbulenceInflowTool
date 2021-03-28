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

#include "turbulentDFMInletFvPatchVectorField.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "IFstream.H"
#include "OFstream.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::turbulentDFMInletFvPatchVectorField::createFiles()
{
    if (nOutputFace_ > 0)
    {
        fileName baseDir = db().time().path();

        word outputPrefix = "postProcessing";
        word prefix = "turbulentDFMInlet";

        if (Pstream::parRun())
        {
            baseDir = baseDir/".."/outputPrefix;
        }
        else
        {
            baseDir = baseDir/outputPrefix;
        }

        const word startTimeName = db().time().timeName(db().time().startTime().value());

        label start = 0;

        for (label i=0; i<Pstream::myProcNo(); i++)
        {
            start = start + patchSize_[i];
        }

        forAll(outputFaceIndices_, i)
        {
            if (!filePtrs_.set(i))
            {
                fileName outputDir(baseDir/prefix/startTimeName);

                mkDir(outputDir);

                word fName(db().time().timeName(outputFaceIndices_[i]+start));

                // Check if file already exists
                IFstream is(outputDir/(fName + ".dat"));

                filePtrs_.set(i, new OFstream(outputDir/(fName + ".dat")));

                initStream(filePtrs_[i]);

                writeFileHeader(i);
            }
        }
    }
}

void Foam::turbulentDFMInletFvPatchVectorField::writeFileHeader(const label i)
{
    label start = 0;

    for (label i=0; i<Pstream::myProcNo(); i++)
    {
        start = start + patchSize_[i];
    }

    writeHeaderValue(filePtrs_[i], "face index", outputFaceIndices_[i]+start);
    writeHeaderValue(filePtrs_[i], "face centroid", patch().Cf()[outputFaceIndices_[i]]);
    writeHeaderValue(filePtrs_[i], "mean velocity magnitude", U_[outputFaceIndices_[i]]);
    writeHeaderValue(filePtrs_[i], "Reynolds stress tensor", R_[outputFaceIndices_[i]]);
    writeHeaderValue(filePtrs_[i], "integral length scales", L_[outputFaceIndices_[i]]);
    writeHeaderValue(filePtrs_[i], "synthetic method", "digitial filtering");
    writeHeaderValue(filePtrs_[i], "filter function", filterType_);

    writeCommented(filePtrs_[i], "Time");

    writeTabbed(filePtrs_[i], "ux");
    writeTabbed(filePtrs_[i], "uy");
    writeTabbed(filePtrs_[i], "uz");

    filePtrs_[i] << endl;
}

void Foam::turbulentDFMInletFvPatchVectorField::writeValues(const label i, const vector u)
{
    writeTime(filePtrs_[i]);

    filePtrs_[i]
        << tab << u.component(vector::X)
        << tab << u.component(vector::Y)
        << tab << u.component(vector::Z) << endl;
}


Foam::label Foam::turbulentDFMInletFvPatchVectorField::charWidth() const
{
    label addChars = 8;
    return IOstream::defaultPrecision() + addChars;
}

void Foam::turbulentDFMInletFvPatchVectorField::initStream(Ostream& os) const
{
    os.setf(ios_base::scientific, ios_base::floatfield);
    os.width(charWidth());
}

void Foam::turbulentDFMInletFvPatchVectorField::writeCommented
(
    Ostream& os,
    const string& str
) const
{
    os  << setw(1) << "#" << setw(1) << ' '
        << setf(ios_base::left) << setw(charWidth() - 2) << str.c_str();
}

void Foam::turbulentDFMInletFvPatchVectorField::writeTabbed
(
    Ostream& os,
    const string& str
) const
{
    os  << tab << setw(charWidth()) << str.c_str();
}

void Foam::turbulentDFMInletFvPatchVectorField::writeHeader
(
    Ostream& os,
    const string& str
) const
{
    os  << setw(1) << "#" << setw(1) << ' '
        << setf(ios_base::left) << setw(charWidth() - 2) << str.c_str() << nl;
}

template<class Type>
void Foam::turbulentDFMInletFvPatchVectorField::writeHeaderValue
(
    Ostream& os,
    const string& property,
    const Type& value
) const
{
    os  << setw(1) << '#' << setw(1) << ' '
        << setf(ios_base::left) << setw(charWidth() - 2) << property.c_str()
        << setw(1) << ':' << setw(1) << ' ' << value << nl;
}

void Foam::turbulentDFMInletFvPatchVectorField::writeTime(Ostream& os) const
{
    os  << setw(charWidth()) << db().time().timeName();
}

void Foam::turbulentDFMInletFvPatchVectorField::initialiseOutput()
{
    sort(outputFaceIndices_);

    if (!Pstream::parRun())
    {
        filePtrs_.setSize(nOutputFace_);
    }
    else
    {
        label start = 0;

        for (label i=0; i<Pstream::myProcNo(); i++)
        {
            start = start + patchSize_[i];
        }

        label currentNumb = 0;
        label currentStar = 0;

        forAll(outputFaceIndices_, i)
        {
            if (outputFaceIndices_[i]>=start&&outputFaceIndices_[i]<start+patchSize_[Pstream::myProcNo()])
            {
                if (currentNumb == 0)
                {
                    currentStar = i;
                }
                currentNumb += 1;
            }
        }

        labelList temp = outputFaceIndices_;

        if (currentNumb > 0)
        {
            nOutputFace_ = currentNumb;

            filePtrs_.setSize(nOutputFace_);
            outputFaceIndices_.setSize(nOutputFace_);

            forAll(outputFaceIndices_, i)
            {
                outputFaceIndices_[i] = temp[currentStar+i]-start;
            }
        }
        else
        {
            filePtrs_.setSize(0);
            outputFaceIndices_.setSize(0);
        }
    }

    createFiles();
}

void Foam::turbulentDFMInletFvPatchVectorField::initialisePatch()
{
    const vectorField nf(patch().nf());

    // Patch normal points into domain
    patchNormal_ = -gAverage(nf);

    // Check that patch is planar
    scalar error = max(magSqr(patchNormal_+nf));

    if (error > SMALL)
    {
        WarningInFunction << "Patch " << patch().name() << " is not planar" << endl;
    }

    forAll (patchNormal_, label)
    {
        if (mag(patchNormal_[label]) <= SMALL)
        {
            patchNormal_[label] = 0;
        }
    }
    // Set number of patch faces for each processor
    patchSize_[Pstream::myProcNo()] = patch().Cf().size();

    Pstream::gatherList(patchSize_);
    Pstream::scatterList(patchSize_);
}

void Foam::turbulentDFMInletFvPatchVectorField::initialiseParameters()
{   
    const vectorField Cf(patch().Cf());
    
    forAll(U_, label)
    {
        if (U_[label] < 0)
        {
            Pout << "error: the patch-normal velocity magnitude at the point " << Cf[label]
                 << " is no larger than 0, please modify the input parameters for U" << endl;
        }
        
        if (R_[label].component(symmTensor::XX) < 0)
        {
            Pout << "error: the Reynolds stress component R_XX at the point " << Cf[label] 
                 << " is no larger than 0, please modify the input parameters for R" << endl;
        }
        else
        {
            Lund_[label].component(tensor::XX) = sqrt(R_[label].component(symmTensor::XX));
            Lund_[label].component(tensor::YX) = R_[label].component(symmTensor::XY)/Lund_[label].component(tensor::XX);
            Lund_[label].component(tensor::ZX) = R_[label].component(symmTensor::XZ)/Lund_[label].component(tensor::XX);
            
            const scalar sqrLundYY = R_[label].component(symmTensor::YY)-sqr(Lund_[label].component(tensor::YX));
            
            if (sqrLundYY < 0)
            {
                Pout << "error: the Reynolds stress component R_YY at the point " << Cf[label]
                     << " is no larger than the square of Lund_YX,"
                     << " Please modify the input parameters for R" << endl;
            }
            else
            {
                Lund_[label].component(tensor::YY) = sqrt(sqrLundYY);
                Lund_[label].component(tensor::ZY) = (R_[label].component(symmTensor::YZ)-Lund_[label].component(tensor::YX)*Lund_[label].component(tensor::ZX))/Lund_[label].component(tensor::YY);
                
                const scalar sqrLundZZ = R_[label].component(symmTensor::ZZ)-sqr(Lund_[label].component(tensor::ZX))-sqr(Lund_[label].component(tensor::ZY));
                
                if (sqrLundZZ < 0)
                {
                    Pout << "error: the Reynolds stress component R_ZZ at the point " << Cf[label]
                         << " is no larger than sum of the squares of Lund_ZX and Lund_ZY,"
                         << " Please modify the input parameters for R" << endl;
                }
                else
                {
                    Lund_[label].component(tensor::ZZ) = sqrt(sqrLundZZ);
                    
                    const tensor sqrLund = cmptMultiply(Lund_[label], Lund_[label]);
                    const vector ex = sqrLund.x()/R_[label].xx();
                    const vector ey = sqrLund.y()/R_[label].yy();
                    const vector ez = sqrLund.z()/R_[label].zz();

                    const tensor E = tensor(ex, ey, ez);

                    L0_[label] = inv(E)&L_[label];
        
                    forAll(L0_[label], subLabel)
                    {
                        if (L0_[label][subLabel] < 0)
                        {
                            Pout << "error: the " << subLabel+1 << "-th component of the converted length scales at the point " << Cf[label] 
                                 << " is no larger than 0, please modify the input parameters for L" << endl;
                        }
                    }
                }
            }
        }
    }
}

void Foam::turbulentDFMInletFvPatchVectorField::initialiseVirtualGrid()
{
    // Set the origin and spacing of the virtual grid points
    boundBox patchBounds(patch().Cf());

    origin_ = patchBounds.min();

    delta_ = Foam::sqrt(gMin(patch().magSf()))/gridFactor_;
    
    scalar maxSpan = patchBounds.span().component(vector::Y);

    if (maxSpan < patchBounds.span().component(vector::Z))
    {
        maxSpan = patchBounds.span().component(vector::Z);
    }

    My_ = ceil(patchBounds.span().component(vector::Y)/delta_) + 1;
    Mz_ = ceil(patchBounds.span().component(vector::Z)/delta_) + 1;

    // Get the patch face centres from all processors
    vectorField Cf = gatherProc(patch().Cf());
    Cf -= origin_;

    // Set the indices of virtual grid points
    yindices_.setSize(Cf.size());
    zindices_.setSize(Cf.size());

    forAll(Cf, faceI)
    {
        scalar distY = Cf[faceI].component(vector::Y);
        scalar distZ = Cf[faceI].component(vector::Z);

        yindices_[faceI] = floor(distY/delta_);
        zindices_[faceI] = floor(distZ/delta_);

        scalar remY = distY - floor(distY/delta_) * delta_;
        scalar remZ = distZ - floor(distZ/delta_) * delta_;

        if (remY > delta_/2) 
        {
            yindices_[faceI] += 1;
        }

        if (remZ > delta_/2) 
        {
            zindices_[faceI] += 1;
        }        
    }

    // Set the number of grids for each processsor
    indicesPerProc_ = 0;
    rest_ = 0;

    if (Pstream::master())
    {
        indicesPerProc_ = floor(scalar(Cf.size())/(Pstream::nProcs()));
        rest_ = Cf.size()-(indicesPerProc_*(Pstream::nProcs()));

        Info << nl << "Generating Inflow for " << Cf.size() << " faces" << endl;
        Info << "Distributing Indices per Proc: " << indicesPerProc_ << endl;

        if (rest_ > 0)
        {
            Info << "First " << rest_ << " Procs will do +1 Indices" << nl << endl;
        }
    }

    Pstream::scatter(indicesPerProc_);
    Pstream::scatter(rest_);
}

void Foam::turbulentDFMInletFvPatchVectorField::initialiseFilterCoeff()
{
    // Gather the length scale field from each processors
    tensorField L(gatherProc(L0_));

    ny_.setSize(L.size());
    nz_.setSize(L.size());

    forAll(L, label)
    {
        ny_[label] = vector(ceil(L[label].xy()/delta_), ceil(L[label].yy()/delta_), ceil(L[label].zy()/delta_));
        nz_[label] = vector(ceil(L[label].xz()/delta_), ceil(L[label].yz()/delta_), ceil(L[label].zz()/delta_));
    }

    const labelVector nyMax = gMax(ny_);
    const labelVector nzMax = gMax(nz_);

    // Set size of virtual grid
    rndSize_.component(0) = (My_+2*nfK_*nyMax.component(0))*(Mz_+2*nfK_*nzMax.component(0));
    rndSize_.component(1) = (My_+2*nfK_*nyMax.component(1))*(Mz_+2*nfK_*nzMax.component(1));
    rndSize_.component(2) = (My_+2*nfK_*nyMax.component(2))*(Mz_+2*nfK_*nzMax.component(2));

    if (Pstream::master())
    {
        Info<< "Generating " << rndSize_.component(0)+rndSize_.component(1)+rndSize_.component(2) << " Random Numbers" << endl;
    }

    //For current processor, get start and end index in virtual list.
    label start;
    label size;

    if (Pstream::myProcNo() < rest_)
    {
        start = Pstream::myProcNo()*indicesPerProc_ + Pstream::myProcNo();
        size = indicesPerProc_+1;
    }
    else
    {
        start = Pstream::myProcNo()*indicesPerProc_ + rest_;
        size = indicesPerProc_;
    }

    //only fill element of current proc
    filterCoeffProcx.setSize(size);
    filterCoeffProcy.setSize(size);
    filterCoeffProcz.setSize(size);

    //loop through all indices in one proc, and set filter kernel size and fill with 0.0
    for (label subI = 0; subI < size; subI++)
    {
        //index in full array
        label I = subI+start;

        filterCoeffProcx[subI].setSize((2*nfK_*ny_[I].component(0)+1)*(2*nfK_*nz_[I].component(0)+1), 0.0);
        filterCoeffProcy[subI].setSize((2*nfK_*ny_[I].component(1)+1)*(2*nfK_*nz_[I].component(1)+1), 0.0);
        filterCoeffProcz[subI].setSize((2*nfK_*ny_[I].component(2)+1)*(2*nfK_*nz_[I].component(2)+1), 0.0);

        get2DFilterCoeff(filterCoeffProcx[subI], ny_[I].component(0), nz_[I].component(0));
        get2DFilterCoeff(filterCoeffProcy[subI], ny_[I].component(1), nz_[I].component(1));
        get2DFilterCoeff(filterCoeffProcz[subI], ny_[I].component(2), nz_[I].component(2));
    }
}

inline Foam::label 
Foam::turbulentDFMInletFvPatchVectorField::get1DIndex(label x, label y, label yDim)
{
    return x * yDim + y;
}


void Foam::turbulentDFMInletFvPatchVectorField::get1DFilterCoeff(scalarList& b, label n)
{
    const scalar pi = constant::mathematical::pi;

    if (n == 0)
    {
        b.setSize(1, 0.0);
    }
    else
    {
        scalar sum = 0.0;

        label N = nfK_*n;

        b.setSize(2*N+1);

        for (label k = 0; k < 2*N+1; k++)
        {
            if (filterType_ == "exponential")
            {
                b[k] = Foam::exp(-fabs(2.0*(k-N)/n));
            }
            else if (filterType_ == "gaussian")
            {
                b[k] = Foam::exp(-pi*Foam::sqr(scalar(k-N))/(2.0*Foam::sqr(scalar(n))));
            }
            else if (filterType_ == "bessel")
            {
                b[k] = 0.1505857*bessk0(fabs(1.0228626*(k-N)/n)+0.0178114);
            }
            else
            {
                Info << "filter coefficient function" << filterType_ << " does not exist (ERROR)" << endl;
            }

            sum += b[k]*b[k];
        }

        sum = Foam::sqrt(sum);
        b = b/sum;
    }
}

void Foam::turbulentDFMInletFvPatchVectorField::get2DFilterCoeff(scalarList& filter,label ny, label nz)
{
    scalarList by;
    scalarList bz;

    get1DFilterCoeff(by, ny);
    get1DFilterCoeff(bz, nz);

    for (label i = 0; i < 2*nfK_*ny+1; i++)
    {
        for (label j = 0; j < 2*nfK_*nz+1; j++)
        {
            filter[get1DIndex(i, j, 2*nfK_*nz+1)] = by[i]*bz[j];
        }
    }
}

Foam::scalarField
Foam::turbulentDFMInletFvPatchVectorField::getRandomField(label totalSize)
{
    List<scalarField> virtualRandomFieldProc;

    virtualRandomFieldProc.setSize(Pstream::nProcs());

    label size = floor(totalSize/Pstream::nProcs());
    label rest = totalSize-size*Pstream::nProcs();

    if (Pstream::myProcNo() < rest)
    {
        virtualRandomFieldProc[Pstream::myProcNo()].setSize(size+1);
    }
    else
    {
        virtualRandomFieldProc[Pstream::myProcNo()].setSize(size);
    }

    forAll (virtualRandomFieldProc[Pstream::myProcNo()], i)
    {
        virtualRandomFieldProc[Pstream::myProcNo()][i] = rndGen_.scalarNormal();
    }

    Pstream::gatherList(virtualRandomFieldProc);
    Pstream::scatterList(virtualRandomFieldProc);

    return ListListOps::combine<scalarField>(virtualRandomFieldProc, accessOp<scalarField>());
}

void Foam::turbulentDFMInletFvPatchVectorField::spatialCorr()
{
    Info<< "Generating spatial correlation" << endl;
    
    scalarField virtualRandomFieldx = getRandomField(rndSize_.component(0));
    scalarField virtualRandomFieldy = getRandomField(rndSize_.component(1));
    scalarField virtualRandomFieldz = getRandomField(rndSize_.component(2));

    const labelVector nyMax = gMax(ny_);
    const labelVector nzMax = gMax(nz_);

    if (periodicInY_)
    {
        scalarField virtualRandomFieldxp = virtualRandomFieldx;
        scalarField virtualRandomFieldyp = virtualRandomFieldy;
        scalarField virtualRandomFieldzp = virtualRandomFieldz;

        for (label i=0; i<nfK_*nyMax.component(0); i++)
        {
            label size_rnd = Mz_+2*nfK_*nzMax.component(0);

            label start1_rnd = (Mz_+2*nfK_*nzMax.component(0))*i;
            label start2_rnd = (Mz_+2*nfK_*nzMax.component(0))*(My_+i);

            SubField<scalar>(virtualRandomFieldxp, size_rnd, start1_rnd) = 
            SubField<scalar>(virtualRandomFieldxp, size_rnd, start2_rnd);

            start1_rnd = (Mz_+2*nfK_*nzMax.component(0))*(My_+2*nfK_*nyMax.component(0)-i-1);
            start2_rnd = (Mz_+2*nfK_*nzMax.component(0))*(2*nfK_*nyMax.component(0)-i-1);

            SubField<scalar>(virtualRandomFieldxp, size_rnd, start1_rnd) = 
            SubField<scalar>(virtualRandomFieldxp, size_rnd, start2_rnd);
        }

        for (label i=0; i<nfK_*nyMax.component(1); i++)
        {
            label size_rnd = Mz_+2*nfK_*nzMax.component(1);

            label start1_rnd = (Mz_+2*nfK_*nzMax.component(1))*i;
            label start2_rnd = (Mz_+2*nfK_*nzMax.component(1))*(My_+i);

            SubField<scalar>(virtualRandomFieldyp, size_rnd, start1_rnd) = 
            SubField<scalar>(virtualRandomFieldyp, size_rnd, start2_rnd);

            start1_rnd = (Mz_+2*nfK_*nzMax.component(1))*(My_+2*nfK_*nyMax.component(1)-i-1);
            start2_rnd = (Mz_+2*nfK_*nzMax.component(1))*(2*nfK_*nyMax.component(1)-i-1);

            SubField<scalar>(virtualRandomFieldyp, size_rnd, start1_rnd) = 
            SubField<scalar>(virtualRandomFieldyp, size_rnd, start2_rnd);
        }


        for (label i=0; i<nfK_*nyMax.component(2); i++)
        {
            label size_rnd = Mz_+2*nfK_*nzMax.component(2);

            label start1_rnd = (Mz_+2*nfK_*nzMax.component(2))*i;
            label start2_rnd = (Mz_+2*nfK_*nzMax.component(2))*(My_+i);

            SubField<scalar>(virtualRandomFieldzp, size_rnd, start1_rnd) = 
            SubField<scalar>(virtualRandomFieldzp, size_rnd, start2_rnd);

            start1_rnd = (Mz_+2*nfK_*nzMax.component(2))*(My_+2*nfK_*nyMax.component(2)-i-1);
            start2_rnd = (Mz_+2*nfK_*nzMax.component(2))*(2*nfK_*nyMax.component(2)-i-1);

            SubField<scalar>(virtualRandomFieldzp, size_rnd, start1_rnd) = 
            SubField<scalar>(virtualRandomFieldzp, size_rnd, start2_rnd);
        }

        virtualRandomFieldx = virtualRandomFieldxp;
        virtualRandomFieldy = virtualRandomFieldyp;
        virtualRandomFieldz = virtualRandomFieldzp;
    }

    if (periodicInZ_)
    {
        scalarField virtualRandomFieldxp = virtualRandomFieldx;
        scalarField virtualRandomFieldyp = virtualRandomFieldy;
        scalarField virtualRandomFieldzp = virtualRandomFieldz;

        for (label i=0; i<My_+2*nfK_*nyMax.component(0); i++)
        {
            label size_rnd = nfK_*nzMax.component(0);

            label start1_rnd = (Mz_+2*nfK_*nzMax.component(0))*i;
            label start2_rnd = (Mz_+2*nfK_*nzMax.component(0))*i+Mz_;

            SubField<scalar>(virtualRandomFieldxp, size_rnd, start1_rnd) = 
            SubField<scalar>(virtualRandomFieldxp, size_rnd, start2_rnd);

            start1_rnd = (Mz_+2*nfK_*nzMax.component(0))*i+Mz_+nfK_*nzMax.component(0);
            start2_rnd = (Mz_+2*nfK_*nzMax.component(0))*i+nfK_*nzMax.component(0);

            SubField<scalar>(virtualRandomFieldxp, size_rnd, start1_rnd) = 
            SubField<scalar>(virtualRandomFieldxp, size_rnd, start2_rnd);
        }

        for (label i=0; i<My_+2*nfK_*nyMax.component(1); i++)
        {
            label size_rnd = nfK_*nzMax.component(1);

            label start1_rnd = (Mz_+2*nfK_*nzMax.component(1))*i;
            label start2_rnd = (Mz_+2*nfK_*nzMax.component(1))*i+Mz_;

            SubField<scalar>(virtualRandomFieldyp, size_rnd, start1_rnd) = 
            SubField<scalar>(virtualRandomFieldyp, size_rnd, start2_rnd);

            start1_rnd = (Mz_+2*nfK_*nzMax.component(1))*i+Mz_+nfK_*nzMax.component(1);
            start2_rnd = (Mz_+2*nfK_*nzMax.component(1))*i+nfK_*nzMax.component(1);

            SubField<scalar>(virtualRandomFieldyp, size_rnd, start1_rnd) = 
            SubField<scalar>(virtualRandomFieldyp, size_rnd, start2_rnd);
        }

        for (label i=0; i<My_+2*nfK_*nyMax.component(2); i++)
        {
            label size_rnd = nfK_*nzMax.component(2);

            label start1_rnd = (Mz_+2*nfK_*nzMax.component(2))*i;
            label start2_rnd = (Mz_+2*nfK_*nzMax.component(2))*i+Mz_;

            SubField<scalar>(virtualRandomFieldzp, size_rnd, start1_rnd) = 
            SubField<scalar>(virtualRandomFieldzp, size_rnd, start2_rnd);

            start1_rnd = (Mz_+2*nfK_*nzMax.component(2))*i+Mz_+nfK_*nzMax.component(2);
            start2_rnd = (Mz_+2*nfK_*nzMax.component(2))*i+nfK_*nzMax.component(2);

            SubField<scalar>(virtualRandomFieldzp, size_rnd, start1_rnd) = 
            SubField<scalar>(virtualRandomFieldzp, size_rnd, start2_rnd);
        }

        virtualRandomFieldx = virtualRandomFieldxp;
        virtualRandomFieldy = virtualRandomFieldyp;
        virtualRandomFieldz = virtualRandomFieldzp;
    }

    label start, size;

    if (Pstream::myProcNo() < rest_)
    {
        start = Pstream::myProcNo()*indicesPerProc_ + Pstream::myProcNo();
        size = indicesPerProc_ + 1;
    }
    else
    {
        start = Pstream::myProcNo()*indicesPerProc_ + rest_;
        size = indicesPerProc_;
    }
    
    List<vectorField> virtualFilteredFieldProc;

    virtualFilteredFieldProc.setSize(Pstream::nProcs());
    virtualFilteredFieldProc[Pstream::myProcNo()].setSize(size, vector::zero);

    labelVector yOffset = nfK_*nyMax;
    labelVector zOffset = nfK_*nzMax;

    forAll (virtualFilteredFieldProc[Pstream::myProcNo()], subI)
    {
        label I = subI+start;

        label i = yindices_[I]; // i = yindices on virtual Grid
        label j = zindices_[I]; // j = zindices on virtual Grid

        vector u = vector::zero;

        for (label ii = 0; ii < 2*nfK_*ny_[I].component(0)+1; ii++)
        {
            label start_rnd = get1DIndex
            (
                i+yOffset.component(0)-nfK_*ny_[I].component(0)+ii,
                j+zOffset.component(0)-nfK_*nz_[I].component(0),
                Mz_+2*nfK_*nzMax.component(0)
            );

            label size_rnd = 2*nfK_*nz_[I].component(0)+1;
            label start_filt = get1DIndex(ii, 0, 2*nfK_*nz_[I].component(0)+1);

            scalarField rnd = SubField<scalar>(virtualRandomFieldx, size_rnd, start_rnd);
            scalarField filt = SubField<scalar>(filterCoeffProcx[subI], size_rnd, start_filt);

            u.component(0) += sumProd(rnd,filt);
        }

        for (label ii = 0; ii < 2*nfK_*ny_[I].component(1)+1; ii++)
        {
            label start_rnd = get1DIndex
            (
                i+yOffset.component(1)-nfK_*ny_[I].component(1)+ii,
                j+zOffset.component(1)-nfK_*nz_[I].component(1),
                Mz_+2*nfK_*nzMax.component(1)
            );

            label size_rnd = 2*nfK_*nz_[I].component(1)+1;
            label start_filt = get1DIndex(ii, 0, 2*nfK_*nz_[I].component(1)+1);

            scalarField rnd = SubField<scalar>(virtualRandomFieldy, size_rnd, start_rnd);
            scalarField filt = SubField<scalar>(filterCoeffProcy[subI], size_rnd, start_filt);

            u.component(1) += sumProd(rnd,filt);
        }

        for (label ii = 0; ii < 2*nfK_*ny_[I].component(2)+1; ii++)
        {
            label start_rnd = get1DIndex
            (
                i+yOffset.component(2)-nfK_*ny_[I].component(2)+ii,
                j+zOffset.component(2)-nfK_*nz_[I].component(2),
                Mz_+2*nfK_*nzMax.component(2)
            );

            label size_rnd = 2*nfK_*nz_[I].component(2)+1;
            label start_filt = get1DIndex(ii, 0, 2*nfK_*nz_[I].component(2)+1);

            scalarField rnd = SubField<scalar>(virtualRandomFieldz, size_rnd, start_rnd);
            scalarField filt = SubField<scalar>(filterCoeffProcz[subI], size_rnd, start_filt);

            u.component(2) += sumProd(rnd,filt);
        }

        virtualFilteredFieldProc[Pstream::myProcNo()][subI] = u;
    }

    Pstream::gatherList(virtualFilteredFieldProc);
    Pstream::scatterList(virtualFilteredFieldProc);

    vectorField virtualFilteredField = ListListOps::combine<vectorField>(virtualFilteredFieldProc, accessOp<vectorField>());

    start = 0;

    if (Pstream::myProcNo() > 0)
    {
        for (label i = 0; i < Pstream::myProcNo(); i++)
        {
            start = start + patchSize_[i];
        }
    }

    if (patchSize_[Pstream::myProcNo()] > 0)
    {
        uFluctFiltered_ = SubField<vector>(virtualFilteredField, patchSize_[Pstream::myProcNo()], start);
    }
    
    Info<< "Spatial correlation generated" << endl;
}

void Foam::turbulentDFMInletFvPatchVectorField::temporalCorr()
{
    Info<< "Generating temporal correlation" << endl;
    
    scalar dt = db().time().deltaT().value();

    forAll(uFluctTemporal_, faceI)
    {
        const vector L = vector(L0_[faceI].xx(),L0_[faceI].yx(),L0_[faceI].zx());
        const vector T = L/U_[faceI];

        for (label ii = 0; ii <= 2; ii++)
        {
            uFluctTemporal_[faceI].component(ii) = uFluctTemporalOld_[faceI].component(ii)*Foam::exp(-dt/T.component(ii))
                                                 + uFluctFiltered_[faceI].component(ii)*Foam::sqrt(1.0-Foam::exp(-2.0*dt/T.component(ii)));
        }
    }

    uFluctTemporalOld_ = uFluctTemporal_;
    
    Info<< "Temporal correlation generated" << endl;
}


Foam::scalar Foam::turbulentDFMInletFvPatchVectorField::bessi0(const scalar x)
{
   scalar ax, ans;
   scalar y;

   if ((ax=fabs(x)) < 3.75) 
   {
      y = x/3.75;
      y = y*y;
      ans = 1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
          + y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
   } 
   else 
   {
      y = 3.75/ax;
      ans = (exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
          + y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
          + y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
          + y*0.392377e-2))))))));
   }

   return ans;
}

Foam::scalar Foam::turbulentDFMInletFvPatchVectorField::bessk0(const scalar x)
{
   scalar y, ans;

   if (x <= 2.0) 
   {
      y=x*x/4.0;
      ans=(-log(x/2.0)*bessi0(x))+(-0.57721566+y*(0.42278420
         +y*(0.23069756+y*(0.3488590e-1+y*(0.262698e-2
         +y*(0.10750e-3+y*0.74e-5))))));
   } 
   else 
   {
      y = 2.0/x;
      ans = (exp(-x)/sqrt(x))*(1.25331414+y*(-0.7832358e-1
          + y*(0.2189568e-1+y*(-0.1062446e-1+y*(0.587872e-2
          + y*(-0.251540e-2+y*0.53208e-3))))));
   }

   return ans;
}

const Foam::pointToPointPlanarInterpolation&
Foam::turbulentDFMInletFvPatchVectorField::patchMapper() const
{
    // Initialise interpolation (2D planar interpolation by triangulation)
    if (mapperPtr_.empty())
    {
        // Reread values and interpolate
        fileName samplePointsFile
        (
            this->db().time().path()
           /this->db().time().caseConstant()
           /"boundaryData"
           /this->patch().name()
           /"points"
        );

        pointField samplePoints((IFstream(samplePointsFile)()));

        if (debug)
        {
            InfoInFunction
                << " Read " << samplePoints.size() << " sample points from "
                << samplePointsFile << endl;
        }

        // tbd: run-time selection
        bool nearestOnly =
        (
           !mapMethod_.empty()
         && mapMethod_ != "planarInterpolation"
        );

        // Allocate the interpolator
        mapperPtr_.reset
        (
            new pointToPointPlanarInterpolation
            (
                samplePoints,
                patch().Cf(),
                perturb_,
                nearestOnly
            )
        );
    }

    return *mapperPtr_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbulentDFMInletFvPatchVectorField::
turbulentDFMInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),

    curTimeIndex_(-1),
    patchSize_(Pstream::nProcs(), 0),
    patchNormal_(vector::zero),
    periodicInY_(false),
    periodicInZ_(false),

    uFluctFiltered_(p.size(), pTraits<vector>::zero),
    uFluctTemporalOld_(p.size(), pTraits<vector>::zero),
    uFluctTemporal_(p.size(), pTraits<vector>::zero),

    perturb_(1e-5),
    mapMethod_("planarInterpolation"),
    mapperPtr_(nullptr),

    interpolateU_(false),
    interpolateR_(false),
    interpolateL_(false),
    U_(p.size(), 0.0),
    R_(p.size(), pTraits<symmTensor>::zero),
    L_(p.size(), pTraits<tensor>::zero),
    L0_(p.size(), pTraits<tensor>::zero),
    Lund_(p.size(), pTraits<tensor>::zero),

    isInitialized_(false),
    isCleanRestart_(false),
    isRestart_(false),
    gridFactor_(1.0),
    origin_(vector::zero),
    My_(0),
    Mz_(0),
    delta_(0),
    ny_(),
    nz_(),
    nfK_(2),
    yindices_(),
    zindices_(),

    indicesPerProc_(0),
    rest_(0),

    rndGen_((Pstream::myProcNo()+1)*time(NULL)),
    filterType_("exponential"),
    rndSize_(vector::zero),
    filterCoeffProcx(),
    filterCoeffProcy(),
    filterCoeffProcz(),

    nOutputFace_(0),
    outputFaceIndices_(),
    filePtrs_()
{}


Foam::turbulentDFMInletFvPatchVectorField::
turbulentDFMInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),

    curTimeIndex_(-1),
    patchSize_(Pstream::nProcs(), 0),
    patchNormal_(vector::zero),
    periodicInY_(dict.lookupOrDefault<bool>("periodicInY", false)),
    periodicInZ_(dict.lookupOrDefault<bool>("periodicInZ", false)),

    uFluctFiltered_(p.size(), pTraits<vector>::zero),
    uFluctTemporalOld_(p.size(), pTraits<vector>::zero),
    uFluctTemporal_(p.size(), pTraits<vector>::zero),

    perturb_(dict.lookupOrDefault<scalar>("perturb", 1e-5)),
    mapMethod_(dict.lookupOrDefault<word>("mapMethod", "nearestCell")),
    mapperPtr_(nullptr),

    interpolateU_(false),
    interpolateR_(false),
    interpolateL_(false),
    U_(interpolateOrRead<scalar>("U", dict, interpolateU_)),
    R_(interpolateOrRead<symmTensor>("R", dict, interpolateR_)),
    L_(interpolateOrRead<tensor>("L", dict, interpolateL_)),
    L0_(p.size(), pTraits<tensor>::zero),
    Lund_(p.size(), pTraits<tensor>::zero),

    isInitialized_(false),
    isCleanRestart_(dict.lookupOrDefault<bool>("cleanRestart", false)),
    isRestart_(false),
    gridFactor_(dict.lookupOrDefault<scalar>("gridFactor", 1.0)),
    origin_(vector::zero),
    My_(0),
    Mz_(0),
    delta_(0),
    ny_(),
    nz_(),
    nfK_(dict.lookupOrDefault<label>("filterFactor", 2)),
    yindices_(),
    zindices_(),

    indicesPerProc_(0),
    rest_(0),

    rndGen_((Pstream::myProcNo()+1)*time(NULL)),
    filterType_(dict.lookupOrDefault<word>("filterType", "exponential")),
    rndSize_(vector::zero),
    filterCoeffProcx(),
    filterCoeffProcy(),
    filterCoeffProcz(),

    nOutputFace_(dict.lookupOrDefault<label>("nOutputFace", 0)),
    outputFaceIndices_(dict.lookupOrDefault<labelList>("outputFaceIndices", labelList(nOutputFace_, 0))),
    filePtrs_()
{
    if (dict.found("uFluctTemporal") && !isCleanRestart_)
    {
        isRestart_ = true;
        uFluctTemporal_ = vectorField("uFluctTemporal", dict, p.size());
    }
}


Foam::turbulentDFMInletFvPatchVectorField::
turbulentDFMInletFvPatchVectorField
(
    const turbulentDFMInletFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),

    curTimeIndex_(ptf.curTimeIndex_),
    patchSize_(ptf.patchSize_),
    patchNormal_(ptf.patchNormal_),
    periodicInY_(ptf.periodicInY_),
    periodicInZ_(ptf.periodicInZ_),

    uFluctFiltered_(mapper(ptf.uFluctFiltered_)),
    uFluctTemporalOld_(mapper(ptf.uFluctTemporalOld_)),
    uFluctTemporal_(mapper(ptf.uFluctTemporal_)),

    perturb_(ptf.perturb_),
    mapMethod_(ptf.mapMethod_),
    mapperPtr_(nullptr),

    interpolateU_(ptf.interpolateU_),
    interpolateR_(ptf.interpolateR_),
    interpolateL_(ptf.interpolateL_),
    U_(mapper(ptf.U_)),
    R_(mapper(ptf.R_)),
    L_(mapper(ptf.L_)),
    L0_(mapper(ptf.L0_)),
    Lund_(mapper(ptf.Lund_)),

    isInitialized_(ptf.isInitialized_),
    isCleanRestart_(ptf.isCleanRestart_),
    isRestart_(ptf.isRestart_),
    gridFactor_(ptf.gridFactor_),
    origin_(ptf.origin_),
    My_(ptf.My_),
    Mz_(ptf.Mz_),
    delta_(ptf.delta_),
    ny_(ptf.ny_),
    nz_(ptf.nz_),
    nfK_(ptf.nfK_),
    yindices_(ptf.yindices_),
    zindices_(ptf.zindices_),

    indicesPerProc_(ptf.indicesPerProc_),
    rest_(ptf.rest_),

    rndGen_(ptf.rndGen_),
    filterType_(ptf.filterType_),
    rndSize_(ptf.rndSize_),
    filterCoeffProcx(ptf.filterCoeffProcx),
    filterCoeffProcy(ptf.filterCoeffProcy),
    filterCoeffProcz(ptf.filterCoeffProcz),

    nOutputFace_(ptf.nOutputFace_),
    outputFaceIndices_(ptf.outputFaceIndices_),
    filePtrs_()
{}


Foam::turbulentDFMInletFvPatchVectorField::
turbulentDFMInletFvPatchVectorField
(
    const turbulentDFMInletFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),

    curTimeIndex_(ptf.curTimeIndex_),
    patchSize_(ptf.patchSize_),
    patchNormal_(ptf.patchNormal_),
    periodicInY_(ptf.periodicInY_),
    periodicInZ_(ptf.periodicInZ_),

    uFluctFiltered_(ptf.uFluctFiltered_),
    uFluctTemporalOld_(ptf.uFluctTemporalOld_),
    uFluctTemporal_(ptf.uFluctTemporal_),

    perturb_(ptf.perturb_),
    mapMethod_(ptf.mapMethod_),
    mapperPtr_(nullptr),

    interpolateU_(ptf.interpolateU_),
    interpolateR_(ptf.interpolateR_),
    interpolateL_(ptf.interpolateL_),
    U_(ptf.U_),
    R_(ptf.R_),
    L_(ptf.L_),
    L0_(ptf.L0_),
    Lund_(ptf.Lund_),

    isInitialized_(ptf.isInitialized_),
    isCleanRestart_(ptf.isCleanRestart_),
    isRestart_(ptf.isRestart_),
    gridFactor_(ptf.gridFactor_),
    origin_(ptf.origin_),
    My_(ptf.My_),
    Mz_(ptf.Mz_),
    delta_(ptf.delta_),
    ny_(ptf.ny_),
    nz_(ptf.nz_),
    nfK_(ptf.nfK_),
    yindices_(ptf.yindices_),
    zindices_(ptf.zindices_),

    indicesPerProc_(ptf.indicesPerProc_),
    rest_(ptf.rest_),

    rndGen_(ptf.rndGen_),
    filterType_(ptf.filterType_),
    rndSize_(ptf.rndSize_),
    filterCoeffProcx(ptf.filterCoeffProcx),
    filterCoeffProcy(ptf.filterCoeffProcy),
    filterCoeffProcz(ptf.filterCoeffProcz),

    nOutputFace_(ptf.nOutputFace_),
    outputFaceIndices_(ptf.outputFaceIndices_),
    filePtrs_()
{}


Foam::turbulentDFMInletFvPatchVectorField::
turbulentDFMInletFvPatchVectorField
(
    const turbulentDFMInletFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),

    curTimeIndex_(ptf.curTimeIndex_),
    patchSize_(ptf.patchSize_),
    patchNormal_(ptf.patchNormal_),
    periodicInY_(ptf.periodicInY_),
    periodicInZ_(ptf.periodicInZ_),

    uFluctFiltered_(ptf.uFluctFiltered_),
    uFluctTemporalOld_(ptf.uFluctTemporalOld_),
    uFluctTemporal_(ptf.uFluctTemporal_),

    perturb_(ptf.perturb_),
    mapMethod_(ptf.mapMethod_),
    mapperPtr_(nullptr),

    interpolateU_(ptf.interpolateU_),
    interpolateR_(ptf.interpolateR_),
    interpolateL_(ptf.interpolateL_),
    U_(ptf.U_),
    R_(ptf.R_),
    L_(ptf.L_),
    L0_(ptf.L0_),
    Lund_(ptf.Lund_),

    isInitialized_(ptf.isInitialized_),
    isCleanRestart_(ptf.isCleanRestart_),
    isRestart_(ptf.isRestart_),
    gridFactor_(ptf.gridFactor_),
    origin_(ptf.origin_),
    My_(ptf.My_),
    Mz_(ptf.Mz_),
    delta_(ptf.delta_),
    ny_(ptf.ny_),
    nz_(ptf.nz_),
    nfK_(ptf.nfK_),
    yindices_(ptf.yindices_),
    zindices_(ptf.zindices_),

    indicesPerProc_(ptf.indicesPerProc_),
    rest_(ptf.rest_),

    rndGen_(ptf.rndGen_),
    filterType_(ptf.filterType_),
    rndSize_(ptf.rndSize_),
    filterCoeffProcx(ptf.filterCoeffProcx),
    filterCoeffProcy(ptf.filterCoeffProcy),
    filterCoeffProcz(ptf.filterCoeffProcz),

    nOutputFace_(ptf.nOutputFace_),
    outputFaceIndices_(ptf.outputFaceIndices_),
    filePtrs_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::turbulentDFMInletFvPatchVectorField::initialise()
{
    initialisePatch();
    initialiseParameters();
    initialiseVirtualGrid();
    initialiseFilterCoeff();
    initialiseOutput();

    if (isRestart_ && !isCleanRestart_)
    {
        uFluctTemporalOld_ = uFluctTemporal_;
    }
    else
    {
        spatialCorr();
        uFluctTemporalOld_ = uFluctFiltered_;
    }
    
    Info<< "initialization complete" << endl;
}


void Foam::turbulentDFMInletFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    //initialize virtual grid points
    if(!isInitialized_)
    {
        initialise();
        isInitialized_ = true;
    }

    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        Info<< "updating velocity at Time = " << db().time().value() << endl;

        curTimeIndex_ = db().time().timeIndex();

        // Filter random field
        spatialCorr();

        // create new temporally correlated slice
        temporalCorr();

        // Set final velocity field
        vectorField& U = *this;
        U = (U_*patchNormal_) + (Lund_&uFluctTemporal_);

        // Re-scale to ensure correct flow rate
        scalar fCorr = gSum(U_*patch().magSf())/gSum(-U&patch().Sf());

        U *= fCorr;

        if (Pstream::master())
        {
            Info<< "mass flow correction coefficient: " << fCorr << endl;
        }

        forAll(outputFaceIndices_, i)
        {
            writeValues(i, U[outputFaceIndices_[i]]);
        }

        isRestart_ = true;
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::turbulentDFMInletFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "value", *this);

    writeEntryIfDifferent<bool>(os, "periodicInY", false, periodicInY_);
    writeEntryIfDifferent<bool>(os, "periodicInZ", false, periodicInZ_);

    writeEntry(os, "uFluctTemporal", uFluctTemporal_);

    writeEntry(os, "U", U_);
    writeEntry(os, "R", R_);
    writeEntry(os, "L", L_);

    writeEntryIfDifferent<scalar>(os, "gridFactor", 1.0, gridFactor_);
    writeEntryIfDifferent<label>(os, "filterFactor", 2, nfK_);
    writeEntryIfDifferent<word>(os, "filterType", "exponential", filterType_);

    if (nOutputFace_ > 0)
    {
        os.writeKeyword("nOutputFace") << nOutputFace_ << token::END_STATEMENT << nl;
        writeEntry(os, "outputFaceIndices", outputFaceIndices_);
    }
}

void Foam::turbulentDFMInletFvPatchVectorField::autoMap(const fvPatchFieldMapper& mapper)
{
    fixedValueFvPatchVectorField::autoMap(mapper);

    mapper(uFluctTemporal_, uFluctTemporal_);
    mapper(U_, U_);
    mapper(R_, R_);
    mapper(L_, L_);
}

void Foam::turbulentDFMInletFvPatchVectorField::rmap(const fvPatchField<vector>& ptf, const labelList& addr)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    const turbulentDFMInletFvPatchVectorField& tiptf = refCast<const turbulentDFMInletFvPatchVectorField>(ptf);

    uFluctTemporal_.rmap(tiptf.uFluctTemporal_, addr);

    U_.rmap(tiptf.U_, addr);
    R_.rmap(tiptf.R_, addr);
    L_.rmap(tiptf.L_, addr);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        turbulentDFMInletFvPatchVectorField
    );
}

// ************************************************************************* //
