/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    dynamicAlphaContactAngleKistlerFvPatchScalarField

Description
    Evaluation of the dynamic contact angle using the Hoffman function,
    based on Kistler (1993)

SourceFiles
    dynamicAlphaContactAngleKistlerFvPatchScalarField.C

Author
    Edin Berberovic 2008 - Updated by Alexander Rattner, 2012

\*---------------------------------------------------------------------------*/

#ifndef dynamicAlphaContactAngleKistlerFvPatchScalarField_H
#define dynamicAlphaContactAngleKistlerFvPatchScalarField_H

#include "alphaContactAngleTwoPhaseFvPatchScalarField.H"
#include "mathematicalConstants.H"
#include "surfaceFields.H"
#include "RiddersRoot.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class dynamicAlphaContactAngleFvPatch Declaration
\*---------------------------------------------------------------------------*/

class dynamicAlphaContactAngleKistlerFvPatchScalarField
:
    public alphaContactAngleTwoPhaseFvPatchScalarField
{
    // Private data

        //- Reference contact angle
        scalar theta0_;
        
        //- Limiting advancing contact angle
        scalar thetaA_;

        //- Limiting receeding contact angle
        scalar thetaR_;

        //eb Name of the dynamic viscosity field
        //word muName_;

        //eb Name of the surface tension
       // word sigmaName_;

public:

    //eb - Function object class, inverse of the Hoffman function
    class InverseHoffmanFunction : public FuncBase
    {

        // Private data

            //- Equilibrium contact angle (advancing or receding), as parameter
            scalar thetaE_;

    public:

        // Constructors

            //- Construct from data
            InverseHoffmanFunction(const scalar& thetaE)
            :
                thetaE_(thetaE)
            {}

        // operator()
            scalar operator()(scalar fHI) const
            {
                return
                  pow(fHI,0.706)
                -   pow((1+1.31*pow(fHI,0.99)),0.706)
                   *(1.0/pow(5.16,0.706))*atanh((1-cos(thetaE_))/2);
            }
    };



public:

    //- Runtime type information
    TypeName("dynamicAlphaContactAngleKistler");

    // Static data members

        //eb - Conversion factor for radians into degrees
        static const scalar convertToDeg;

        //eb - Conversion factor for degrees into radians
        static const scalar convertToRad;

        //eb - Contact angle corresponding to zero contact line velocity
        //static const scalar theta0;

    // Constructors

        //- Construct from patch and internal field
        dynamicAlphaContactAngleKistlerFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        dynamicAlphaContactAngleKistlerFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  dynamicAlphaContactAngleKistlerFvPatchScalarField
        //  onto a new patch
        dynamicAlphaContactAngleKistlerFvPatchScalarField
        (
            const dynamicAlphaContactAngleKistlerFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        dynamicAlphaContactAngleKistlerFvPatchScalarField
        (
            const dynamicAlphaContactAngleKistlerFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new dynamicAlphaContactAngleKistlerFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        dynamicAlphaContactAngleKistlerFvPatchScalarField
        (
            const dynamicAlphaContactAngleKistlerFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new dynamicAlphaContactAngleKistlerFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Evaluate and return dynamic contact-angle
        virtual tmp<scalarField> theta
        (
            const fvPatchVectorField& Up,
            const fvsPatchVectorField& nHat
        ) const;

        //eb - Evaluate the value of the Hoffman function
        scalar HoffmanFunction
        (
            const scalar& x
        ) const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
