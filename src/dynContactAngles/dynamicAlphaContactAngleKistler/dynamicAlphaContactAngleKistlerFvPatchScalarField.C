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

\*---------------------------------------------------------------------------*/

#include "dynamicAlphaContactAngleKistlerFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "fvPatchFields.H"
#include "volMesh.H"

#include "volFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * //

const scalar dynamicAlphaContactAngleKistlerFvPatchScalarField::convertToDeg =
    180.0/constant::mathematical::pi;

const scalar dynamicAlphaContactAngleKistlerFvPatchScalarField::convertToRad =
    constant::mathematical::pi/180.0;

//const scalar dynamicAlphaContactAngleKistlerFvPatchScalarField::theta0 = 90.0;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

dynamicAlphaContactAngleKistlerFvPatchScalarField::
dynamicAlphaContactAngleKistlerFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(p, iF),
    theta0_(0.0),
    thetaA_(0.0),
    thetaR_(0.0)
    //muName_("undefined"),
    //sigmaName_("undefined")
{}

dynamicAlphaContactAngleKistlerFvPatchScalarField::
dynamicAlphaContactAngleKistlerFvPatchScalarField
(
    const dynamicAlphaContactAngleKistlerFvPatchScalarField& acpsf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(acpsf, p, iF, mapper),
    theta0_(acpsf.theta0_),
    thetaA_(acpsf.thetaA_),
    thetaR_(acpsf.thetaR_)
    //muName_(acpsf.muName_),
    //sigmaName_(acpsf.sigmaName_)
{}


dynamicAlphaContactAngleKistlerFvPatchScalarField::
dynamicAlphaContactAngleKistlerFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(p, iF),
    theta0_(readScalar(dict.lookup("theta0"))),
    thetaA_(readScalar(dict.lookup("thetaA"))),
    thetaR_(readScalar(dict.lookup("thetaR")))
    //muName_(dict.lookup("muEffKistler")),
    //sigmaName_(dict.lookup("sigmaKistler"))
{
    evaluate();
}


dynamicAlphaContactAngleKistlerFvPatchScalarField::
dynamicAlphaContactAngleKistlerFvPatchScalarField
(
    const dynamicAlphaContactAngleKistlerFvPatchScalarField& acpsf
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(acpsf),
    theta0_(acpsf.theta0_),
    thetaA_(acpsf.thetaA_),
    thetaR_(acpsf.thetaR_)
    //muName_(acpsf.muName_),
    //sigmaName_(acpsf.sigmaName_)
{}


dynamicAlphaContactAngleKistlerFvPatchScalarField::
dynamicAlphaContactAngleKistlerFvPatchScalarField
(
    const dynamicAlphaContactAngleKistlerFvPatchScalarField& acpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(acpsf, iF),
    theta0_(acpsf.theta0_),
    thetaA_(acpsf.thetaA_),
    thetaR_(acpsf.thetaR_)
    //muName_(acpsf.muName_),
    //sigmaName_(acpsf.sigmaName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> dynamicAlphaContactAngleKistlerFvPatchScalarField::theta
(
    const fvPatchVectorField& Up,
    const fvsPatchVectorField& nHat
) const
{
    /*
    //eb - Lookup and return the patchField of dynamic viscosity of mixture
    //     and surface tension
    if((muName_ != "muEffKistler") || (sigmaName_ != "sigmaKistler"))
    {
        FatalErrorIn
        (
            "dynamicAlphaContactAngleKistlerFvPatchScalarField"
        )   << " muEffKistler or sigma set inconsitently, muEffKistler = "
            << muName_ << ", sigmaKistler = " << sigmaName_ << '.' << nl
            << "    Set both muEffKistler and sigmaKistler according to the "
            << "definition of dynamicAlphaContactAngleKistler"
            << "\n    on patch " << this->patch().name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalError);
    }
    */

    //RATTNER to change, we probably want to switch mu to muEff
    //const fvPatchField<scalar>& mup =
    //    patch().lookupPatchField<volScalarField, scalar>(muName_);


    //const fvPatchField<scalar>& sigmap =
    //    patch().lookupPatchField<volScalarField, scalar>(sigmaName_);

    const vectorField nf = patch().nf();

    // Calculate the component of the velocity parallel to the wall
    vectorField Uwall = Up.patchInternalField() - Up;
    Uwall -= (nf & Uwall)*nf;

    // Find the direction of the interface parallel to the wall
    vectorField nWall = nHat - (nf & nHat)*nf;

    // Normalise nWall
    nWall /= (mag(nWall) + SMALL);

    // Calculate Uwall resolved normal to the interface parallel to
    // the wall
    scalarField uwall = - nWall & Uwall;

    //eb - Calculate local Capillary number
    //scalarField Ca = mup*mag(uwall)/sigmap;
    
    
    
    const label patchi = this->patch().index();
    
    const volScalarField& nu1 =
        this->db().objectRegistry::lookupObject<volScalarField>("nu1");

    const dictionary& transportProperties =
    this->db().objectRegistry::lookupObject<IOdictionary>
    (
        "transportProperties"
    );

    word phase1Name (wordList(transportProperties.lookup("phases"))[0]);
    word phase2Name (wordList(transportProperties.lookup("phases"))[1]);

    const volScalarField& alpha =
        this->db().objectRegistry::lookupObject<volScalarField>
        (
            IOobject::groupName("alpha", phase1Name)
        );

    dimensionedScalar rho1(phase1Name + ":rho", transportProperties.subDict(phase1Name));
    dimensionedScalar sigmap("sigma", transportProperties);
    
    //dimensionedScalar rho1(transportProperties.subDict(phase1Name).lookup("rho"));
    //dimensionedScalar sigmap(transportProperties.lookup("sigma"));

    const fvPatchScalarField&  nu1p = nu1.boundaryField()[patchi];
    
    scalarField Ca(nu1p*rho1.value()*uwall/(sigmap.value() + SMALL));
    
    
    
    //eb - Instantiate function object InverseHoffmanFunction for thetaA and thetaR
    dynamicAlphaContactAngleKistlerFvPatchScalarField::InverseHoffmanFunction
    InvHoffFuncThetaA
    (
        convertToRad*thetaA_
    );

    dynamicAlphaContactAngleKistlerFvPatchScalarField::InverseHoffmanFunction
    InvHoffFuncThetaR
    (
        convertToRad*thetaR_
    );

    //eb - Calculate InverseHoffmanFunction for thetaA and thetaR using
    // RiddersRoot
    RiddersRoot RRInvHoffFuncThetaA(InvHoffFuncThetaA, 1.e-10);
    scalar InvHoffFuncThetaAroot = RRInvHoffFuncThetaA.root(0,65);

    RiddersRoot RRInvHoffFuncThetaR(InvHoffFuncThetaR, 1.e-10);
    scalar InvHoffFuncThetaRroot = RRInvHoffFuncThetaR.root(0,65);

    //eb - Calculate and return the value of contact angle on patch faces,
    //     a general approach: the product of Uwall and nWall is negative
    //     for advancing and positiv for receding motion.
    //     thetaDp is initialized to 90 degrees corresponding to no wall
    //     adhesion
    
    scalarField thetaDp(patch().size(), convertToRad*theta0_);
    forAll(uwall, pfacei)
    {
        if(uwall[pfacei] > 0.0)
        {
            thetaDp[pfacei] = HoffmanFunction(   mag(Ca[pfacei])
                                               + InvHoffFuncThetaAroot);
        }
        else if (uwall[pfacei] < 0.0)
        {
            thetaDp[pfacei] = HoffmanFunction(   mag(Ca[pfacei])
                                               + InvHoffFuncThetaRroot);
        }
    }
    return convertToDeg*thetaDp;
}


scalar dynamicAlphaContactAngleKistlerFvPatchScalarField::HoffmanFunction
(
    const scalar& x
) const
{
    return acos(1 - 2*tanh(pow(5.16*(x/(1+1.31*pow(x,0.99))),0.706)));
}


void dynamicAlphaContactAngleKistlerFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("theta0") << theta0_ << token::END_STATEMENT << nl; 
    os.writeKeyword("thetaA") << thetaA_ << token::END_STATEMENT << nl;
    os.writeKeyword("thetaR") << thetaR_ << token::END_STATEMENT << nl;
    //os.writeKeyword("muEffKistler") << muName_ << token::END_STATEMENT << nl;
    //os.writeKeyword("sigmaKistler") << sigmaName_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    dynamicAlphaContactAngleKistlerFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
