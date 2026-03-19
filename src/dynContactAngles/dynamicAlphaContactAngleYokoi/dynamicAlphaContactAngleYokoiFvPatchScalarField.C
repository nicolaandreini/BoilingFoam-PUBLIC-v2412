/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    
    Dynamic contact angle based on Eq. (58) in Yokoi et al., Phys. Fluids 21, 072102 (2009); doi: 10.1063/1.3158468.
    Implemented by Mirco Magnini, University of Nottingham

\*---------------------------------------------------------------------------*/

#include "dynamicAlphaContactAngleYokoiFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volMesh.H"

#include "volFields.H"
#include "mathematicalConstants.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicAlphaContactAngleYokoiFvPatchScalarField::
dynamicAlphaContactAngleYokoiFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(p, iF),
    theta0_(0.0),
    thetaA_(0.0),
    thetaR_(0.0),
    ka_(0.0),
    kr_(0.0)
{}


Foam::dynamicAlphaContactAngleYokoiFvPatchScalarField::
dynamicAlphaContactAngleYokoiFvPatchScalarField
(
    const dynamicAlphaContactAngleYokoiFvPatchScalarField& gcpsf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(gcpsf, p, iF, mapper),
    theta0_(gcpsf.theta0_),
    thetaA_(gcpsf.thetaA_),
    thetaR_(gcpsf.thetaR_),
    ka_(gcpsf.ka_),
    kr_(gcpsf.kr_)
{}


Foam::dynamicAlphaContactAngleYokoiFvPatchScalarField::
dynamicAlphaContactAngleYokoiFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(p, iF, dict),
    theta0_(readScalar(dict.lookup("theta0"))),
    thetaA_(readScalar(dict.lookup("thetaA"))),
    thetaR_(readScalar(dict.lookup("thetaR"))),
    ka_(readScalar(dict.lookup("ka"))),
    kr_(readScalar(dict.lookup("kr")))
{
    evaluate();
}


Foam::dynamicAlphaContactAngleYokoiFvPatchScalarField::
dynamicAlphaContactAngleYokoiFvPatchScalarField
(
    const dynamicAlphaContactAngleYokoiFvPatchScalarField& gcpsf
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(gcpsf),
    theta0_(gcpsf.theta0_),
    thetaA_(gcpsf.thetaA_),
    thetaR_(gcpsf.thetaR_),
    ka_(gcpsf.ka_),
    kr_(gcpsf.kr_)
{}


Foam::dynamicAlphaContactAngleYokoiFvPatchScalarField::
dynamicAlphaContactAngleYokoiFvPatchScalarField
(
    const dynamicAlphaContactAngleYokoiFvPatchScalarField& gcpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(gcpsf, iF),
    theta0_(gcpsf.theta0_),
    thetaA_(gcpsf.thetaA_),
    thetaR_(gcpsf.thetaR_),
    ka_(gcpsf.ka_),
    kr_(gcpsf.kr_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::dynamicAlphaContactAngleYokoiFvPatchScalarField::theta
(
    const fvPatchVectorField& Up,
    const fvsPatchVectorField& nHat
) const
{

    const vectorField nf(patch().nf());

    // Calculated the component of the velocity parallel to the wall
    vectorField Uwall(Up.patchInternalField() - Up);
    Uwall -= (nf & Uwall)*nf;

    // Find the direction of the interface parallel to the wall
    vectorField nWall(nHat - (nf & nHat)*nf);

    // Normalise nWall
    nWall /= (mag(nWall) + SMALL);

    // Calculate Uwall resolved normal to the interface parallel to
    // the interface
    scalarField uwall(-nWall & Uwall);

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
  //  dimensionedScalar rho1(transportProperties.subDict(phase1Name).lookup("rho"));

    // dimensionedScalar sigmap(transportProperties.lookup("sigma"));

    const fvPatchScalarField&  nu1p = nu1.boundaryField()[patchi];

    scalarField Ca(nu1p*rho1.value()*uwall/sigmap.value());    
   
    return pos(uwall)*min(theta0_+pow(mag(Ca)/ka_,1.0/3.0),thetaA_)+(1.0-pos(uwall))*max(theta0_-pow(mag(Ca)/kr_,1.0/3.0),thetaR_);        
}


void Foam::dynamicAlphaContactAngleYokoiFvPatchScalarField::write(Ostream& os) const
{
    alphaContactAngleTwoPhaseFvPatchScalarField::write(os);
    os.writeKeyword("theta0") << theta0_ << token::END_STATEMENT << nl;
    os.writeKeyword("thetaA") << thetaA_ << token::END_STATEMENT << nl;
    os.writeKeyword("thetaR") << thetaR_ << token::END_STATEMENT << nl;
    os.writeKeyword("ka") << ka_ << token::END_STATEMENT << nl;
    os.writeKeyword("kr") << kr_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        dynamicAlphaContactAngleYokoiFvPatchScalarField
    );
}


// ************************************************************************* //
