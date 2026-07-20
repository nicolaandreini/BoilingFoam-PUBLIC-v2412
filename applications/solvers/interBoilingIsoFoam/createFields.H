Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

Info<< "Reading field mDot\n" << endl;
volScalarField mDot
(
    IOobject
    (
        "mDot",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField mDot0
(
    IOobject
    (
        "mDot0",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mDot
);


Info<< "Reading transportProperties\n" << endl;
incompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();


const uniformDimensionedVectorField& g = meshObjects::gravity::New(runTime);

uniformDimensionedScalarField hRef
(
    IOobject
    (
        "hRef",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    dimensionedScalar("hRef", dimLength, Zero)
);

const dimensionedScalar ghRef
(
    mag(g.value()) > SMALL
  ? g & (cmptMag(g.value())/mag(g.value()))*hRef
  : dimensionedScalar("ghRef", g.dimensions()*dimLength, 0)
);

volScalarField gh("gh", (g & mesh.C()) - ghRef);
surfaceScalarField ghf("ghf", (g & mesh.Cf()) - ghRef);


Info<< "Creating turbulence model\n" << endl;
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, mixture)
);


//- This also builds the surface reconstruction scheme
Boiling::interfaceProperties interface(alpha1, U, phi, mixture);


// * * * * * * * * * * * Phase and phase-change properties * * * * * * * * * //

const dictionary phase1(mixture.subDict(mixture.phase1Name()));
const dictionary phase2(mixture.subDict(mixture.phase2Name()));
const dictionary pcProperties(mixture.subDict("PhaseChangeProperties"));

const dimensionedScalar k1
(
    "kappa",
    dimensionSet(1,1,-3,-1,0,0,0),
    readScalar(phase1.lookup("kappa"))
);

const dimensionedScalar cp1
(
    "Cp",
    dimensionSet(0,2,-2,-1,0,0,0),
    readScalar(phase1.lookup("Cp"))
);

const dimensionedScalar k2
(
    "kappa",
    dimensionSet(1,1,-3,-1,0,0,0),
    readScalar(phase2.lookup("kappa"))
);

const dimensionedScalar cp2
(
    "Cp",
    dimensionSet(0,2,-2,-1,0,0,0),
    readScalar(phase2.lookup("Cp"))
);

const dimensionedScalar Tsat
(
    "Tsat",
    dimTemperature,
    readScalar(pcProperties.lookup("Tsat"))
);

const dimensionedScalar R
(
    "R",
    dimensionSet(0,2,-2,-1,0,0,0),
    readScalar(pcProperties.lookup("R"))
);

const dimensionedScalar sigmaEvap
(
    "sigmaEvap",
    dimless,
    readScalar(pcProperties.lookup("sigmaEvap"))
);

const dimensionedScalar DmDot
(
    "DmDot",
    dimArea,
    readScalar(pcProperties.lookup("DmDot"))
);

const dimensionedScalar hf
(
    "hf",
    dimensionSet(0,2,-2,0,0,0,0),
    readScalar(pcProperties.lookup("hf"))
);

const dimensionedScalar Prt
(
    "Prt",
    dimless,
    mixture.lookupOrDefault<scalar>("Prt", 0.85)
);


// * * * * * * * * * * * * * * * Derived fields  * * * * * * * * * * * * * * //

Info<< "Calculating field kappa\n" << endl;
volScalarField kappa
(
    IOobject
    (
        "kappa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    k1*alpha1 + k2*alpha2
);

Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

scalar maxT = gMax(T);

Info<< "Calculating field rho\n" << endl;
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    alpha1*rho1 + alpha2*rho2
);

volScalarField rhoCp
(
    IOobject
    (
        "rhoCp",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1*cp1 + alpha2*rho2*cp2,
    alpha1.boundaryField().types()
);
rhoCp.oldTime();

// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

// Energy flux
surfaceScalarField rhoCpPhi
(
    IOobject
    (
        "rhoPhiCp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rhoCp)*phi
);

volScalarField continuitySource
(
    IOobject
    (
        "continuitySource",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "continuitySource",
        dimensionSet(1,-3,-1,0,0,0,0),
        0.0
    )
);

volScalarField energySource
(
    IOobject
    (
        "energySource",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "energySource",
        dimensionSet(1,-1,-3,-1,0,0,0),
        0.0
    )
);

dimensionedScalar evapRate("evapRate", dimensionSet(1,0,-1,0,0,0,0), 0.0);
dimensionedScalar totMdot("totMdot", dimensionSet(1,0,-1,0,0,0,0), 0.0);


Info<< "Calculating field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = -1;
scalar pRefValue = 0.0;
setRefCell(p, p_rgh, pimple.dict(), pRefCell, pRefValue);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());


volScalarField rAU
(
    IOobject
    (
        "rAU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rAU", dimTime/dimDensity, 1)
);


#include "createMRF.H"
#include "createFvOptions.H"


//- Geometric interface advection. This solver has no MULES fall-back.
Boiling::isoAdvection advector(alpha1, phi, U, interface.surf());


//- Optional debug field for the microlayer contact-line length
autoPtr<volScalarField> lcField;

if (mixture.found("microlayerModel"))
{
    const Switch debugML(mixture.lookupOrDefault<Switch>("debug", false));

    if (debugML)
    {
        lcField.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "lc",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("lc", dimLength, 0)
            )
        );
    }
}

//- Optional wall heat flux field
autoPtr<volVectorField> qwField;

if (mixture.found("writeQwall"))
{
    qwField.reset
    (
        new volVectorField
        (
            IOobject
            (
                "qw",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            -kappa*fvc::grad(T)
        )
    );
}


// * * * * * * * * * * * * * * Solver controls * * * * * * * * * * * * * * * //

bool frozenFlow(pimple.dict().getOrDefault("frozenFlow", false));

Switch computeFluxes
(
    pimple.dict().getOrDefault<Switch>("computeFluxes", false)
);

Switch stopWhenConverged
(
    pimple.dict().getOrDefault<Switch>("stopWhenConverged", false)
);

scalar fluxTol(pimple.dict().getOrDefault<scalar>("fluxTolerance", 1e-8));
scalar TTol(pimple.dict().getOrDefault<scalar>("TTolerance", 1e-8));

bool converged(false);
