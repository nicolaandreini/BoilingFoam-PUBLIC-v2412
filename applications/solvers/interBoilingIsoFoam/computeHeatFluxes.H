scalar sumConduction = 0.0;
scalar sumConvection = 0.0;
scalar sumTot = 0.0;

surfaceScalarField heatFluxD(-fvc::interpolate(kappa)*fvc::snGrad(T));

const surfaceScalarField::Boundary& patchHeatFluxD = heatFluxD.boundaryField();

if (qwField)
{
    qwField.ref() = -kappa*fvc::grad(T);
}

Info<< "\nWall heat fluxes [W]" << endl;

forAll(patchHeatFluxD, patchi)
{
    if (isA<processorFvPatchScalarField>(T.boundaryField()[patchi]))
    {
        continue;
    }

    scalar conduction = gSum
    (
        mesh.magSf().boundaryField()[patchi]
       *heatFluxD.boundaryField()[patchi]
    );

    scalar convection = gSum
    (
        rhoCp.boundaryField()[patchi]
       *T.boundaryField()[patchi]
       *phi.boundaryField()[patchi]
    );

    Info<< mesh.boundary()[patchi].name()
        << " "
        << conduction
        << " "
        << convection
        << " "
        << conduction + convection
        << endl;

    sumConduction += conduction;
    sumConvection += convection;
}

sumTot = sumConduction + sumConvection;

Info<< "sum "
    << sumConduction
    << " "
    << sumConvection
    << " "
    << sumTot
    << nl << endl;

//- Compute maximum temperature
scalar maxTnew(gMax(T));

scalar Tres(mag(maxT - maxTnew)/(maxT + SMALL));

maxT = maxTnew;

Info<< "maxT = " << maxT << "   residuals on T = " << Tres << nl << endl;

if (Pstream::master())
{
    std::ofstream outfile;

    outfile.open("heatFluxBalance.csv", std::ios_base::app);
    outfile << runTime.value() << ", " << sumTot << "\n";
    outfile.close();

    outfile.open("maxResT.csv", std::ios_base::app);
    outfile << runTime.value() << ", " << maxT << ", " << Tres << "\n";
    outfile.close();
}

converged = (mag(sumTot) < fluxTol && Tres < TTol);
