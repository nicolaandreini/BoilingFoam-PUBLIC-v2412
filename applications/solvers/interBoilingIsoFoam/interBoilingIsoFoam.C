/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd.
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

Application
    interBoilingIsoFoam

Group
    grpMultiphaseSolvers

Description
    Single-region transient solver for two incompressible, isothermally
    immiscible fluids with phase change (boiling/condensation), based on the
    Schrage interfacial mass-transfer model.

    The interface is captured geometrically with isoAdvector; no MULES
    limiter is used.

    Single-region counterpart of icoBoilingFoam: heat conduction in the solid
    is not solved, so a wall thermal boundary condition must be prescribed on
    the fluid patches.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "incompressibleTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"
#include "gravityMeshObject.H"
#include "processorFvPatchFields.H"
#include "isoAdvection.H"
#include "cutFaceIso.H"
#include "interfacePropertiesBoiling.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for two incompressible, immiscible fluids with"
        " phase change, using isoAdvector interface capturing."
    );

    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createUfIfPresent.H"

    turbulence->validate();

    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"
        #include "readSolverControls.H"

        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- PIMPLE loop
        while (pimple.loop())
        {
            const bool firstIter = pimple.firstIter();
            const bool finalIter = pimple.finalIter();

            if (firstIter || moveMeshOuterCorrectors)
            {
                // With phase change div(U) != 0, so the flux correction after
                // mesh motion must be driven by the divergence carried over
                // from the unrefined mesh rather than by zero.
                volScalarField absDivU
                (
                    "divU0",
                    fvc::div(fvc::absolute(phi, U))
                );

                interface.correct();

                mesh.update();

                if (mesh.changing())
                {
                    gh = (g & mesh.C()) - ghRef;
                    ghf = (g & mesh.Cf()) - ghRef;

                    MRF.update();

                    interface.surf().mapAlphaField();
                    alpha2 = 1.0 - alpha1;
                    alpha2.correctBoundaryConditions();
                    rho == alpha1*rho1 + alpha2*rho2;
                    rho.correctBoundaryConditions();
                    rho.oldTime() = rho;
                    alpha2.oldTime() = alpha2;

                    rhoCp == alpha1*rho1*cp1 + alpha2*rho2*cp2;
                    rhoCp.correctBoundaryConditions();
                    rhoCp.oldTime() = rhoCp;

                    if (correctPhi)
                    {
                        // Calculate absolute flux from the mapped surface
                        // velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);
                        mixture.correct();
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            if (frozenFlow)
            {
                #include "TEqn.H"
            }
            else
            {
                //- Correct interface before advection: this also updates the
                //  reconstruction scheme.
                interface.correct();
                mixture.correct();

                #include "alphaControls.H"
                #include "alphaEqnSubCycle.H"
                #include "hardBoundAlpha.H"

                #include "TEqn.H"

                #include "UEqn.H"

                // --- PISO loop
                while (pimple.correct())
                {
                    #include "pEqn.H"
                }

                if (pimple.turbCorr())
                {
                    turbulence->correct();
                }
            }

            if (computeFluxes)
            {
                #include "computeHeatFluxes.H"
            }
        }

        runTime.write();

        runTime.printExecutionTime(Info);

        if (stopWhenConverged && converged)
        {
            Info<< "All fluxes are converged."
                << " The simulation ends here." << endl;

            runTime.writeAndEnd();
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
