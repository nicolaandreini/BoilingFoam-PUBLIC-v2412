//- Update mDot to account for microlayer models at the walls

if (mixture.found("microlayerModel"))
{
    const dictionary microDict(mixture.subDict("microlayerModel"));
    const Switch debugML(mixture.lookupOrDefault<Switch>("debug", false));

    //- Coefficients of the polynomial q'(deltaT) [W/m] fitted on the
    //  microlayer heat flux per unit contact-line length
    scalarField coeffs(microDict.lookup("coefficients"));

    if (coeffs.size() != 4)
    {
        Info<< "\nThe microlayer model requires 4 coefficients!\n"
            << abort(FatalError);
    }

    wordList activePatches(microDict.lookup("patches"));

    scalarField ap(volPointInterpolation::New(mesh).interpolate(alpha1));

    cutFaceIso faceCut(mesh, ap);

    forAll(activePatches, patchi)
    {
        label patchID(mesh.boundaryMesh().findPatchID(activePatches[patchi]));

        if (patchID == -1)
        {
            continue;
        }

        const fvPatchScalarField& TBf(T.boundaryFieldRef()[patchID]);
        const fvPatchScalarField& aBf(alpha1.boundaryFieldRef()[patchID]);
        scalarField& kappaBf(kappa.boundaryFieldRef()[patchID]);
        const labelUList& fCells(TBf.patch().faceCells());

        const label startFace(TBf.patch().start());

        forAll(aBf, facei)
        {
            label cellI = fCells[facei];

            mDot0[cellI] = 0;
            energySource[cellI] = 0;
            energyExp[cellI] = 0;

            //- Select interface cells only
            if (aBf[facei] < 1e-3 || aBf[facei] > 0.999)
            {
                continue;
            }

            if (faceCut.calcSubFace(startFace + facei, 0.5) != 0)
            {
                continue;
            }

            //- Compute the contact length
            point p0(faceCut.surfacePoints()[0]);
            point p1(faceCut.surfacePoints()[1]);

            dimensionedScalar lc("lc", dimLength, mag(p0 - p1));

            //- Set conduction to zero: the microlayer carries the heat flux
            kappaBf[facei] = 0.;

            scalar deltaTmp(TBf[facei] - Tsat.value());

            //- Re-compute mDot0 at the interface
            mDot0[cellI] =
            (
                lc.value()/(hf.value()*mesh.V()[cellI])
              * (
                    (coeffs[0]*deltaTmp)
                  + (coeffs[1]*pow(deltaTmp,2))
                  + (coeffs[2]*pow(deltaTmp,3))
                  + (coeffs[3]*pow(deltaTmp,4))
                )
            );

            energySource[cellI] = 0;

            energyExp[cellI] =
            (
                -lc.value()/mesh.V()[cellI]
              * (
                    (coeffs[0]*(deltaTmp))
                  + (coeffs[1]*pow(deltaTmp,2))
                  + (coeffs[2]*pow(deltaTmp,3))
                  + (coeffs[3]*pow(deltaTmp,4))
                )
            );

            if (debugML)
            {
                fvPatchScalarField& lcBf
                (
                    lcField->boundaryFieldRef()[patchID]
                );

                lcBf[facei] = lc.value();
            }
        }
    }
}
