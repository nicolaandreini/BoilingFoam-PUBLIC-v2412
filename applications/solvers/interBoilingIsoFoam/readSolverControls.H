frozenFlow = pimple.dict().getOrDefault("frozenFlow", frozenFlow);

computeFluxes =
    pimple.dict().getOrDefault<Switch>("computeFluxes", computeFluxes);

stopWhenConverged =
    pimple.dict().getOrDefault<Switch>("stopWhenConverged", stopWhenConverged);

fluxTol = pimple.dict().getOrDefault<scalar>("fluxTolerance", fluxTol);

TTol = pimple.dict().getOrDefault<scalar>("TTolerance", TTol);
