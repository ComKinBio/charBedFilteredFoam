Info<< "\nConstructing coalCloud "<< endl;
coalCloud coalParcels
(
    "coalCloud",
    rhoc,
    Uc,
    g,
    slgThermo,
    Ucavg,
    rhoavg,
    Tavg,
    muavg,
    kappaavg,
    Cpavg,
    alphacavg,
    O2avg,
    CO2avg,
    H2Oavg,
    H2avg,
    &sourceFilterModel
);

//default value for DEMFlag is 1, if DEM flag is not actived, the solver is normal
//may not need, just for sure, if run from modified cases

coalParcels.DEMFlagSet(1);
