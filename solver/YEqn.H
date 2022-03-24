tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        alphaRhoPhic,
        mesh.divScheme("div(alphaRhoPhi,Yi_h)")
    )
);


{
    combustion->correct();
    Qdot = combustion->Qdot();
    volScalarField Yt(0.0*Y[0]);

    const volScalarField muEff(turbulence->muEff());

    forAll(Y, i)
    {
        if (i != inertIndex && composition.active(i))
        {
            volScalarField& Yi = Y[i];
            
            fvScalarMatrix YiEqn
            (
                fvm::ddt(alphacavg, rhoc, Yi)
              + mvConvection->fvmDiv(alphaRhoPhic, Yi)
              - fvm::laplacian
                (
                    fvc::interpolate(alphacavg)
                    *fvc::interpolate(muEff),
                    Yi
                )
              ==
                coalParcels.SYi(i, Yi)
              + combustion->R(Yi)
              + fvOptions(rhoc, Yi)
            );

            YiEqn.relax();

            fvOptions.constrain(YiEqn);

//             YiEqn.solve(mesh.solver("Yi"));
            YiEqn.solve();

            fvOptions.correct(Yi);

            Yi.max(0.0);
            Yt += Yi;
        }
    }

    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);

}
