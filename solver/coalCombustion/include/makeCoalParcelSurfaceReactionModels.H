/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#ifndef makeCoalParcelSurfaceReactionModels_H
#define makeCoalParcelSurfaceReactionModels_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "NoSurfaceReaction.H"
#include "COxidationDiffusionLimitedRate.H"
#include "COxidationIntrinsicRate.H"
#include "COxidationKineticDiffusionLimitedRate.H"
#include "COxidationKineticDiffusionLimitedRateWithAsh.H"
#include "COxidationKineticDiffusionLimitedRateWithAshSteam.H"
#include "COxidationHurtMitchell.H"
#include "COxidationMurphyShaddix.H"
#include "COxidationKineticDiffusionLimitedRateWithAshSteamMurphyShaddix.H"
#include "COxidationKineticDuffyMurphyShaddix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeCoalParcelSurfaceReactionModels(CloudType)                         \
                                                                               \
    makeSurfaceReactionModelType(COxidationDiffusionLimitedRate, CloudType);   \
    makeSurfaceReactionModelType                                               \
    (                                                                          \
        COxidationKineticDiffusionLimitedRate,                                 \
        CloudType                                                              \
    );                                                                         \
    makeSurfaceReactionModelType                                               \
    (                                                                          \
        COxidationKineticDiffusionLimitedRateWithAsh,                          \
        CloudType                                                              \
    );                                                                         \
    makeSurfaceReactionModelType                                               \
    (                                                                          \
        COxidationKineticDiffusionLimitedRateWithAshSteam,                     \
        CloudType                                                              \
    );                                                                         \
    makeSurfaceReactionModelType                                               \
    (                                                                          \
        COxidationKineticDiffusionLimitedRateWithAshSteamMurphyShaddix,        \
        CloudType                                                              \
    );                                                                         \
    makeSurfaceReactionModelType                                               \
    (                                                                          \
        COxidationKineticDuffyMurphyShaddix,                                   \
        CloudType                                                              \
    );                                                                         \
    makeSurfaceReactionModelType(COxidationIntrinsicRate, CloudType);          \
    makeSurfaceReactionModelType(COxidationHurtMitchell, CloudType);           \
    makeSurfaceReactionModelType(COxidationMurphyShaddix, CloudType);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
