/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

Class
    Foam::radiationModels::absorptionEmissionModels::cloud

Description
    Retrieves absorption/emission data from a cloud object, particle properties are averaged into space

SourceFiles
    cloudAverageAbsorptionEmission.C

\*---------------------------------------------------------------------------*/

#ifndef cloudAverageAbsorptionEmission_H
#define cloudAverageAbsorptionEmission_H

#include "absorptionEmissionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
namespace absorptionEmissionModels
{

/*---------------------------------------------------------------------------*\
                            Class cloudAvg Declaration
\*---------------------------------------------------------------------------*/

class cloudAvg
:
    public absorptionEmissionModel
{
    // Private Data

        //- Coefficients dictionary
        dictionary coeffsDict_;

        //- Cloud name(s)
        const wordList cloudNames_;


public:

    //- Runtime type information
    TypeName("cloudAvg");


    // Constructors

        //- Construct from components
        cloudAvg(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~cloudAvg();


    // Member Functions

        //- Absorption coefficient for dispersed phase
        virtual tmp<volScalarField> aDisp(const label bandI = 0) const;

        //- Emission coefficient for dispersed phase
        virtual tmp<volScalarField> eDisp(const label bandI = 0) const;

        //- Return emission contribution for dispersed phase
        virtual tmp<volScalarField> EDisp(const label bandI = 0) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace absorptionEmissionModels
} // End namespace radiationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
