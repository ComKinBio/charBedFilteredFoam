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
    Foam::CollidingCloud

Description
    Adds coolisions to kinematic clouds

SourceFiles
    CollidingCloudI.H
    CollidingCloud.C

\*---------------------------------------------------------------------------*/

#ifndef CollidingCloud_H
#define CollidingCloud_H

#include "particle.H"
#include "Cloud.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "fvMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class CloudType>
class CollisionModel;

class sourceFilter;
/*---------------------------------------------------------------------------*\
                       Class CollidingCloud Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class CollidingCloud
:
    public CloudType
{
public:

    // Public Typedefs

        //- Type of cloud this cloud was instantiated for
        typedef CloudType cloudType;

        //- Type of parcel the cloud was instantiated for
        typedef typename CloudType::particleType parcelType;

        //- Convenience typedef for this cloud type
        typedef CollidingCloud<CloudType> collidingCloudType;
        
        //test Tian 
        label DEMFlag;


private:

    // Private Data

        //- Cloud copy pointer
        autoPtr<CollidingCloud<CloudType>> cloudCopyPtr_;


protected:

    // Protected data

        //- Thermo parcel constant properties
        typename parcelType::constantProperties constProps_;


        // References to the cloud sub-models

            //- Collision model
            autoPtr<CollisionModel<CollidingCloud<CloudType>>>
                collisionModel_;


        // Initialisation

            //- Set cloud sub-models
            void setModels();


        // Cloud evolution functions
            
            //- resolve-collide particles conversion
            template<class TrackCloudType>
            void calcParticleCollide
            (
                TrackCloudType& cloud,
                typename parcelType::trackingData& td,
                const scalar deltaT
            );

            //- Move-collide particles
            template<class TrackCloudType>
            void moveCollide
            (
                TrackCloudType& cloud,
                typename parcelType::trackingData& td,
                const scalar deltaT
            );

            //- Reset state of cloud
            void cloudReset(CollidingCloud<CloudType>& c);


public:

    // Constructors

        //- Construct given carrier gas fields
        CollidingCloud
        (
            const word& cloudName,
            const volScalarField& rho,
            const volVectorField& U,
            const volScalarField& mu,
            const dimensionedVector& g,
            bool readFields = true
        );
        
        //- Construct given carrier gas fields
        CollidingCloud
        (
            const word& cloudName,
            const volScalarField& rho,
            const volVectorField& U,
            const volScalarField& mu,
            const dimensionedVector& g,
            const volVectorField& Uavg,
            const volScalarField& rhoavg,
            const volScalarField& muavg,
            const volScalarField& alphaavg,
            sourceFilter* sourceFilterModel,
            bool readFields = true
        );

        //- Copy constructor with new name
        CollidingCloud
        (
            CollidingCloud<CloudType>& c,
            const word& name
        );

        //- Copy constructor with new name - creates bare cloud
        CollidingCloud
        (
            const fvMesh& mesh,
            const word& name,
            const CollidingCloud<CloudType>& c
        );

        //- Disallow default bitwise copy construction
        CollidingCloud(const CollidingCloud&) = delete;

        //- Construct and return clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType>> clone(const word& name)
        {
            return autoPtr<Cloud<parcelType>>
            (
                new CollidingCloud(*this, name)
            );
        }

        //- Construct and return bare clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType>> cloneBare(const word& name) const
        {
            return autoPtr<Cloud<parcelType>>
            (
                new CollidingCloud(this->mesh(), name, *this)
            );
        }


    //- Destructor
    virtual ~CollidingCloud();


    // Member Functions

        // Access

            //- Return a reference to the cloud copy
            inline const CollidingCloud& cloudCopy() const;

            //- Return the constant properties
            inline const typename parcelType::constantProperties&
                constProps() const;


            // Sub-models

                //- Return const access to the collision model
                inline const CollisionModel<CollidingCloud<CloudType>>&
                    collision() const;

                //- Return reference to the collision model
                inline CollisionModel<CollidingCloud<CloudType>>&
                    collision();

        // Check

            //- Total rotational kinetic energy in the system
            inline scalar rotationalKineticEnergyOfSystem() const;


        // Cloud evolution functions

            //- Store the current cloud state
            void storeState();

            //- Reset the current cloud to the previously stored state
            void restoreState();

            //- Evolve the cloud
            void evolve();
            
            //- Particle conversion
            template<class TrackCloudType>
            void conversion
            (
                TrackCloudType& cloud,
                typename parcelType::trackingData& td
            );

            //- Particle motion
            template<class TrackCloudType>
            void motion
            (
                TrackCloudType& cloud,
                typename parcelType::trackingData& td
            );


        // I-O

            //- Print cloud information
            void info();

    //- Return the optional particle forces
            inline label DEMFlagSet(label flag);
            
    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const CollidingCloud&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CollidingCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CollidingCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
