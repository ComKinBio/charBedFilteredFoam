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
    Foam::ThermoParcel

Description
    Thermodynamic parcel class with one/two-way coupling with the continuous
    phase. Includes Kinematic parcel sub-models, plus:
    - heat transfer

SourceFiles
    ThermoParcelI.H
    ThermoParcel.C
    ThermoParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef ThermoParcel_H
#define ThermoParcel_H

#include "particle.H"
#include "SLGThermo.H"
#include "demandDrivenEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ParcelType>
class ThermoParcel;

template<class ParcelType>
Ostream& operator<<
(
    Ostream&,
    const ThermoParcel<ParcelType>&
);

/*---------------------------------------------------------------------------*\
                       Class ThermoParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class ThermoParcel
:
    public ParcelType
{
    // Private Data

        //- Size in bytes of the fields
        static const std::size_t sizeofFields_;


public:

    //- Class to hold thermo particle constant properties
    class constantProperties
    :
        public ParcelType::constantProperties
    {

        // Private Data

            //- Particle initial temperature [K]
            demandDrivenEntry<scalar> T0_;

            //- Minimum temperature [K]
            demandDrivenEntry<scalar> TMin_;

            //- Maximum temperature [K]
            demandDrivenEntry<scalar> TMax_;

            //- Particle specific heat capacity [J/kg/K]
            demandDrivenEntry<scalar> Cp0_;
            
            //- Particle specific heat conductivity [W/m/K]
            demandDrivenEntry<scalar> kp0_;

            //- Particle emissivity [] (radiation)
            demandDrivenEntry<scalar> epsilon0_;

            //- Particle scattering factor [] (radiation)
            demandDrivenEntry<scalar> f0_;
            
            //- use heat conduction between particles flag
            demandDrivenEntry<bool> ppConductionFlag_;
            
            //- enable conductive radiaon
            demandDrivenEntry<bool> ppRadiationFlag_;
            
            //- enable convection correction
            demandDrivenEntry<bool> convectionCorrection_;


    public:

        // Constructors

            //- Null constructor
            constantProperties();

            //- Copy constructor
            constantProperties(const constantProperties& cp);

            //- Construct from dictionary
            constantProperties(const dictionary& parentDict);


        // Member Functions

            // Access

                //- Return const access to the particle initial temperature [K]
                inline scalar T0() const;

                //- Return const access to minimum temperature [K]
                inline scalar TMin() const;

                //- Return const access to maximum temperature [K]
                inline scalar TMax() const;

                //- Set the maximum temperature [K]
                inline void setTMax(const scalar TMax);

                //- Return const access to the particle specific heat capacity
                //  [J/kg/K]
                inline scalar Cp0() const;
                
                //- Return const access to the particle heat conductivity
                //  [W/m/K]
                inline scalar kp0() const;

                //- Return const access to the particle emissivity []
                //  Active for radiation only
                inline scalar epsilon0() const;

                //- Return const access to the particle scattering factor []
                //  Active for radiation only
                inline scalar f0() const;
                
                //- Return const access to heat conduction between particles flag
                inline bool ppConductionFlag() const;
                
                //- enable conductive radiaon
                inline bool ppRadiationFlag() const;
                
                //- enable convection correction
                inline bool convectionCorrection() const;
        };


    class trackingData
    :
        public ParcelType::trackingData
    {
    private:

        // Private Data

            //- Local copy of carrier specific heat field
            //  Cp not stored on carrier thermo, but returned as tmp<...>
            const volScalarField Cp_;

            //- Local copy of carrier thermal conductivity field
            //  kappa not stored on carrier thermo, but returned as tmp<...>
            const volScalarField kappa_;


            // Interpolators for continuous phase fields

                //- Temperature field interpolator
                autoPtr<interpolation<scalar>> TInterp_;

                //- Specific heat capacity field interpolator
                autoPtr<interpolation<scalar>> CpInterp_;

                //- Thermal conductivity field interpolator
                autoPtr<interpolation<scalar>> kappaInterp_;
                
                //- Averave Temperature field interpolator
                autoPtr<interpolation<scalar>> TavgInterp_;

                //- Averave Specific heat capacity field interpolator
                autoPtr<interpolation<scalar>> CpavgInterp_;

                //- Averave Thermal conductivity field interpolator
                autoPtr<interpolation<scalar>> kappaavgInterp_;

                //- Radiation field interpolator
                autoPtr<interpolation<scalar>> GInterp_;


            // Cached continuous phase properties

                //- Temperature [K]
                scalar Tc_;

                //- Specific heat capacity [J/kg/K]
                scalar Cpc_;
                
                //- Specific thermal conductivity [W/m/K]
                scalar kappac_;


    public:

        typedef typename ParcelType::trackingData::trackPart trackPart;

        // Constructors

            //- Construct from components
            template <class TrackCloudType>
            inline trackingData
            (
                const TrackCloudType& cloud,
                trackPart part = ParcelType::trackingData::tpLinearTrack
            );


        // Member Functions

            //- Return access to the locally stored carrier Cp field
            inline const volScalarField& Cp() const;

            //- Return access to the locally stored carrier kappa field
            inline const volScalarField& kappa() const;

            //- Return const access to the interpolator for continuous
            //  phase temperature field
            inline const interpolation<scalar>& TInterp() const;

            //- Return const access to the interpolator for continuous
            //  phase specific heat capacity field
            inline const interpolation<scalar>& CpInterp() const;

            //- Return const access to the interpolator for continuous
            //  phase thermal conductivity field
            inline const interpolation<scalar>& kappaInterp() const;
            
            //- Return const access to the interpolator for continuous
            //  phase average temperature field
            inline const interpolation<scalar>& TavgInterp() const;

            //- Return const access to the interpolator for continuous
            //  phase average specific heat capacity field
            inline const interpolation<scalar>& CpavgInterp() const;

            //- Return const access to the interpolator for continuous
            //  phase average thermal conductivity field
            inline const interpolation<scalar>& kappaavgInterp() const;

            //- Return const access to the interpolator for continuous
            //  radiation field
            inline const interpolation<scalar>& GInterp() const;

            //- Return the continuous phase temperature
            inline scalar Tc() const;

            //- Access the continuous phase temperature
            inline scalar& Tc();

            //- Return the continuous phase specific heat capacity
            inline scalar Cpc() const;

            //- Access the continuous phase specific heat capacity
            inline scalar& Cpc();
            
            //- Return the continuous phase specific heat conductivity
            inline scalar kappac() const;

            //- Access the continuous phase specific heat conductivity
            inline scalar& kappac();
    };


protected:

    // Protected data

        // Parcel properties

            //- Temperature [K]
            scalar T_;

            //- Specific heat capacity [J/kg/K]
            scalar Cp_;
            
            //- Heat conductivity [W/m/K]
            scalar kp_;
            
            //used to calculate particle to particle heat transfer dT_ by conduction, rest for radiation
            //- At the end of the claculation updte T_ by T_ + dT_*deltT in ReactingMultiphaseParcel.C [K/s] 
            scalar dT_;
            
            //particle radiation info, T4 sum
            scalar ppT4sum_;
            
            //particle radiation info, neighbour particle number
            scalar neighbourNum_;
            
            //particle radiation info, neighbour max temperature
            scalar neighbourMax_;
            
            //particle radiation info, neighbour min temperature
            scalar neighbourMin_;
        
            //- diamter before DEM update [m]
            scalar dBeforeDEM_;


    // Protected Member Functions

        //- Calculate new particle temperature
        template<class TrackCloudType>
        scalar calcHeatTransfer
        (
            TrackCloudType& cloud,
            trackingData& td,
            const scalar dt,           // timestep
            const scalar Re,           // Reynolds number
            const scalar Pr,           // Prandtl number - surface
            const scalar kappa,        // Thermal conductivity - surface
            const scalar NCpW,         // Sum of N*Cp*W of emission species
            const scalar Sh,           // explicit particle enthalpy source
            scalar& dhsTrans,          // sensible enthalpy transfer to carrier
            scalar& Sph                // linearised heat transfer coefficient
        );


public:
    
                
    // Static data members

        //- Runtime type information
        TypeName("ThermoParcel");

        //- String representation of properties
        AddToPropertyList
        (
            ParcelType,
            " T"
          + " Cp"
          + " kp"
          + " dT"
          + " ppT4sum"
          + " neighbourNum"
          + " neighbourMax"
          + " neighbourMin"
          + " dBeforeDEM"
        );


    // Constructors

        //- Construct from mesh, coordinates and topology
        //  Other properties initialised as null
        inline ThermoParcel
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti
        );

        //- Construct from a position and a cell, searching for the rest of the
        //  required topology. Other properties are initialised as null.
        inline ThermoParcel
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli
        );

        //- Construct from components
        inline ThermoParcel
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti,
            const label typeId,
            const scalar nParticle0,
            const scalar d0,
            const scalar dTarget0,
            const vector& U0,
            const vector& f0,
            const vector& angularMomentum0,
            const vector& torque0,
            const constantProperties& constProps
        );

        //- Construct from Istream
        ThermoParcel
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true
        );

        //- Construct as a copy
        ThermoParcel(const ThermoParcel& p);

        //- Construct as a copy
        ThermoParcel(const ThermoParcel& p, const polyMesh& mesh);

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new ThermoParcel(*this));
        }

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone(const polyMesh& mesh) const
        {
            return autoPtr<particle>(new ThermoParcel(*this, mesh));
        }

        //- Factory class to read-construct particles used for
        //  parallel transfer
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<ThermoParcel<ParcelType>> operator()(Istream& is) const
            {
                return autoPtr<ThermoParcel<ParcelType>>
                (
                    new ThermoParcel<ParcelType>(mesh_, is, true)
                );
            }
        };


    // Member Functions

        // Access

            //- Return const access to temperature
            inline scalar T() const;

            //- Return const access to specific heat capacity
            inline scalar Cp() const;
            
            //- Return const access to heat conductivity
            inline scalar kp() const;

            //- Return the parcel sensible enthalpy
            inline scalar hs() const;
            
            //- Return const access to conducti heat
            inline scalar dT() const;
            
            //- Return const access to particle radiation info, T4 sum
            inline scalar ppT4sum() const;
            
            //- Return const access to particle radiation info, neighbour particle number
            inline scalar neighbourNum() const;
            
            //- Return const access to particle radiation info, neighbour max temperature
            inline scalar neighbourMax() const;
            
            //- Return const access to particle radiation info, neighbour min temperature
            inline scalar neighbourMin() const;
            
            //- Return const access to particle diamter before DEM update [m]
            inline scalar dBeforeDEM() const;

        // Edit

            //- Return access to temperature
            inline scalar& T();

            //- Return access to specific heat capacity
            inline scalar& Cp();
            
            //- Return access to specific heat conductivity
            inline scalar& kp();

            //- Return access to conduction heat
            inline scalar& dT();
            
            //- Return access to particle radiation info, T4 sum
            inline scalar& ppT4sum();
            
            //- Return access to particle radiation info, neighbour particle number
            inline scalar& neighbourNum();
            
            //- Return access to particle radiation info, neighbour max temperature
            inline scalar& neighbourMax();
            
            //- Return access to particle radiation info, neighbour min temperature
            inline scalar& neighbourMin();
            
            //- Return const access to particle diamter before DEM update [m]
            inline scalar& dBeforeDEM();
            
            //- Return const access to the particle heat conductivity
            //  [J/kg/K]
            inline scalar kpp(const scalar T) const;

            

        // Main calculation loop

            //- Set cell values
            template<class TrackCloudType>
            void setCellValues(TrackCloudType& cloud, trackingData& td);

            //- Correct cell values using latest transfer information
            template<class TrackCloudType>
            void cellValueSourceCorrection
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar dt
            );

            //- Calculate surface thermo properties
            template<class TrackCloudType>
            void calcSurfaceValues
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar T,
                scalar& Ts,
                scalar& rhos,
                scalar& mus,
                scalar& Pr,
                scalar& kappas
            ) const;

            //- Update parcel properties over the time interval due to heat transfer
            template<class TrackCloudType>
            void calc
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar dt
            );
            
            //- Update parcel properties over the time interval due to motion
            template<class TrackCloudType>
            void calcV
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar dt
            );


        // I-O

            //- Read
            template<class CloudType>
            static void readFields(CloudType& c);

            //- Write
            template<class CloudType>
            static void writeFields(const CloudType& c);


    // Ostream Operator

        friend Ostream& operator<< <ParcelType>
        (
            Ostream&,
            const ThermoParcel<ParcelType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ThermoParcelI.H"
#include "ThermoParcelTrackingDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ThermoParcel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
