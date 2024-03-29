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

#include "ReactingMultiphaseParcel.H"
#include "mathematicalConstants.H"


#define CP_GAS 1.18698e3

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
const Foam::label Foam::ReactingMultiphaseParcel<ParcelType>::GAS(0);

template<class ParcelType>
const Foam::label Foam::ReactingMultiphaseParcel<ParcelType>::LIQ(1);

template<class ParcelType>
const Foam::label Foam::ReactingMultiphaseParcel<ParcelType>::SLD(2);


// * * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackCloudType>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::CpEff
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar p,
    const scalar T,
    const label idG,
    const label idL,
    const label idS
) const
{
    return
        this->Y_[GAS]*cloud.composition().Cp(idG, YGas_, p, T)
      + this->Y_[LIQ]*cloud.composition().Cp(idL, YLiquid_, p, T)
      + this->Y_[SLD]*cloud.composition().Cp(idS, YSolid_, p, T);
}


template<class ParcelType>
template<class TrackCloudType>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::HsEff
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar p,
    const scalar T,
    const label idG,
    const label idL,
    const label idS
) const
{
    return
        this->Y_[GAS]*cloud.composition().Hs(idG, YGas_, p, T)
      + this->Y_[LIQ]*cloud.composition().Hs(idL, YLiquid_, p, T)
      + this->Y_[SLD]*cloud.composition().Hs(idS, YSolid_, p, T);
}


template<class ParcelType>
template<class TrackCloudType>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::LEff
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar p,
    const scalar T,
    const label idG,
    const label idL,
    const label idS
) const
{
    return
        this->Y_[GAS]*cloud.composition().L(idG, YGas_, p, T)
      + this->Y_[LIQ]*cloud.composition().L(idL, YLiquid_, p, T)
      + this->Y_[SLD]*cloud.composition().L(idS, YSolid_, p, T);
}


template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::updateMassFractions
(
    const scalar mass0,
    const scalarField& dMassGas,
    const scalarField& dMassLiquid,
    const scalarField& dMassSolid
)
{
    scalarField& YMix = this->Y_;

    scalar massGas =
        this->updateMassFraction(mass0*YMix[GAS], dMassGas, YGas_);
    scalar massLiquid =
        this->updateMassFraction(mass0*YMix[LIQ], dMassLiquid, YLiquid_);
    scalar massSolid =
        this->updateMassFraction(mass0*YMix[SLD], dMassSolid, YSolid_);

    scalar massNew = max(massGas + massLiquid + massSolid, rootVSmall);

    YMix[GAS] = massGas/massNew;
    YMix[LIQ] = massLiquid/massNew;
    YMix[SLD] = 1.0 - YMix[GAS] - YMix[LIQ];

    return massNew;
}

template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::kpp 
(
    const scalar Tp
)
{
    scalar kgas = (-0.000000000003493)*Foam::pow3(Tp) + 0.000000003319264*Foam::sqr(Tp) + 0.000060059499759*Tp + 0.008533051948052;
    
    const scalar k_char = 0.1;
    
    scalar k_ash =   1.03*(1.0-0.65) + 0.65*kgas;  /* Sven */
    
    //id 0 = char, 1 = ash
    scalar k = YSolid()[0]*k_char + YSolid()[1]*k_ash;
    
    return k;
}


// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackCloudType>
void Foam::ReactingMultiphaseParcel<ParcelType>::setCellValues
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    ParcelType::setCellValues(cloud, td);
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::ReactingMultiphaseParcel<ParcelType>::cellValueSourceCorrection
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    // Re-use correction from reacting parcel
    ParcelType::cellValueSourceCorrection(cloud, td, dt);
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::ReactingMultiphaseParcel<ParcelType>::calc
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    typedef typename TrackCloudType::reactingCloudType reactingCloudType;
    const CompositionModel<reactingCloudType>& composition =
        cloud.composition();


    // Define local properties at beginning of timestep
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const scalar np0 = this->nParticle_;
    const scalar d0 = this->d_;
    const vector& U0 = this->U_;
    const scalar T0 = this->T_;
    const scalar mass0 = this->mass();

    const scalar pc = td.pc();

    const scalarField& YMix = this->Y_;
    const label idG = composition.idGas();
    const label idL = composition.idLiquid();
    const label idS = composition.idSolid();
    
    const label idChar= composition.localId(idS,"C");
    const label idAsh= composition.localId(idS,"Ash");
    
    const scalar rhoChar = composition.solids().properties()[idChar].rho();
    const scalar rhoAsh = composition.solids().properties()[idAsh].rho();
    
    scalar mChar = mass0*composition.YMixture0()[idS]*composition.Y0(idS)[idChar];
    if(mChar > small)
    {
        canCombust_ = 1;
    }
    scalar dChar = cbrt(6.0*mChar/rhoChar/3.141593);
    
    const scalar e_bedMin = (1.0 - cloud.constProps().alphaMax());
    const scalar charSphericity = cloud.constProps().charSphericity();
//     const scalar ashSphericity = cloud.constProps().ashSphericity();
    const scalar ashInFuel = cloud.constProps().ashInFuel();
    const scalar shrinkageFactorAsh = cloud.constProps().shrinkageFactorAsh();
//     const bool usePpRadiation = cloud.constProps().ppRadiationFlag(); 
    const bool coarseGridSurfaceCombustion = cloud.constProps().coarseGridSurfaceCombustion();
    const bool useGasFilter = coarseGridSurfaceCombustion;
    label particleShape = 1;
    
    // 0-O2, 1-h2o, 2-co2, 3-h2, 4-e_bed
    scalarField Csavg(6, 0.0); 
    Csavg[0] = cloud.o2avg()[this->cell()];
//     Csavg[1] = cloud.h2oavg()[this->cell()];
    Csavg[2] = cloud.co2avg()[this->cell()];
//     Csavg[3] = cloud.h2avg()[this->cell()];
    
    // Calc surface values
    scalar Tcavg, Ts, rhos, mus, Prs, Cpcs, kappas, e_bed, Res;
    
    if (useGasFilter)
    {
        
        Tcavg = cloud.Tavg()[this->cell()];

        // Surface temperature using two thirds rule
        Ts = (2.0*T0 + Tcavg)/3.0;

        if (Ts < cloud.constProps().TMin())
        {
            if (debug)
            {
                WarningInFunction
                    << "Limiting parcel surface temperature to "
                    << cloud.constProps().TMin() <<  nl << endl;
            }

            Ts = cloud.constProps().TMin();
        }

        // Assuming thermo props vary linearly with T for small d(T)
        const scalar TRatio = Tcavg/Ts;

        rhos = cloud.rhoavg()[this->cell()]*TRatio;

        Cpcs = cloud.cpavg()[this->cell()];
        mus = cloud.muavg()[this->cell()];
        kappas = cloud.kappaavg()[this->cell()];

        Prs = Cpcs*mus/kappas;
        Prs = max(rootVSmall, Prs);
        
        Res = this->Re(rhos, U0, cloud.Uavg()[this->cell()], this->d_, mus);
        
        e_bed = max(cloud.alphaavg()[this->cell()], e_bedMin);
    }
    else
    {
        // Calc surface values
        this->calcSurfaceValues(cloud, td, T0, Ts, rhos, mus, Prs, kappas);
        
        Tcavg = td.Tc();
        
        // Reynolds number
        Res = this->Re(rhos, U0, td.Uc(), this->d_, mus);
        
        e_bed = max(1 - this->volume()/cloud.pMesh().cellVolumes()[this->cell()], e_bedMin);
    }
        
    Csavg[4] = e_bed;
    
    Csavg[5] = shrinkageFactorAsh; 


    // Sources
    //~~~~~~~~

//     // Explicit momentum source for particle
//     vector Su = Zero;

    // Linearised momentum source coefficient
//     scalar Spu = 0.0;

    // Momentum transfer from the particle to the carrier phase
//     vector dUTrans = Zero;

    // Explicit enthalpy source for particle
    scalar Sh = 0.0;

    // Linearised enthalpy source coefficient
    scalar Sph = 0.0;

    // Sensible enthalpy transfer from the particle to the carrier phase
    scalar dhsTrans = 0.0;


    // 1. Compute models that contribute to mass transfer - U, T held constant
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Phase change in liquid phase
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Mass transfer due to phase change
    scalarField dMassPC(YLiquid_.size(), 0.0);

    // Molar flux of species emitted from the particle (kmol/m^2/s)
    scalar Ne = 0.0;

    // Sum Ni*Cpi*Wi of emission species
    scalar NCpW = 0.0;

    // Surface concentrations of emitted species
//     scalarField Cs(composition.carrier().species().size(), 0.0);

    // Calc mass and enthalpy transfer due to phase change
//     this->calcPhaseChange
//     (
//         cloud,
//         td,
//         dt,
//         Res,
//         Prs,
//         Ts,
//         mus/rhos,
//         d0,
//         T0,
//         mass0,
//         idL,
//         YMix[LIQ],
//         YLiquid_,
//         dMassPC,
//         Sh,
//         Ne,
//         NCpW,
//         Cs
//     );


    // Devolatilisation
    // ~~~~~~~~~~~~~~~~

    // Mass transfer due to devolatilisation
//     scalarField dMassDV(YGas_.size(), 0.0);
//     scalarField dMassSOLID(YGas_.size(), 0.0);

    // Calc mass and enthalpy transfer due to devolatilisation
//     calcDevolatilisation
//     (
//         cloud,
//         td,
//         dt,
//         this->age_,
//         Ts,
//         d0,
//         T0,
//         mass0,
//         this->mass0_,
//         YMix[GAS]*YGas_,
//         YMix[LIQ]*YLiquid_,
//         YMix[SLD]*YSolid_,
//         canCombust_,
//         dMassDV,
//         dMassSOLID,//added
//         Sh,
//         Ne,
//         NCpW,
//         Cs
//     );


    // Surface reactions
    // ~~~~~~~~~~~~~~~~~

    // Change in carrier phase composition due to surface reactions
    scalarField dMassSRGas(YGas_.size(), 0.0);
    scalarField dMassSRLiquid(YLiquid_.size(), 0.0);
    scalarField dMassSRSolid(YSolid_.size(), 0.0);
    scalarField dMassSRCarrier(composition.carrier().species().size(), 0.0);
    
    scalar defaultScalar = 0;
    scalarField defaultScalarField(YGas_.size(), 0.0);
//     label defaultLabel = 0;

    // Calc mass and enthalpy transfer due to surface reactions
    calcSurfaceReactions
    (
        cloud,
        td,
        dt,
        d0,
        dChar,//di, inner layer radius
        T0,
        mass0,
        canCombust_,
        Ne,
        YMix,
        YGas_,
        YLiquid_,
        YSolid_,
        dMassSRGas,
        dMassSRLiquid,
        dMassSRSolid,
        dMassSRCarrier,
        Sh,
        dhsTrans,
        defaultScalar,//QComb, record combustion heat
        Res, //Re, Re number
        Tcavg, //Tc, continue phase temperature
        rhos, //rhoc, continue phase
        mus, //muc, continue phase viscosity
        charSphericity, //Xi, cylinder shape const, changed to charSphericity
        defaultScalarField, //dHeatSRCarrier, SR heat
        particleShape, //particleShape
        Csavg, //Csavg, carriers surface concentrations and bed voidage
        useGasFilter, //coarseGrid,
        ashInFuel, //deq, cylinder eq d, changed to ashInFuel
        false, //heatRatio 
        false //limitedCombustion
    );


    // 2. Update the parcel properties due to change in mass
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//     scalarField dMassGas(dMassDV + dMassSRGas);
    scalarField dMassLiquid(dMassPC + dMassSRLiquid);
    scalarField dMassSolid(dMassSRSolid);
    scalar mass1 =
        updateMassFractions(mass0, dMassSRGas, dMassLiquid, dMassSolid);

    this->Cp_ = CpEff(cloud, td, pc, T0, idG, idL, idS);

    // Update particle density or diameter
    scalar Vchar = mass1*composition.YMixture0()[idS]*composition.Y0(idS)[idChar]/rhoChar;
    scalar Vash = mass1*composition.YMixture0()[idS]*composition.Y0(idS)[idAsh]/rhoAsh/(1-shrinkageFactorAsh);
    this->d_ = cbrt(6.0*(Vchar + Vash)/3.141593);    
    this->rho_ = mass1/this->volume();
    
//     if (cloud.constProps().constantVolume())
//     {
//         this->rho_ = mass1/this->volume();
//     }
//     else
//     {
//         this->d_ = cbrt(mass1/this->rho_*6.0/pi);
//     }

    // Remove the particle when mass falls below minimum threshold
//     if (np0*mass1 < cloud.constProps().minParcelMass())
//     {
//         td.keepParticle = false;
// 
//         if (cloud.solution().coupled())
//         {
//             scalar dm = np0*mass0;
// 
//             // Absorb parcel into carrier phase
//             forAll(YGas_, i)
//             {
//                 label gid = composition.localToCarrierId(GAS, i);
//                 cloud.rhoTrans(gid)[this->cell()] += dm*YMix[GAS]*YGas_[i];
//             }
//             forAll(YLiquid_, i)
//             {
//                 label gid = composition.localToCarrierId(LIQ, i);
//                 cloud.rhoTrans(gid)[this->cell()] += dm*YMix[LIQ]*YLiquid_[i];
//             }
// 
//             // No mapping between solid components and carrier phase
//             /*
//             forAll(YSolid_, i)
//             {
//                 label gid = composition.localToCarrierId(SLD, i);
//                 cloud.rhoTrans(gid)[this->cell()] += dm*YMix[SLD]*YSolid_[i];
//             }
//             */
// 
//             cloud.UTrans()[this->cell()] += dm*U0;
// 
//             cloud.hsTrans()[this->cell()] +=
//                 dm*HsEff(cloud, td, pc, T0, idG, idL, idS);
// 
//             cloud.phaseChange().addToPhaseChangeMass(np0*mass1);
//         }
// 
//         return;
//     }

    // Correct surface values due to emitted species
//     this->correctSurfaceValues(cloud, td, Ts, Cs, rhos, mus, Prs, kappas);
//     Res = this->Re(rhos, U0, td.Uc(), this->d_, mus);
    
    if (useGasFilter)
    {
        Res = this->Re(rhos, U0, cloud.Uavg()[this->cell()], this->d_, mus);
    }
    else
    {
        Res = this->Re(rhos, U0, td.Uc(), this->d_, mus);
    }


    // 3. Compute heat- and momentum transfers
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Heat transfer
    // ~~~~~~~~~~~~~
//debug
// scalar dhstransReaction = dhsTrans;
// scalar oldTp = this->T_;
    // Calculate new particle temperature
    this->T_ =
        this->calcHeatTransfer
        (
            cloud,
            td,
            dt,
            Res,
            Prs,
            kappas,
            NCpW,
            Sh,
            dhsTrans,
            Sph
        );
        
// if(this->T_ < 600)
// {
//     Info<<"reaction heat: "<<dhstransReaction<<endl;
//     Info<<"convection heat: "<<dhsTrans - dhstransReaction<<endl;
//     Info<<"old Tp: "<< oldTp<< ", Tg: "<<td.Tc()<<", Tp: "<<this->T_<<endl;
//     Info<<"dMassSRCarrier: "<<dMassSRCarrier<<endl;
    
// }


    this->Cp_ = CpEff(cloud, td, pc, this->T_, idG, idL, idS);


    // Motion
    // ~~~~~~
    // moved to calcV

    


    // 4. Accumulate carrier phase source terms
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (cloud.solution().coupled())
    {
        // Transfer mass lost to carrier mass, momentum and enthalpy sources
//         forAll(YGas_, i)
//         {
//             scalar dm = np0*dMassGas[i];
//             label gid = composition.localToCarrierId(GAS, i);
//             scalar hs = composition.carrier().Hs(gid, pc, T0);
//             cloud.rhoTrans(gid)[this->cell()] += dm;
//             cloud.UTrans()[this->cell()] += dm*U0;
//             cloud.hsTrans()[this->cell()] += dm*hs;
//         }
//         forAll(YLiquid_, i)
//         {
//             scalar dm = np0*dMassLiquid[i];
//             label gid = composition.localToCarrierId(LIQ, i);
//             scalar hs = composition.carrier().Hs(gid, pc, T0);
//             cloud.rhoTrans(gid)[this->cell()] += dm;
//             cloud.UTrans()[this->cell()] += dm*U0;
//             cloud.hsTrans()[this->cell()] += dm*hs;
//         }

        forAll(dMassSRCarrier, i)
        {
            scalar dm = np0*dMassSRCarrier[i];
            scalar hs = composition.carrier().Hs(i, pc, T0);
            cloud.rhoTrans(i)[this->cell()] += dm;
            cloud.UTrans()[this->cell()] += dm*U0;
            cloud.hsTrans()[this->cell()] += dm*hs;
        }

        // Update momentum transfer
//         cloud.UTrans()[this->cell()] += np0*dUTrans;
//         cloud.UCoeff()[this->cell()] += np0*Spu;

        // Update sensible enthalpy transfer
        cloud.hsTrans()[this->cell()] += np0*dhsTrans;
        cloud.hsCoeff()[this->cell()] += np0*Sph;

        // Update radiation fields
        if (cloud.radiation())
        {
            const scalar ap = this->areaP();
            const scalar T4 = pow4(T0);
            cloud.radAreaP()[this->cell()] += dt*np0*ap;
            cloud.radT4()[this->cell()] += dt*np0*T4;
            cloud.radAreaPT4()[this->cell()] += dt*np0*ap*T4;
        }
    }
}

template<class ParcelType>
template<class TrackCloudType>
void Foam::ReactingMultiphaseParcel<ParcelType>::calcV
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
//     const bool coarseGrid = cloud.constProps().coarseGrid();
//     
//     word UScheme = cloud.solution().interpolationSchemes().lookup("U");
//     word TScheme = cloud.solution().interpolationSchemes().lookup("T");
//     
//     scalarField mp(4);
//     
    const vector& U0 = this->U_;
    const scalar T0 = this->T_;
    const scalar mass0 = this->mass();

    const bool useGasFilter = cloud.constProps().coarseGridSurfaceCombustion();
//     
//     scalar TInifinit = td.Tc();
    // Explicit momentum source for particle
    vector Su = Zero;

    // Linearised momentum source coefficient
    scalar Spu = 0.0;

    // Momentum transfer from the particle to the carrier phase
    vector dUTrans = Zero;
    
// Calc surface values
    scalar Tcavg, Ts, rhos, mus, Prs, Cpcs, kappas, Res;
    
    if (useGasFilter)
    {
        
        Tcavg = cloud.Tavg()[this->cell()];

        // Surface temperature using two thirds rule
        Ts = (2.0*T0 + Tcavg)/3.0;

        if (Ts < cloud.constProps().TMin())
        {
            if (debug)
            {
                WarningInFunction
                    << "Limiting parcel surface temperature to "
                    << cloud.constProps().TMin() <<  nl << endl;
            }

            Ts = cloud.constProps().TMin();
        }

        // Assuming thermo props vary linearly with T for small d(T)
        const scalar TRatio = Tcavg/Ts;

        rhos = cloud.rhoavg()[this->cell()]*TRatio;

        Cpcs = cloud.cpavg()[this->cell()];
        mus = cloud.muavg()[this->cell()];
        kappas = cloud.kappaavg()[this->cell()];

        Prs = Cpcs*mus/kappas;
        Prs = max(rootVSmall, Prs);
        
        Res = this->Re(rhos, U0, cloud.Uavg()[this->cell()], this->d_, mus);
    }
    else
    {
        // Calc surface values
        this->calcSurfaceValues(cloud, td, T0, Ts, rhos, mus, Prs, kappas);
        
        // Reynolds number
        Res = this->Re(rhos, U0, td.Uc(), this->d_, mus);
    }
        
    // Calculate new particle velocity
    this->U_ = this->calcVelocity(cloud, td, dt, Res, mus, mass0, Su, dUTrans, Spu);
}
// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class ParcelType>
template<class TrackCloudType>
void Foam::ReactingMultiphaseParcel<ParcelType>::calcDevolatilisation
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt,
    const scalar age,
    const scalar Ts,
    const scalar d,
    const scalar T,
    const scalar mass,
    const scalar mass0,
    const scalarField& YGasEff,
    const scalarField& YLiquidEff,
    const scalarField& YSolidEff,
    label& canCombust,
    scalarField& dMassDV,
    scalarField& dMassDVSoild,
    scalar& Sh,
    scalar& N,
    scalar& NCpW,
    scalarField& Cs
) const
{
    // Check that model is active
    if (!cloud.devolatilisation().active())
    {
        if (canCombust != -1)
        {
            canCombust = 1;
        }
        return;
    }

    // Initialise demand-driven constants
    (void)cloud.constProps().TDevol();
    (void)cloud.constProps().LDevol();

    // Check that the parcel temperature is within necessary limits for
    // devolatilisation to occur
    if (T < cloud.constProps().TDevol() || canCombust == -1)
    {
        return;
    }

    typedef typename TrackCloudType::reactingCloudType reactingCloudType;
    const CompositionModel<reactingCloudType>& composition =
        cloud.composition();


    // Total mass of volatiles evolved
    cloud.devolatilisation().calculate
    (
        dt,
        age,
        mass0,
        mass,
        T,
        YGasEff,
        YLiquidEff,
        YSolidEff,
        canCombust,
        dMassDV,
        dMassDVSoild
    );

    scalar dMassTot = sum(dMassDV);

    cloud.devolatilisation().addToDevolatilisationMass
    (
        this->nParticle_*dMassTot
    );

    Sh -= dMassTot*cloud.constProps().LDevol()/dt;

    // Update molar emissions
    if (cloud.heatTransfer().BirdCorrection())
    {
        // Molar average molecular weight of carrier mix
        const scalar Wc = max(small, td.rhoc()*RR*td.Tc()/td.pc());

        // Note: hardcoded gaseous diffusivities for now
        // TODO: add to carrier thermo
        const scalar beta = sqr(cbrt(15.0) + cbrt(15.0));

        forAll(dMassDV, i)
        {
            const label id = composition.localToCarrierId(GAS, i);
            const scalar Cp = composition.carrier().Cp(id, td.pc(), Ts);
            const scalar W = composition.carrier().Wi(id);
            const scalar Ni = dMassDV[i]/(this->areaS(d)*dt*W);

            // Dab calc'd using API vapour mass diffusivity function
            const scalar Dab =
                3.6059e-3*(pow(1.8*Ts, 1.75))
               *sqrt(1.0/W + 1.0/Wc)
               /(td.pc()*beta);

            N += Ni;
            NCpW += Ni*Cp*W;
            Cs[id] += Ni*d/(2.0*Dab);
        }
    }
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::ReactingMultiphaseParcel<ParcelType>::calcSurfaceReactions
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt,
    const scalar d,
    const scalar di,
    const scalar T,
    const scalar mass,
    const label canCombust,
    const scalar N,
    const scalarField& YMix,
    const scalarField& YGas,
    const scalarField& YLiquid,
    const scalarField& YSolid,
    scalarField& dMassSRGas,
    scalarField& dMassSRLiquid,
    scalarField& dMassSRSolid,
    scalarField& dMassSRCarrier,
    scalar& Sh,
    scalar& dhsTrans,
    scalar& QComb,
    const scalar Re,
    const scalar Tc,
    const scalar rhoc,
    const scalar muc,
    const scalar charSphericity,
    scalarField& dHeatSRCarrier,
    const label particleShape,
    const scalarField& Csavg,
    const bool coarseGrid,
    const scalar deq,
    const bool heatRatio,
    const bool limitedCombustion
) const
{
    // Check that model is active
    if (!cloud.surfaceReaction().active())
    {
        return;
    }

    // Initialise demand-driven constants
    (void)cloud.constProps().hRetentionCoeff();
    (void)cloud.constProps().TMax();

    // Check that model is active
    if (canCombust != 1)
    {
        return;
    }
    
    scalar hReaction;
    
    // Update surface reactions
    hReaction = cloud.surfaceReaction().calculate
    (
        dt,
        this->cell(),
        d,
        T,
        Tc,
        td.pc(),
        rhoc,
        mass,
        YGas,
        YLiquid,
        YSolid,
        YMix,
        N,
        dMassSRGas,
        dMassSRLiquid,
        dMassSRSolid,
        dMassSRCarrier,
        di,
        muc,
        Re,
        charSphericity,
        dHeatSRCarrier,
        particleShape,
        Csavg,
        coarseGrid,
        deq
    );
    

    cloud.surfaceReaction().addToSurfaceReactionMass
    (
        this->nParticle_
       *(sum(dMassSRGas) + sum(dMassSRLiquid) + sum(dMassSRSolid))
    );

    const scalar xsi = min(T/cloud.constProps().TMax(), 1.0);
    const scalar coeff =
        (1.0 - xsi*xsi)*cloud.constProps().hRetentionCoeff();

    Sh += coeff*hReaction/dt;
// Info<<" cloud.constProps().TMax() :"<< cloud.constProps().TMax() <<", xsi :"<< xsi <<", cloud.constProps().hRetentionCoeff() :"<< cloud.constProps().hRetentionCoeff() <<", coeff :"<< coeff << endl;
    dhsTrans += (1.0 - coeff)*hReaction;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::ReactingMultiphaseParcel<ParcelType>::ReactingMultiphaseParcel
(
    const ReactingMultiphaseParcel<ParcelType>& p
)
:
    ParcelType(p),
    YGas_(p.YGas_),
    YLiquid_(p.YLiquid_),
    YSolid_(p.YSolid_),
    canCombust_(p.canCombust_)
{}


template<class ParcelType>
Foam::ReactingMultiphaseParcel<ParcelType>::ReactingMultiphaseParcel
(
    const ReactingMultiphaseParcel<ParcelType>& p,
    const polyMesh& mesh
)
:
    ParcelType(p, mesh),
    YGas_(p.YGas_),
    YLiquid_(p.YLiquid_),
    YSolid_(p.YSolid_),
    canCombust_(p.canCombust_)
{}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "ReactingMultiphaseParcelIO.C"

// ************************************************************************* //
