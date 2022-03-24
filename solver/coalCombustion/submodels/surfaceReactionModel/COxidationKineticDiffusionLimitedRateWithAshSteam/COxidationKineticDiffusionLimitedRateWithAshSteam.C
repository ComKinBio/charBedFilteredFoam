/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "COxidationKineticDiffusionLimitedRateWithAshSteam.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::COxidationKineticDiffusionLimitedRateWithAshSteam<CloudType>::
COxidationKineticDiffusionLimitedRateWithAshSteam
(
    const dictionary& dict,
    CloudType& owner
)
:
    SurfaceReactionModel<CloudType>(dict, owner, typeName),
    COmega_(readScalar(this->coeffDict().lookup("COmiga"))),
    ep3_(readScalar(this->coeffDict().lookup("ep3"))),
    C1_O2_(readScalar(this->coeffDict().lookup("C1_O2"))),
    C1_H2O_(readScalar(this->coeffDict().lookup("C1_H2O"))),
    C1_CO2_(readScalar(this->coeffDict().lookup("C1_CO2"))),
    C1_H2_(readScalar(this->coeffDict().lookup("C1_H2"))),
    C2_1_(readScalar(this->coeffDict().lookup("C2_1"))),
    C2_2_(readScalar(this->coeffDict().lookup("C2_2"))),
    C2_3_(readScalar(this->coeffDict().lookup("C2_3"))),
    C2_4_(readScalar(this->coeffDict().lookup("C2_4"))),
    E_1_(readScalar(this->coeffDict().lookup("E1"))),
    E_2_(readScalar(this->coeffDict().lookup("E2"))),
    E_3_(readScalar(this->coeffDict().lookup("E3"))),
    E_4_(readScalar(this->coeffDict().lookup("E4"))),
    CsLocalId_(-1),
    AshLocalId_(-1),
    O2GlobalId_(owner.composition().carrierId("O2")),
    CO2GlobalId_(owner.composition().carrierId("CO2")),
    COGlobalId_(owner.composition().carrierId("CO")),
    WC_(0.0),
    WCO2_(0.0),
    WO2_(0.0),
    HcCO2_(0.0),
    HcCO_(0.0)
{
    // Determine Cs ids
    label idSolid = owner.composition().idSolid();
    CsLocalId_ = owner.composition().localId(idSolid, "C");
    AshLocalId_ = owner.composition().localId(idSolid, "Ash");

    // Set local copies of thermo properties
    WO2_ = owner.thermo().carrier().Wi(O2GlobalId_);
    WCO2_ = owner.thermo().carrier().Wi(CO2GlobalId_);
    WC_ = WCO2_ - WO2_;
    
    HcCO2_ = owner.thermo().carrier().Hc(CO2GlobalId_);
    HcCO_ = owner.thermo().carrier().Hc(COGlobalId_);

    const scalar YCloc = owner.composition().Y0(idSolid)[CsLocalId_];
    const scalar YSolidTot = owner.composition().YMixture0()[idSolid];
    Info<< "    C(s): particle mass fraction = " << YCloc*YSolidTot << endl;
}


template<class CloudType>
Foam::COxidationKineticDiffusionLimitedRateWithAshSteam<CloudType>::
COxidationKineticDiffusionLimitedRateWithAshSteam
(
    const COxidationKineticDiffusionLimitedRateWithAshSteam<CloudType>& srm
)
:
    SurfaceReactionModel<CloudType>(srm),
    //Sb_(srm.Sb_),
    COmega_(srm.COmega_),
    ep3_(srm.ep3_),
    C1_O2_(srm.C1_O2_),
    C1_H2O_(srm.C1_H2O_),
    C1_CO2_(srm.C1_CO2_),
    C1_H2_(srm.C1_H2_),
    C2_1_(srm.C2_1_),
    C2_2_(srm.C2_2_),
    C2_3_(srm.C2_3_),
    C2_4_(srm.C2_4_),
    E_1_(srm.E_1_),
    E_2_(srm.E_2_),
    E_3_(srm.E_3_),
    E_4_(srm.E_4_),
    CsLocalId_(srm.CsLocalId_),
    AshLocalId_(srm.AshLocalId_),
    O2GlobalId_(srm.O2GlobalId_),
    CO2GlobalId_(srm.CO2GlobalId_),
    COGlobalId_(srm.COGlobalId_),
    WC_(srm.WC_),
    WCO2_(srm.WCO2_),
    WO2_(srm.WO2_),
    HcCO2_(srm.HcCO2_),
    HcCO_(srm.HcCO_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::COxidationKineticDiffusionLimitedRateWithAshSteam<CloudType>::
~COxidationKineticDiffusionLimitedRateWithAshSteam()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::COxidationKineticDiffusionLimitedRateWithAshSteam<CloudType>::calculate
(
    const scalar dt,
    const label celli,
    const scalar d,
    const scalar T,
    const scalar Tc,
    const scalar pc,
    const scalar rhoc,
    const scalar mass,
    const scalarField& YGas,
    const scalarField& YLiquid,
    const scalarField& YSolid,
    const scalarField& YMixture,
    const scalar N,
    scalarField& dMassGas,
    scalarField& dMassLiquid,
    scalarField& dMassSolid,
    scalarField& dMassSRCarrier,
    const scalar di,
    const scalar muc,
    const scalar Rec,
    const scalar Xi,//changed into charSphericity
    scalarField& dHeatSRCarrier,
    const label particleShape,
    const scalarField& Csavg,
    const bool coarseGrid,
    const scalar deq
    
) const
{
   // Fraction of remaining combustible material
    const label idSolid = CloudType::parcelType::SLD;
    const scalar fComb = YMixture[idSolid]*YSolid[CsLocalId_];

    // Surface combustion active combustible fraction is consumed
    if (fComb < SMALL)
    {
        return 0.0;
    }

    const SLGThermo& thermo = this->owner().thermo();
    
    scalar YO2, YCO2, ebed, epsilon_ash;
    
     // Local mass fraction of O2 in the carrier phase
    if (coarseGrid)
    {
        YO2 = Csavg[0];
        YCO2 = Csavg[2];
    }
    else
    {
        YO2 = thermo.carrier().Y(O2GlobalId_)[celli];
        YCO2 = thermo.carrier().Y(CO2GlobalId_)[celli];
    }
    
    ebed = Csavg[4];
    
    epsilon_ash = Csavg[5];

    // Local mole concentration of O2 in the carrier phase,            
    const scalar C_O2 = YO2*rhoc/WO2_;
    const scalar C_CO2 = YCO2*rhoc/WCO2_;

    // Diffusion rate coefficient
    const scalar D0_O2 = C1_O2_*pow(0.5*(T + Tc), 2.0);
    const scalar D0ash_O2 = C1_O2_*pow(T, 2.0);
    const scalar D0_CO2 = C1_O2_*pow(0.5*(T + Tc), 2.0);
    const scalar D0ash_CO2 = C1_O2_*pow(T, 2.0);

    
    /*const scalar D0_CO2 = C1_CO2_*pow(0.5*(T + Tc), 1.81);
    const scalar D0ash_CO2 = C1_CO2_*pow(T, 1.81);
    const scalar D0_H2O = C1_H2O_*pow(0.5*(T + Tc), 1.81);
    const scalar D0ash_H2O = C1_H2O_*pow(T, 1.81);
    const scalar D0_H2 = C1_H2_*pow(0.5*(T + Tc), 1.81);
    const scalar D0ash_H2 = C1_H2_*pow(T, 1.81);*/

    // Schmidt number
    const scalar Sc_O2 = muc/(rhoc*D0_O2);
    const scalar Sc_CO2 = muc/(rhoc*D0_CO2);


    // Sherwood number
//     const scalar Sh_O2 = 2.0+1.1*pow(Sc_O2, 0.3333)*pow(Rec, 0.6);
//     const scalar Sh_CO2 = 2.0+1.1*pow(Sc_CO2, 0.3333)*pow(Rec, 0.6);
//     const scalar Sh_H2O = 2.0+1.1*pow(Sc_H2O, 0.3333)*pow(Rec, 0.6);
//     const scalar Sh_H2 = 2.0+1.1*pow(Sc_H2, 0.3333)*pow(Rec, 0.6);kappas/(this->d_)*(2.0*e_bed+0.69*pow(Res/e_bed,0.5)*cbrt(Prs));
    scalar ReTerm = 0.69*sqrt(Rec/ebed);
    const scalar Sh_O2 = 2.0*ebed+ReTerm*cbrt(Sc_O2);
    const scalar Sh_CO2 = 2.0*ebed+ReTerm*cbrt(Sc_CO2);

    
//     scalar ReTerm = 0.29*pow(ebed,-0.81)*pow(Rec,0.73);
//     const scalar Sh_O2 = 1.77+ReTerm*sqrt(Sc_O2);
//     const scalar Sh_CO2 = 1.77+ReTerm*sqrt(Sc_CO2);
//     const scalar Sh_H2O = 1.77+ReTerm*sqrt(Sc_H2O);
//     const scalar Sh_H2 = 1.77+ReTerm*sqrt(Sc_H2);


    // Mass convection coefficient hma
    const scalar hma_O2 = 2.0*D0_O2*Sh_O2/(d);
    const scalar hma_CO2 = 2.0*D0_CO2*Sh_CO2/(d);


    // Mass diffusion coefficient through the ash layer hmi
    scalar deltad = max((d-di)/2.0, small);
    const scalar hmi_O2 = D0ash_O2*pow(epsilon_ash, 2.0)/deltad;
    const scalar hmi_CO2 = D0ash_CO2*pow(epsilon_ash, 2.0)/deltad;
 

    // General mass transfer rate for combustion
    const scalar beta_d_O2 = hmi_O2*hma_O2/(hmi_O2+hma_O2);
    const scalar beta_d_CO2 = hmi_CO2*hma_CO2/(hmi_CO2+hma_CO2);


    // Kinetic rate beta_r
    const scalar beta_r_O2 = C2_1_*T*exp(-E_1_/T);
    const scalar beta_r_CO2 = C2_2_*T*exp(-E_2_/T);

// if (T>700)
// {
//     Info<<"T : "<<T<<endl;
// Info<<"hma_O2: "<<hma_O2<<endl;
// Info<<"hmi_O2: "<<hmi_O2<<endl;
// Info<<"beta_d_O2: "<<beta_d_O2<<endl;
// Info<<"beta_r_O2: "<<beta_r_O2<<endl; 
// Info<<"-------------------------------------"<<endl;
// }
    // Correlation for the formation of CO and CO2 omega_c
    const scalar omega_c = 2.0*(1.0+4.3*exp(-COmega_/T))/(2.0+4.3*exp(-COmega_/T));
    
    // Particle surface area for the char layer for cylinder particle shape = 1
    scalar Ap;
    if (particleShape == 1)
    {
        Ap = constant::mathematical::pi*sqr(di)/Xi;
    }
    else
    {
        Ap = constant::mathematical::pi*sqr(di);
    }

    // Change in C mass [kg]
    scalar dmC_O2 = WC_*omega_c*C_O2*Ap*(beta_r_O2*beta_d_O2/(beta_r_O2+beta_d_O2))*dt;
    scalar dmC_CO2 = WC_*C_CO2*Ap*(beta_r_CO2*beta_d_CO2/(beta_r_CO2+beta_d_CO2))*dt;


    // Limit mass transfer by availability of C
    dmC_O2 = min(mass*fComb, dmC_O2);
    dmC_CO2 = min(mass*fComb, dmC_CO2);


    // Molar consumption
    const scalar dOmega_O2 = dmC_O2/omega_c/WC_;
    const scalar dOmega_CO2 = dmC_CO2/WC_;

    // Change in O2 mass [kg]
    const scalar dmO2 = dOmega_O2*WO2_;

    // Mass of newly created CO2 [kg]
    const scalar dmCO2 = (dOmega_O2*(2.0-omega_c)-dOmega_CO2)*(WC_ + WO2_);

    // Mass of newly created CO [kg]
    const scalar dmCO = (dOmega_O2*(2.0*(omega_c-1.0))+dOmega_CO2*2.0)*(WC_ + WO2_/2.0);

    
//     Info<<"dmO2: "<<dmO2<<endl;
//     Info<<"dmCO2: "<<dmCO2<<endl;
//     Info<<"dmCO: "<<dmCO<<endl;
//     Info<<"dmC: "<<(dOmega_O2*omega_c+dOmega_CO2)*WC_<<endl;

    // Update local particle C mass
    dMassSolid[CsLocalId_] = (dOmega_O2*omega_c+dOmega_CO2)*WC_;
    dMassSolid[AshLocalId_] = -(dOmega_O2*omega_c+dOmega_CO2)*WC_*deq;//deq is changed to AshInFuel
    
    // Update carrier O2 and CO2 mass
    dMassSRCarrier[O2GlobalId_] -= dmO2;
    dMassSRCarrier[CO2GlobalId_] += dmCO2;
    dMassSRCarrier[COGlobalId_] += dmCO;

    //dMassSRCarrier[H2GlobalId_] = 0.0;
    const scalar LHV_index_CO2 = 0.0;               /* per definition */
    const scalar LHV_index_CO = 10.25e6;            /* J/kg D&D p. 26 */
    const scalar LHV_index_C = 32.79e6;         /* J/kg Alberto's code */

//     const scalar HsC = thermo.solids().properties()[CsLocalId_].Hs(T);

//     Correlation for heatRatio CO and CO2
    const scalar COCO2Ratio = 4.3*exp(-COmega_/T);
    const scalar heatRatio = (COCO2Ratio+0.3)/(COCO2Ratio+1.0);

    // carrier sensible enthalpy exchange handled via change in mass
//  Info<<"dmC:   "<<dmC<<endl; 
    // Heat of reaction [J]
    const scalar Heat_O2 = (dOmega_O2*(2.0-omega_c)*(WC_ + WO2_))*LHV_index_CO2 + (dOmega_O2*(2.0*(omega_c-1.0))*(WC_ + WO2_/2.0))*LHV_index_CO - dmC_O2*LHV_index_C;
    const scalar Heat_CO2 = dOmega_CO2*(WC_ + WO2_/2.0)*2.0*LHV_index_CO-dOmega_CO2*LHV_index_C*WC_;
    
    dHeatSRCarrier[0] = Heat_O2;
    dHeatSRCarrier[1] = Heat_CO2;

//     const scalar HsC = thermo.solids().properties()[CsLocalId_].Hs(T); 
//debug information
    
//     dHeatSRCarrier[4] = beta_r_O2;
//     dHeatSRCarrier[5] = beta_d_O2;
//     dHeatSRCarrier[6] = hma_O2;
//     dHeatSRCarrier[7] = Sh_O2;
    

/*    if (T > 1000.0 || T < 0.0)
{
    Info<<"Sh_O2 "<<Sh_O2<<endl;    
    
    scalar Sh_O2_old = 2.0+0.6*pow(Sc_O2, 0.3333)*pow(Rec, 0.5);
    
    Info<<"Sh_O2_old  "<<Sh_O2_old <<endl;
    
    Info<<"WC_  "<<WC_ <<endl; 
    
    Info<<"heatRatio  "<<heatRatio <<endl; 
} */   
    
    return -(Heat_O2 + Heat_CO2);
}


// ************************************************************************* //
