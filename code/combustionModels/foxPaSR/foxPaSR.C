/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "foxPaSR.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::combustionModels::foxPaSR<ReactionThermo>::foxPaSR
(
    const word& modelType,
    ReactionThermo& thermo,
    const compressibleTurbulenceModel& turb,
    const word& combustionProperties
)
:
    laminar<ReactionThermo>(modelType, thermo, turb, combustionProperties),
    Cmix_(this->coeffs().getScalar("Cmix")),
    Y_(this->thermo().composition().Y()),
    reactions_
    (
        dynamic_cast<const reactingMixture<gasHThermoPhysics>&>(this->thermo())
    ),
    modRR_(this->chemistryPtr_->nSpecie())
{
    forAll(modRR_, i)
    {
        modRR_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    thermo.phasePropertyName(
                        typeName + ":modRR_" + Y_[i].name()
                    ),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh(),
                dimensionedScalar("modRR",  dimensionSet(1, -3, -1, 0, 0), 0.0)
            )
        );
        Info << "Creating field " << modRR_[i].name() << endl;
    } 
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::combustionModels::foxPaSR<ReactionThermo>::~foxPaSR()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ReactionThermo>
void Foam::combustionModels::foxPaSR<ReactionThermo>::correct()
{
    if (this->active())
    {
        // Update chemistry
        laminar<ReactionThermo>::correct();

        // Reference to turbulent dissipation
        tmp<volScalarField> tepsilon(this->turbulence().epsilon());
        const scalarField& epsilon = tepsilon();

        // Reference to effective viscosity
        tmp<volScalarField> tmuEff(this->turbulence().muEff());
        const scalarField& muEff = tmuEff();

        // Reference to turbulent kinetic energy
        tmp<volScalarField> ttke(this->turbulence().k());
        const scalarField& tke = ttke();

        // Reference to density
        tmp<volScalarField> trho(this->rho());
        const scalarField& rho = trho();

        // Reset modified production rates of all species in all cells
        forAll(modRR_, i)
        {
            forAll(epsilon, celli)
            {
                modRR_[i][celli] = 0.0;
            }
        }

        // Declare scalars used when scaling reaction rates
        scalar tc, tk, tp, tm, kappa, tcSpecies;

        // Loop over all reactions
        forAll(reactions_, ri)
        {
            // Select the first species on the LHS as the reference species
            label refSpecies = reactions_[ri].lhs()[0].index;

            // Retrieve the molecular weight of the reference species
            const scalar& refW =
                this->thermo().composition().W(refSpecies);

            // Compute the reference reaction rate [kmol/s] of the
            // reaction. In other words, the production/consumption
            // rate of the first species on the LHS. We assume that
            // its stoichiometric coefficient is unity.
            scalarField refRR =
                this->chemistryPtr_->calculateRR(ri, refSpecies)/refW;

            // Loop over all cells
            forAll(epsilon, celli)
            {
                // Declare references to local cell values
                const scalar& muC = muEff[celli];
                const scalar& rhoC = rho[celli];
                const scalar& epsC = epsilon[celli];
                const scalar& tkeC = tke[celli];
                const scalar& refRRC = refRR[celli];

                // Initialize the local chemical time scale
                tc = GREAT;

                // If the net reaction rate is positive, the reaction
                // is treated as an irreversible reaction with reactants
                // on the LHS and products on the RHS. If it is negative,
                // the reaction is treated as an irreversible reaction
                // with reactants on the RHS and products on the RHS.
                if (refRRC > SMALL)
                {
                    // Loop over all species on LHS
                    forAll(reactions_[ri].lhs(), si)
                    {
                        // Species index
                        const label& i = reactions_[ri].lhs()[si].index;

                        // Stoichiometric coefficient
                        const scalar& stoichCoeff =
                            reactions_[ri].lhs()[si].stoichCoeff;

                        // Molecular weight
                        const scalar& W =
                            this->thermo().composition().W(i);

                        // Compute the residence time of this species
                        // (concentration / consumption)
                        tcSpecies =
                            mag(rhoC*Y_[i][celli]/(refRRC*W*stoichCoeff));

                        // The shortest residence time is the chemical
                        // time scale
                        if (tcSpecies < tc)
                        {
                            tc = tcSpecies;
                        }
                    }
                }
                else if (refRRC < -SMALL)
                {
                    // Loop over all species on RHS
                    forAll(reactions_[ri].rhs(), si)
                    {
                        // Species index
                        const label& i = reactions_[ri].rhs()[si].index;

                        // Stoichiometric coefficient
                        const scalar& stoichCoeff =
                            reactions_[ri].rhs()[si].stoichCoeff;

                        // Molecular weight
                        const scalar& W = this->thermo().composition().W(i);

                        // Compute the residence time of this species
                        // (concentration / consumption)
                        tcSpecies =
                            mag(rhoC*Y_[i][celli]/(refRRC*W*stoichCoeff));

                        // The shortest residence time is the chemical
                        // time scale
                        if (tcSpecies < tc)
                        {
                            tc = tcSpecies;
                        }
                    }
                }

                // Estimate the Kolmogorov time scale
                tk =
                    sqrt(max(muC/rhoC/(epsC + SMALL), 0));

                tp =
                    max(tkeC/(epsC + SMALL), 0);

                // Compute mixing time scale
                tm = sqrt(tk*tp);

                // If the mixing time scale is very short or
                // the chemical time scale very long, assume
                // perfect mixing. Otherwise, use partial mixing.
                if (tm > SMALL && tc < GREAT)
                {
                    // Reaction rate scaling factor, a.k.a. the
                    // volume fraction of reacting fine structures
                    kappa = tc/(tc + tm);
                }
                else
                {
                    // Perfect mixing
                    kappa = 1.0;
                }
                // Update the modified production rates of involved
                // species by first looping over LHS species
                // and then over RHS species
                forAll(reactions_[ri].lhs(), si)
                {
                    // Species index
                    const label& i = reactions_[ri].lhs()[si].index;

                    // Stoichiometric coefficient
                    const scalar& stoichCoeff =
                        reactions_[ri].lhs()[si].stoichCoeff;

                    // Molecular weight
                    const scalar& W = this->thermo().composition().W(i);

                    // Compute the production rate of the species
                    // and scale it with kappa. Add it to the total
                    // modified production rate of the species.
                    // Note that a positive refRRC means that
                    // species on the LHS are consumed.
                    modRR_[i][celli] -=
                        kappa*refRRC*W*stoichCoeff;
                }
                forAll(reactions_[ri].rhs(), si)
                {
                    // Species index
                    const label& i = reactions_[ri].rhs()[si].index;

                    // Stoichiometric coefficient
                    const scalar& stoichCoeff =
                        reactions_[ri].rhs()[si].stoichCoeff;

                    // Molecular weight
                    const scalar& W = this->thermo().composition().W(i);

                    // Compute the production rate of the species
                    // and scale it with kappa. Add it to the total
                    // modified production rate of the species.
                    // Note that a positive refRRC means that 
                    // species on the RHS are produced.
                    modRR_[i][celli] +=
                        kappa*refRRC*W*stoichCoeff;
                }
            }
        }

        // Overwrite previously computed production rates with modified
        // production rates
        forAll(modRR_, i)
        {
            this->chemistryPtr_->RR(i) = modRR_[i];
        }
    }
}


template<class ReactionThermo>
Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::foxPaSR<ReactionThermo>::R(volScalarField& Y) const
{
    return laminar<ReactionThermo>::R(Y);
}


template<class ReactionThermo>
Foam::tmp<Foam::volScalarField>
Foam::combustionModels::foxPaSR<ReactionThermo>::Qdot() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            this->thermo().phasePropertyName(typeName + ":Qdot"),
            laminar<ReactionThermo>::Qdot()
        )
    );
}


template<class ReactionThermo>
bool Foam::combustionModels::foxPaSR<ReactionThermo>::read()
{
    if (laminar<ReactionThermo>::read())
    {
        this->coeffs().readEntry("Cmix", Cmix_);
        return true;
    }

    return false; 
}

// ************************************************************************* //
