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

#include "renPaSR.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::combustionModels::renPaSR<ReactionThermo>::renPaSR
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
    Ce_(this->coeffs().getScalar("Ce")),
    kappa_
    (
        IOobject
        (
            thermo.phasePropertyName(typeName + ":kappa"),
            this->mesh().time().timeName(),
        this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("kappa", dimless, 0.0)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::combustionModels::renPaSR<ReactionThermo>::~renPaSR()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ReactionThermo>
void Foam::combustionModels::renPaSR<ReactionThermo>::correct()
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

        // Declare scalars used when scaling reaction rates
        scalar tc, tcSpecies, tk, tp, tm;   

        forAll(epsilon, celli)
        {
            // Ren model: compute the chemical time scale
            // based on the shortest species residence time
            tc = GREAT;
            forAll(Y_, si)
            {
                scalarField& Ri = this->chemistryPtr_->RR(si);
                if (Y_[si][celli] > 1e-16 && Ri[celli] < -1e-16)
                {
                    tcSpecies = Y_[si][celli]/mag(Ri[celli]);
                    if (tcSpecies < tc)
                    {
                        tc = tcSpecies;
                    }
                }
            }

            tc = Ce_*tc;

            // Estimate the Kolmogorov time-scale
            tk = sqrt(max(muEff[celli]/rho[celli]/(epsilon[celli] + SMALL), 0));

            tp = max(tke[celli]/(epsilon[celli] + SMALL), 0);

            // Compute the mixing time scale
            tm = Cmix_*sqrt(tk*tp);

            // If the mixing time scale is very short or
            // the chemical time scale very long, assume
            // perfect mixing. Otherwise, use partial mixing.
            if (tm > SMALL && tc < GREAT) 
            {
                // Reaction rate scaling factor, a.k.a. the
                // volume fraction of reacting fine structures
                kappa_[celli] = tc/(tc + tm);     
            } 
            else
            {
                // Perfect mixing
                kappa_[celli] = 1.0;
            }
        } 
    }
}


template<class ReactionThermo>
Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::renPaSR<ReactionThermo>::R(volScalarField& Y) const
{
    return kappa_*laminar<ReactionThermo>::R(Y);
}


template<class ReactionThermo>
Foam::tmp<Foam::volScalarField>
Foam::combustionModels::renPaSR<ReactionThermo>::Qdot() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            this->thermo().phasePropertyName(typeName + ":Qdot"),
            kappa_*laminar<ReactionThermo>::Qdot()
        )
    );
}


template<class ReactionThermo>
bool Foam::combustionModels::renPaSR<ReactionThermo>::read()
{
    if (laminar<ReactionThermo>::read())
    {
        this->coeffs().readEntry("Cmix", Cmix_);
        return true;
    }

    return false; 
}

// ************************************************************************* //
