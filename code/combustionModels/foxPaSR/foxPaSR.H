/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::combustionModels::foxPaSR

Group
    grpCombustionModels

Description
    Partially stirred reactor turbulent combustion model with individual
    reaction scaling.

    This model calculates finite rates based on both turbulent and chemical
    time scales.  Depending on mesh resolution, the Cmix parameter can be used
    to adjust the turbulent mixing time scale. A separate chemical time scale
    is computed for each reaction, and each reaction is scaled separately.
    The turbulent mixing time scale is a harmonic average between the
    Kolmogorov time scale and the time scale of sub-grid velocity stretch.

SourceFiles
    foxPaSR.C

\*---------------------------------------------------------------------------*/

#ifndef foxPaSR_H
#define foxPaSR_H

#include "laminar.H"
#include "thermoPhysicsTypes.H"
#include "StandardChemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                           Class foxPaSR Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo>
class foxPaSR
:
    public laminar<ReactionThermo>
{
    // Private data

        //- Mixing parameter
        scalar Cmix_;  

        //- Number of species
        //label nSpecie_;

        //- Species mass fractions
        const PtrList<volScalarField>& Y_;

        //- List of reactions
        const PtrList<Reaction<gasHThermoPhysics>>& reactions_;

        //- Scaling factors
        PtrList<volScalarField> modRR_;

    // Private Member Functions

        //- No copy construct
        foxPaSR(const foxPaSR&) = delete;

        //- No copy assignment
        void operator=(const foxPaSR&) = delete;


public:

    //- Runtime type information
    TypeName("foxPaSR");


    // Constructors

        //- Construct from components
        foxPaSR
        (
            const word& modelType,
            ReactionThermo& thermo,
            const compressibleTurbulenceModel& turb,
            const word& combustionProperties
  	    );


    //- Destructor
    virtual ~foxPaSR();


    // Member Functions

        //- Correct combustion rate
        virtual void correct();

        //- Fuel consumption rate matrix
        virtual tmp<fvScalarMatrix> R(volScalarField& Y) const;

        //- Heat release rate [kg/m/s3]
        virtual tmp<volScalarField> Qdot() const;

        //- Update properties from given dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "foxPaSR.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
