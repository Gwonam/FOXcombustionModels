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
    Foam::combustionModels::renPaSR

Group
    grpCombustionModels

Description
    Partially stirred reactor turbulent combustion model.

    This model calculates finite rates based on both turbulent and chemical
    time scales.  Depending on mesh resolution, the Cmix parameter can be used
    to adjust the turbulent mixing time scale. The chemical time scale is
    computed locally based on the shortest species residence time, following
    Ren & Golding (10.1016/j.combustflame.2011.02.018).
    The turbulent mixing time scale is a harmonic average between the
    Kolmogorov time scale and the time scale of sub-grid velocity stretch.

SourceFiles
    renPaSR.C

\*---------------------------------------------------------------------------*/

#ifndef renPaSR_H
#define renPaSR_H

#include "laminar.H"
#include "thermoPhysicsTypes.H"
#include "StandardChemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                           Class renPaSR Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo>
class renPaSR
:
    public laminar<ReactionThermo>
{
    // Private data

        //- Mixing parameter
        scalar Cmix_;  

        //- Species mass fractions
        const PtrList<volScalarField>& Y_;

        //- Chemical time scale parameter
        scalar Ce_;

        //- Fine structure volume fraction
        volScalarField kappa_;

    // Private Member Functions

        //- No copy construct
        renPaSR(const renPaSR&) = delete;

        //- No copy assignment
        void operator=(const renPaSR&) = delete;


public:

    //- Runtime type information
    TypeName("renPaSR");


    // Constructors

        //- Construct from components
        renPaSR
        (
            const word& modelType,
            ReactionThermo& thermo,
            const compressibleTurbulenceModel& turb,
            const word& combustionProperties
  	    );


    //- Destructor
    virtual ~renPaSR();


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
    #include "renPaSR.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
