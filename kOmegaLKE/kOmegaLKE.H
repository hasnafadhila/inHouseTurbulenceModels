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

Class
    Foam::RASModels::kOmegaLKE

Group
    grpRASTurbulence

Description
    This is an in-house implementation of the k-omega LKE model. The model uses
    a new approach for modelling the laminar kinetic energy using a laminar
    diffusion "eddy" viscosity, coupled with the classic k-omega transport
    equations.

    References:
    \verbatim
        Medina, H., Beechook, A., Fadhila, H., Aleksandrova, S., & Benjamin,
        S. (2018). A novel laminar kinetic energy model for the prediction of
        pretransitional velocity fluctuations and boundary layer transition.
        International Journal of Heat and Fluid Flow, 69, 150-163.
    \endverbatim

    The default model coefficients are
    \verbatim
        kOmegaLKECoeffs
        {
            Cmu         0.09;  // Equivalent to betaStar
            alpha       0.52;
            beta        0.072;
            alphak      0.5;
            alphaOmega  0.5;
        }
    \endverbatim

SourceFiles
    kOmegaLKE.C

\*---------------------------------------------------------------------------*/

#ifndef kOmegaLKE_H
#define kOmegaLKE_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class kOmegaLKE Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class kOmegaLKE
:
    public eddyViscosity<RASModel<BasicTurbulenceModel>>
{

protected:

    // Protected data

        // Model coefficients

            dimensionedScalar Cmu_;
            dimensionedScalar cOmega1_;
            dimensionedScalar cOmega2_;
            dimensionedScalar sigmaK_;
            dimensionedScalar sigmaOmega_;
            dimensionedScalar sigmaKl_;
            dimensionedScalar tuInf_;
            dimensionedScalar cCrit_;
            dimensionedScalar cSS_;
            dimensionedScalar cV_;


        // Fields

            volScalarField k_;
            volScalarField omega_;
            volScalarField kl_;

            const volScalarField& y_;


    // Protected Member Functions

        virtual void correctNut
        (
            const volScalarField& S2
            //const volScalarField& Omega
        );

        virtual void correctNut();

        tmp<volScalarField> ReLambda() const;
        tmp<volScalarField> ReUpsilon() const;
        tmp<volScalarField> fSS(const volScalarField& Omega) const;
        tmp<volScalarField> fv() const;



public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("kOmegaLKE");


    // Constructors

        //- Construct from components
        kOmegaLKE
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~kOmegaLKE()
    {}


    // Member Functions

        //- Read RASProperties dictionary
        virtual bool read();
/*
        //- Return the effective diffusivity for k
        //tmp<volScalarField> DkEff(const volScalarField& alphaStar) const
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DkEff",
                    alphaK_*this->nut_ + this->nu()
                    //alphaStar*alphaK_*this->nut_ + this->nu()
                )
            );
        }

        //- Return the effective diffusivity for omega
        //tmp<volScalarField> DomegaEff(const volScalarField& alphaStar) const
        tmp<volScalarField> DomegaEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DomegaEff",
                    alphaOmega_*this->nut_ + this->nu()
                    //alphaStar*alphaOmega_*this->nut_ + this->nu()
                )
            );
        }
*/

//- Return the effective diffusivity for k
        tmp<volScalarField> DkEff(const volScalarField& nutDiff) const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DkEff",
                    nutDiff + this->nu()
                )
            );
        }

        //- Return the effective diffusivity for omega
        tmp<volScalarField> DomegaEff(const volScalarField& nutDiff) const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DomegaEff",
                    nutDiff + this->nu()
                )
            );
        }


        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence specific dissipation rate
        virtual tmp<volScalarField> omega() const
        {
            return omega_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "epsilon",
                        this->mesh_.time().timeName(),
                        this->mesh_
                    ),
                    Cmu_*k_*omega_,
                    omega_.boundaryField().types()
                )
            );
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "kOmegaLKE.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
