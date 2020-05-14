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
    Foam::RASModels::NLkOmega_baseline

Group
    grpRASTurbulence

Description
    This is an in-house implementation of the non-linear (NL) k-omega
    turbulence model. The model uses a classic k-omega transport equation with
    a quadratic expression for the Reynolds stress tensor and	constant values
    for the anisotropy expansion coefficients (Cbeta1 and Cbeta2).

    References:
    \verbatim
        Fadhila, H., Medina, H., Aleksandrova, S., & Benjamin, S. (2020).
        A new non-linear RANS model with enhanced near-wall treatment of
        turbulence anisotropy. Applied Mathematical Modelling, 82, 293-313.
    \endverbatim

    The default model coefficients are
    \verbatim
        NLkOmega_baselineCoeffs
				{
				    betaStar        0.09;
				    Comega1         0.52;
				    Comega2         0.072;
				    sigmaK          0.5;
				    sigmaOmega      0.5;
				    kappa           2.5;
				    Cbeta1          10.2;
				    Cbeta2	    		8.0;
				}
    \endverbatim

SourceFiles
    NLkOmega_baseline.C

\*---------------------------------------------------------------------------*/

#ifndef NLkOmega_baseline_H
#define NLkOmega_baseline_H

#include "RASModel.H"
#include "nonlinearEddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class NLkOmega_baseline Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class NLkOmega_baseline
:
    public nonlinearEddyViscosity<RASModel<BasicTurbulenceModel>>
{

protected:

    // Protected data

        // Model coefficients

            dimensionedScalar betaStar_;
            dimensionedScalar Comega1_;
            dimensionedScalar Comega2_;
            dimensionedScalar sigmaK_;
            dimensionedScalar sigmaOmega_;

						dimensionedScalar kappa_;

	    			dimensionedScalar Cbeta1_;
	    			dimensionedScalar Cbeta2_;


        // Fields

            volScalarField k_;
            volScalarField omega_;
						volScalarField epsilon_;

						const volScalarField& y_;

						tmp<volScalarField> Ret() const;



    // Protected Member Functions

        virtual void correctNut();

        virtual void correctNonlinearStress(const volTensorField& gradU);



				tmp<volScalarField> Cmu
				(
						const volSymmTensorField& S,
						const volTensorField& W
				) const;

				tmp<volScalarField> scalingTerm() const;




public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("NLkOmega_baseline");


    // Constructors

        //- Construct from components
        NLkOmega_baseline
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
    virtual ~NLkOmega_baseline()
    {}


    // Member Functions

        //- Read RASProperties dictionary
        virtual bool read();

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DkEff",
                    sigmaK_*k_/omega_ + this->nu()
                )
            );
        }

        //- Return the effective diffusivity for omega
        tmp<volScalarField> DomegaEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DomegaEff",
                    sigmaOmega_*k_/omega_ + this->nu()
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
            return epsilon_;
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "NLkOmega_baseline.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
