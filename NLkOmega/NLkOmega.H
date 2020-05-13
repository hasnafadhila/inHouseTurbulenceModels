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

Class
    Foam::RASModels::NLkOmega

Description
    This is the final version of the non-linear (NL) k-omega turbulence model.
		The model uses a quadratic expression for the Reynolds stress tensor with
		an added function for turbulence anisotropy.

    References:
    \verbatim
     Fadhila, H., et al. "A new non-linear RANS model with enhanced near-wall
		 treatment of turbulence anisotropy." Applied Mathematical Modelling 82
		 (2020): 293-313.
    \endverbatim

    The default model coefficients are
    \verbatim
        NLkOmegaCoeffs
				{
				    betaStar        0.09;
				    Comega1         0.52;
				    Comega2         0.072;
				    sigmaK          0.5;
				    sigmaOmega      0.5;
				    kappa           2.5;
				    CV1             160.0;
						CB1             25.0;
						CL1             10.2;
				    CV2             122.0;
						CB2             15.0;
						CL2             8.0;
						C1              0.92;
						C2              0.01;
						C3              0.40;
						C4              0.18;
						C5              1.90;
						C6              70.00;
				}

    \endverbatim

SourceFiles
    NLkOmega.C

\*---------------------------------------------------------------------------*/

#ifndef NLkOmega_H
#define NLkOmega_H

#include "RASModel.H"
#include "nonlinearEddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class NLkOmega Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class NLkOmega
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

						dimensionedScalar CV1_;
						dimensionedScalar CB1_;
						dimensionedScalar CL1_;
						dimensionedScalar CV2_;
						dimensionedScalar CB2_;
						dimensionedScalar CL2_;

						dimensionedScalar C1_;
						dimensionedScalar C2_;
						dimensionedScalar C3_;
						dimensionedScalar C4_;
						dimensionedScalar C5_;
						dimensionedScalar C6_;



        // Fields

            volScalarField k_;
            volScalarField omega_;
						volScalarField epsilon_;

						const volScalarField& y_;

						tmp<volScalarField> Ret() const;


    // Protected Member Functions

        virtual void correctNut();

        virtual void correctNonlinearStress(const volTensorField& gradU);

				tmp<volScalarField> f1(const volScalarField& Ret) const;

				tmp<volScalarField> f2(const volScalarField& Ret) const;

				tmp<volScalarField> f3(const volScalarField& Ret) const;

        tmp<volScalarField> fSS() const;

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
    TypeName("NLkOmega");


    // Constructors

        //- Construct from components
        NLkOmega
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
    virtual ~NLkOmega()
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
    #include "NLkOmega.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
