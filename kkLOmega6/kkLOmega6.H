/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::RASModels::kkLOmega6

Group
    grpRASTurbulence

Description
    Low Reynolds-number k-kl-omega turbulence model for incompressible and
    compressible flows with the modification suggested by Lopez and
    Walters (2017)

    This turbulence model is described in:
    \verbatim
        Walters, D. K., & Cokljat, D. (2008).
        A three-equation eddy-viscosity model for Reynolds-averaged
        Navier–Stokes simulations of transitional flow.
        Journal of Fluids Engineering, 130(12), 121401.
    \endverbatim

    however the paper contains several errors which must be corrected for the
    model to operation correctly as explained in

    \verbatim
        Furst, J. (2013).
        Numerical simulation of transitional flows with laminar kinetic energy.
        Engineering MECHANICS, 20(5), 379-388.
    \endverbatim

    All these corrections and updates are included in this implementation.

    The effective length scale for the large scale eddy viscosity is modified
    according to:

    \verbatim
        Lopez, M., & Keith Walters, D. (2017). A recommended correction to the
        kT− kL− ω transition-sensitive eddy-viscosity model. Journal of Fluids
        Engineering, 139(2).
    \endvarbatim

    The default model coefficients are
    \verbatim
        kkLOmegaCoeffs
        {
            A0             4.04
            As             2.12
            Anu             6.75
            Abp            0.6
            Anat           200
            Ats            200
            CbpCrit        1.2
            Cnc            0.1
            CnatCrit       1250
            Cint           0.75
            CtsCrit        1000
            CrNat          0.02
            C11            3.4e-6
            C12            1.0e-10
            CR             0.12
            CalphaTheta    0.035
            Css            1.5
            CTaul          4360
            Comega1        0.44
            Comega2        0.92
            Comega3        0.3
            ComegaR        1.5
            Clambda        2.495
            CmuStd         0.09
            Prtheta        0.85
            sigmaK         1
            sigmaOmega     1.17
        }
    \endverbatim

SourceFiles
    kkLOmega6.C


\*---------------------------------------------------------------------------*/

#ifndef kkLOmega6_H
#define kkLOmega6_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                             Class kkLOmega6 Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class kkLOmega6
:
    public eddyViscosity<RASModel<BasicTurbulenceModel> >
{

protected:

    // Protected data

        // Model coefficients

			dimensionedScalar A0_;
			dimensionedScalar As_;
			dimensionedScalar Anu_;
			dimensionedScalar Abp_;
			dimensionedScalar Anat_;
			dimensionedScalar Ats_;
			dimensionedScalar CbpCrit_;
			dimensionedScalar Cnc_;
			dimensionedScalar CnatCrit_;
			dimensionedScalar Cint_;
			dimensionedScalar CtsCrit_;
			dimensionedScalar CrNat_;
			dimensionedScalar C11_;
			dimensionedScalar C12_;
			dimensionedScalar CR_;
			dimensionedScalar CalphaTheta_;
			dimensionedScalar Css_;
			dimensionedScalar CTaul_;
			dimensionedScalar Comega1_;
			dimensionedScalar Comega2_;
			dimensionedScalar Comega3_;
			dimensionedScalar ComegaR_;
			dimensionedScalar Clambda_;
			dimensionedScalar CmuStd_;
			dimensionedScalar Prtheta_;
			dimensionedScalar sigmaK_;
			dimensionedScalar sigmaOmega_;


        // Fields

            //- Turbulent kinetic energy
            volScalarField kt_;

            //- Laminar kinetic energy
            volScalarField kl_;

            //- Turbulence specific dissipation rate
            volScalarField omega_;

            //- Farfield isotropic dissipation rate
            volScalarField epsilon_;


	    //- NEW OUTPUT VARIABLES
	    volSymmTensorField SOut_;
	    volTensorField OmegaOut_;

	    volScalarField nutSOut_;
	    volScalarField epsilontsOut_;
	    volScalarField ktSOut_;

	    volScalarField PklOut_;
	    volScalarField PktOut_;


	    volScalarField RetOut_;
	    volScalarField fvOut_;


			//- Wall distance
			//  Note: different to wall distance in parent RASModel
			//  which is for near-wall cells only
			const volScalarField& y_;


		// Protected Member Functions

		    virtual void correctNut();

			tmp<volScalarField> lambdaT() const;

			tmp<volScalarField> lambdaEff(const volScalarField& lambdaT) const;

			tmp<volScalarField> Ret(const volScalarField& fw) const;

			tmp<volScalarField> fnu(const volScalarField& Ret) const;

			tmp<volScalarField> fINT() const;

			tmp<volScalarField> fSS(const volScalarField& Omega) const;

			tmp<volScalarField> Cmu(const volScalarField& S) const;

			tmp<volScalarField> BetaTS(const volScalarField& ReOmega) const;

			tmp<volScalarField> fTaul
			(
				const volScalarField& lambdaEff,
				const volScalarField& ktL,
				const volScalarField& Omega
			) const;

			tmp<volScalarField> alphaT
			(
				const volScalarField& lambdaEff,
				const volScalarField& fnu,
				const volScalarField& ktS
			) const;

			tmp<volScalarField> fOmega
			(
				const volScalarField& lambdaEff,
				const volScalarField& lambdaT
			) const;

			tmp<volScalarField> phiBP(const volScalarField& Omega) const;

			tmp<volScalarField> phiNAT
			(
				const volScalarField& ReOmega,
				const volScalarField& fNatCrit
			) const;

			tmp<volScalarField> D(const volScalarField& k) const;

			tmp<volScalarField> fw
			(
				const volScalarField& lambdaEff,
				const volScalarField& lambdaT

			) const;

			tmp<volScalarField> ktS
			(
				const volScalarField& fSS,
				const volScalarField& fw
			) const;

			tmp<volScalarField> nutS
			(
				const volScalarField& fnu,
				const volScalarField& fINT,
				const volScalarField& Cmu,
				const volScalarField& ktS,
				const volScalarField& lambdaEff
			) const;

public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("kkLOmega6");


    // Constructors

        //- Construct from components
        kkLOmega6
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
    virtual ~kkLOmega6()
    {}


    // Member Functions

        //- Read RASProperties dictionary
        virtual bool read();

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff(const volScalarField& alphaT) const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff", alphaT/sigmaK_ + this->nu())
            );
        }

        //- Return the effective diffusivity for omega
        tmp<volScalarField> DomegaEff(const volScalarField& alphaT) const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DomegaEff", alphaT/sigmaOmega_ + this->nu())
            );
        }

        //- Return the laminar kinetic energy
        virtual tmp<volScalarField> kl() const
        {
            return kl_;
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> kt() const
        {
            return kt_;
        }

        //- Return the turbulence specific dissipation rate
        virtual tmp<volScalarField> omega() const
        {
            return omega_;
        }

        //- Return the total fluctuation kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "k",
                        this->mesh_.time().timeName(),
                        this->mesh_
                    ),
                    kt_ + kl_,
                    omega_.boundaryField().types()
                )
            );
        }

        //- Return the total fluctuation kinetic energy dissipation rate
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
#   include "kkLOmega6.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
