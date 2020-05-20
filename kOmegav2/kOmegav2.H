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
    Foam::RASModels::kOmegav2

Group
    grpRASTurbulence

Description
    k-omega-v2 model.

    This turbulence model is described in:
    \verbatim
        Lopez, M., & Walters, D. K. (2016). Prediction of transitional and
        fully turbulent flow using an alternative to the laminar kinetic energy
        approach. Journal of Turbulence, 17(3), 253-273.
    \endverbatim

    The default model coefficients are
    \verbatim
        kOmegav2Coeffs
        {
            A0             4.04
            AS             2.12
            Anu            3.8
            ABP            0.2
            ANAT           200
            ATS            200
            CBPcrit        1.5
            CNC            0.1
            CNATcrit       1450
            CINT           0.95
            CTScrit        1000
            CRNAT          0.02
            C11            3.4e-6
            C12            1.0e-10
            CR             0.32
            CalphaTheta    0.035
            CSS            3.0
            CtauL          4360
            Comega1        0.44
            Comega2        0.92
            ComegaR        1.15
            Clambda        2.495
            BetaStar       0.09
            Prtheta        0.85
            Sigmak         1
            Sigmaomega     1.17
            Sigmaomega2	   1.856
        }
    \endverbatim

SourceFiles
    kOmegav2.C


\*---------------------------------------------------------------------------*/

#ifndef kOmegav2_H
#define kOmegav2_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                             Class kOmegav2 Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class kOmegav2
:
    public eddyViscosity<RASModel<BasicTurbulenceModel> >
{

protected:

    // Protected data

        // Model coefficients

			dimensionedScalar A0_;
			dimensionedScalar AS_;
			dimensionedScalar Anu_;
			dimensionedScalar ABP_;
			dimensionedScalar ANAT_;
			dimensionedScalar ATS_;
			dimensionedScalar CBPcrit_;
			dimensionedScalar CNC_;
			dimensionedScalar CNATcrit_;
			dimensionedScalar CINT_;
			dimensionedScalar CTScrit_;
			dimensionedScalar CRNAT_;
			dimensionedScalar C11_;
			dimensionedScalar C12_;
			dimensionedScalar CR_;
			dimensionedScalar CalphaTheta_;
			dimensionedScalar CSS_;
			dimensionedScalar CtauL_;
			dimensionedScalar Comega1_;
			dimensionedScalar Comega2_;
			dimensionedScalar ComegaR_;
			dimensionedScalar Clambda_;
			dimensionedScalar BetaStar_;
			dimensionedScalar Prtheta_;
			dimensionedScalar Sigmak_;
			dimensionedScalar Sigmaomega_;
			dimensionedScalar Sigmaomega2_;


        // Fields

            //- Turbulent kinetic energy
            volScalarField k_;

            //- Turbulence specific dissipation rate
            volScalarField omega_;

	    //- Turbulence stress normal to streamlines
            volScalarField v2_;

	    // Turbulence kinetic energy dissipation rate
	    volScalarField epsilon_;


	//- Wall distance
	//  Note: different to wall distance in parent RASModel
	//  which is for near-wall cells only
	const volScalarField& y_;


    // Protected Member Functions

        virtual void correctNut();

	// Private member functions

	tmp<volScalarField> F1(const volScalarField& CDkomega) const;

	tmp<volScalarField> ReT(const volScalarField& fw) const;

	tmp<volScalarField> fnu(const volScalarField& ReT) const;

	tmp<volScalarField> fINT() const;

	tmp<volScalarField> fSS(const volScalarField& Omega) const;

	tmp<volScalarField> Cmu(const volScalarField& S) const;

	tmp<volScalarField> BetaTS(const volScalarField& ReOmega) const;

	tmp<volScalarField> ftauL
		(
		const volScalarField& lambdaEff,
		const volScalarField& v2L,
		const volScalarField& Omega
		) const;

	tmp<volScalarField> alphaT
		(
		const volScalarField& lambdaEff,
		const volScalarField& fnu,
		const volScalarField& v2S
		) const;

	tmp<volScalarField> phiBP(const volScalarField& Omega) const;

	tmp<volScalarField> phiNAT
		(
		const volScalarField& ReOmega,
		const volScalarField& fNATcrit
		) const;

	tmp<volScalarField> D(const volScalarField& k) const;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("kOmegav2");


    // Constructors

        //- Construct from components
        kOmegav2
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
    virtual ~kOmegav2()
    {}


    // Member Functions

        //- Read RASProperties dictionary
        virtual bool read();

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff(const volScalarField& alphaT) const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DkEff",
		    (alphaT/Sigmak_) + this->nu()
                )
            );
        }

        //- Return the effective diffusivity for omega
        tmp<volScalarField> DomegaEff(const volScalarField& alphaT) const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DomegaEff",
		    (alphaT/Sigmaomega_) + this->nu()
                )
            );
        }

        //- Return the effective diffusivity for v2
        tmp<volScalarField> Dv2Eff(const volScalarField& alphaT) const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "Dv2Eff",
		    (alphaT/Sigmak_) + this->nu()
                )
            );
        }

	//- Return the turbulence kinetic energy
	virtual tmp<volScalarField> k() const
	{
		return k_;
	}

	//- Return the turbulence frequency
	virtual tmp<volScalarField> omega() const
	{
	   	return omega_;
	}

	//- Return turbulence stress normal to streamlines
	virtual tmp<volScalarField> v2() const
	{
	    	return v2_;
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
#   include "kOmegav2.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
