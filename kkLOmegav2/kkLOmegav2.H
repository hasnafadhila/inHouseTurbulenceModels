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
    Foam::RASModels::kkLOmegav2

Group
    grpRASTurbulence

Description
    Low-Reynolds number k-kl-omega-v2 turbulence model

    This turbulence model is described in:
    \verbatim
        Chitta, V., Dhakal, T. P., & Keith Walters, D. (2015). Sensitization of
        a transition-sensitive linear eddy-viscosity model to rotation and
        curvature effects. Journal of Fluids Engineering, 137(3).
    \endverbatim


    With some functions and model coefficients described in:
    /verbatim
	    Walters, D. K., & Cokljat, D. (2008). A three-equation eddy-viscosity
      model for Reynolds-averaged Navier–Stokes simulations of transitional
      flow. Journal of fluids engineering, 130(12).

	    Dhakal, T. P., & Walters, D. K. (2011). A three-equation variant of the
      SST k-omega model sensitized to rotation and curvature effects. Journal
      of fluids engineering, 133(11).
    /endverbatim

    The default model coefficients are
    \verbatim
        kkLOmegav2Coeffs
        {
            A0             4.04
            AS             2.12
            Anu            6.75
            ABP            0.6
            ANAT           200
            ATS            200
            CBPcrit        1.2
            CNC            0.1
            CNATcrit       1250
            CINT           0.75
            CTScrit        1000
            CRNAT          0.02
            Cl1            3.4e-6
            Cl2            1.0e-10
            CR             0.12
            CalphaTheta    0.035
            CSS            1.5
            CtauL          4360
            Comega1        0.44
            Comega2        0.92
            Comega3        0.3
            ComegaR        1.5
            Clambda        2.495
            Cmustd         0.09
            Prtheta        0.85
            Sigmak         1
            Sigmaomega     1.17

	    CRpsi	   1.8
	    BetaStar	   0.09
	    a0		   1.0
	    a1		   18.57
	    a2		   112.0
	    a3		   331.5
	    a4 		   437.8
	    a5		   147.5
        }
    \endverbatim

SourceFiles
    kkLOmegav2.C

\*---------------------------------------------------------------------------*/

#ifndef kkLOmegav2_H
#define kkLOmegav2_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class kkLOmegav2 Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class kkLOmegav2
:
    public eddyViscosity<RASModel<BasicTurbulenceModel> >
{


protected:

    // Protected data

        // Model coefficients for kkLOmega

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
	dimensionedScalar Cl1_;
	dimensionedScalar Cl2_;
	dimensionedScalar CR_;
	dimensionedScalar CalphaTheta_;
	dimensionedScalar CSS_;
	dimensionedScalar CtauL_;
	dimensionedScalar Comega1_;
	dimensionedScalar Comega2_;
	dimensionedScalar Comega3_;
	dimensionedScalar ComegaR_;
	dimensionedScalar Clambda_;
	dimensionedScalar Cmustd_;
	dimensionedScalar Prtheta_;
	dimensionedScalar Sigmak_;
	dimensionedScalar Sigmaomega_;

	// Model coefficients for the v2 equation
	dimensionedScalar CRpsi_;
	dimensionedScalar BetaStar_;
	dimensionedScalar a0_;
	dimensionedScalar a1_;
	dimensionedScalar a2_;
	dimensionedScalar a3_;
	dimensionedScalar a4_;
	dimensionedScalar a5_;

// 	dimensionedScalar k1_;
// 	dimensionedScalar k2_;
// 	dimensionedScalar k3_;
// 	dimensionedScalar k4_;
// 	dimensionedScalar k5_;
// 	dimensionedScalar k6_;
// 	dimensionedScalar k7_;
// 	dimensionedScalar k8_;

        // Fields

        // Turbulence kinetic energy
	volScalarField kt_;

        // Laminar kinetic energy
 	volScalarField kl_;

	// Turbulence frequency
	volScalarField omega_;

	// Turbulence kinetic energy dissipation rate
	volScalarField epsilon_;

        //- Turbulence stress normal to streamlines
	volScalarField v2_;

// 	//***********************************************************
// 	// Trying original formulation of CmuRC (Dhakal 2011, York 2009, Gatski 1993)
// 	volScalarField CmuRC_;
// 	//************************************************************/

	// New output variables (for monitoring)
	volScalarField xOut_;
	volScalarField etaOut_;
	volScalarField etaEffOut_;
	volScalarField psiOut_;

        // Bounding values
	dimensionedScalar v2Min_;



	//- Wall distance
	//  Note: different to wall distance in parent RASModel
	//  which is for near-wall cells only
	const volScalarField& y_;


    // Protected Member Functions

        virtual void correctNut();

	// Private member functions

	tmp<volScalarField> ReT(const volScalarField& fw) const;

	tmp<volScalarField> fnu(const volScalarField& ReT) const;

	tmp<volScalarField> fINT() const;

	tmp<volScalarField> fSS(const volScalarField& Omega) const;

	tmp<volScalarField> Cmu(const volScalarField& S) const;

	tmp<volScalarField> BetaTS(const volScalarField& ReOmega) const;

	tmp<volScalarField> ftauL
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

	tmp<volScalarField> fomega
		(
		const volScalarField& lambdaEff,
		const volScalarField& lambdaT
		) const;

	tmp<volScalarField> phiBP(const volScalarField& Omega) const;

	tmp<volScalarField> phiNAT
		(
		const volScalarField& ReOmega,
		const volScalarField& fNATcrit
		) const;

	tmp<volScalarField> D(const volScalarField& k) const;


	tmp<volScalarField> omegaStar
		(
		const volScalarField& Dt,
		const volScalarField& Dl) const;

public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("kkLOmegav2");


    // Constructors

        //- Construct from components
        kkLOmegav2
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
    virtual ~kkLOmegav2()
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
	virtual tmp<volScalarField> kt() const
	{
		return kt_;
	}

	//- Return the laminar kinetic energy
	virtual tmp<volScalarField> kl() const
	{
		return kl_;
	}

	//- Return the turbulence frequency
	virtual tmp<volScalarField> omega() const
	{
	   	return omega_;
	}

	//- Return the turbulence kinetic energy dissipation rate
	virtual tmp<volScalarField> epsilon() const
	{
	    	return epsilon_;
	}

	//- Return turbulence stress normal to streamlines
	virtual tmp<volScalarField> v2() const
	{
	    	return v2_;
	}

// 	//***********************************************************
// 	// Trying original formulation of CmuRC (Dhakal 2011, York 2009, Gatski 1993)
// 	virtual tmp<volScalarField> CmuRC() const
// 	{
// 	    	return CmuRC_;
// 	}
// 	//***********************************************************

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

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
#   include "kkLOmegav2.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
