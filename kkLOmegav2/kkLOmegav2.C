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

\*---------------------------------------------------------------------------*/

#include "kkLOmegav2.H"
#include "bound.H"
#include "wallDist.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

	// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //
	
	template<class BasicTurbulenceModel>
	tmp<volScalarField> kkLOmegav2<BasicTurbulenceModel>::ReT(const volScalarField& fw) const
	{
		return ((sqr(fw)* kt_) / this->nu() / (omega_ + this->omegaMin_));
	}

	template<class BasicTurbulenceModel>
	tmp<volScalarField> kkLOmegav2<BasicTurbulenceModel>::fnu(const volScalarField& ReT) const
	{
		return(1.0 - exp(-sqrt(ReT) / Anu_));
	}

	template<class BasicTurbulenceModel>
	tmp<volScalarField> kkLOmegav2<BasicTurbulenceModel>::fINT() const
	{
		return
			(
				min
				(
					(kt_ / CINT_/ (kl_ + kt_ + this->kMin_)),
					dimensionedScalar("1.0", dimless, 1.0)
				)
			);
	}

	template<class BasicTurbulenceModel>
	tmp<volScalarField> kkLOmegav2<BasicTurbulenceModel>::fSS(const volScalarField& Omega) const
	{
		return(exp(-sqr((CSS_*this->nu()*Omega) / (kt_ + this->kMin_))));
	}

	template<class BasicTurbulenceModel>
	tmp<volScalarField> kkLOmegav2<BasicTurbulenceModel>::Cmu(const volScalarField& S) const
	{
		return(1.0 / (A0_ + (AS_*(S / (omega_ + this->omegaMin_)))));
	}

	template<class BasicTurbulenceModel>
	tmp<volScalarField> kkLOmegav2<BasicTurbulenceModel>::BetaTS(const volScalarField& ReOmega) const
	{
		return(scalar(1) - exp(-sqr(max((ReOmega - CTScrit_), scalar(0))) / ATS_));
	}

	template<class BasicTurbulenceModel>
	tmp<volScalarField> kkLOmegav2<BasicTurbulenceModel>::ftauL
	(
		const volScalarField& lambdaEff,
		const volScalarField& ktL,
		const volScalarField& Omega
	) const
	{
		return
			(
				scalar(1)
				- exp
				(
					-CtauL_*ktL
					/
					(
						sqr(lambdaEff*Omega)
						+ dimensionedScalar("ROOTVSMALL", dimLength*dimLength*inv(dimTime)*inv(dimTime), ROOTVSMALL)
					)
				)
			);
	}

	template<class BasicTurbulenceModel>
	tmp<volScalarField> kkLOmegav2<BasicTurbulenceModel>::alphaT
	(
		const volScalarField& lambdaEff,
		const volScalarField& fnu,
		const volScalarField& ktS
	) const
	{
		return(fnu*Cmustd_*sqrt(ktS)*lambdaEff);
	}

	template<class BasicTurbulenceModel>
	tmp<volScalarField> kkLOmegav2<BasicTurbulenceModel>::fomega
	(
		const volScalarField& lambdaEff,
		const volScalarField& lambdaT
	) const
	{
		return
			(
				scalar(1)
				- exp
				(
					-0.41
					*pow4
					(
						lambdaEff
						/ (
						lambdaT + dimensionedScalar("SMALL", dimLength, SMALL)
						)
					)
				)
			);
	}

	template<class BasicTurbulenceModel>
	tmp<volScalarField> kkLOmegav2<BasicTurbulenceModel>::phiBP(const volScalarField& Omega) const
	{
		return
		(
		    //min
		    //(
				max
				(
					(kt_ / this->nu() 
					/ (Omega + dimensionedScalar("ROOTVSMALL", inv(dimTime), ROOTVSMALL))
					)
					- CBPcrit_,
					scalar(0)
				)
				//),
			//scalar(50.0)
		    //)
		);
	}

	template<class BasicTurbulenceModel>
	tmp<volScalarField> kkLOmegav2<BasicTurbulenceModel>::phiNAT
	(
		const volScalarField& ReOmega,
		const volScalarField& fNATcrit
	) const
	{
		return
			(
				max
				(
					(ReOmega
					- (CNATcrit_/ (fNATcrit + dimensionedScalar("ROTVSMALL", dimless, ROOTVSMALL)))),
					scalar(0)
				)
			);
	}

	template<class BasicTurbulenceModel>
	tmp<volScalarField> kkLOmegav2<BasicTurbulenceModel>::D(const volScalarField& k) const
	{
		return this->nu()*magSqr(fvc::grad(sqrt(k)));
	}

	// For calculation of blending function FW (refer to Chitta (2015))
	template<class BasicTurbulenceModel>
	tmp<volScalarField> kkLOmegav2<BasicTurbulenceModel>::omegaStar
	(
		const volScalarField& Dt,
		const volScalarField& Dl
	) const
	{
		return
			(
				(
					omega_
					+ ((Dt + Dl) / (kt_ + kl_ + this->kMin_))
				)
				/ BetaStar_
			);
	}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void kkLOmegav2<BasicTurbulenceModel>::correctNut()
{
    //notImplemented("kkLOmegav2::correctNut()");
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kkLOmegav2<BasicTurbulenceModel>::kkLOmegav2
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<RASModel<BasicTurbulenceModel> >
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),
    
    //Model coefficients for kkLOmega
    A0_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "A0",
            this->coeffDict_,
            4.04
        )
    ),
    AS_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "AS",
            this->coeffDict_,
            2.12
        )
    ),
    Anu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Anu",
            this->coeffDict_,
            6.75
        )
    ),
    ABP_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ABP",
            this->coeffDict_,
            0.6
        )
    ),
    ANAT_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ANAT",
	    this->coeffDict_,
            200
        )
    ),
    ATS_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ATS",
	    this->coeffDict_,
            200
        )
    ),
    CBPcrit_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CBPcrit",
	    this->coeffDict_,
            1.2
        )
    ),
    CNC_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CNC",// 
	    this->coeffDict_,
            0.1
        )
    ),
    CNATcrit_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CNATcrit",
	    this->coeffDict_,
            1250
        )
    ),
    CINT_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CINT",
	    this->coeffDict_,
            0.75
        )
    ),
    CTScrit_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CTScrit",
	    this->coeffDict_,
            1000
        )
    ),
    CRNAT_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CRNAT",
	    this->coeffDict_,
            0.02
        )
    ),
    Cl1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cl1",
	    this->coeffDict_,
            3.4e-6
        )
    ),
    Cl2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cl2",
	    this->coeffDict_,
            1.0e-10
        )
    ),
    CR_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CR",
	    this->coeffDict_,
            0.12
        )
    ),
    CalphaTheta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CalphaTheta",
	    this->coeffDict_,
            0.035
        )
    ),
    CSS_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CSS",
	    this->coeffDict_,
            1.5
        )
    ),
    CtauL_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CtauL",
	    this->coeffDict_,
            4360
        )
    ),
    Comega1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Comega1",
	    this->coeffDict_,
            0.44
        )
    ),
    Comega2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Comega2",
	    this->coeffDict_,
            0.92
        )
    ),
    Comega3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Comega3",
	    this->coeffDict_,
            0.3
        )
    ),
    ComegaR_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ComegaR",
	    this->coeffDict_,
            1.5
        )
    ),
    Clambda_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Clambda",
	    this->coeffDict_,
            2.495
        )
    ),
    Cmustd_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmustd",
	    this->coeffDict_,
            0.09
        )
    ),
    Prtheta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Prtheta",
	    this->coeffDict_,
            0.85
        )
    ),
    Sigmak_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Sigmak",
	    this->coeffDict_,
            1
        )
    ),
    Sigmaomega_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Sigmaomega",
	     this->coeffDict_,
            1.17
        )
    ),
    // Model coefficients for the v2 equation
    CRpsi_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CRpsi",
	    this->coeffDict_,
            1.8
	)
    ),
    BetaStar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "BetaStar",
	    this->coeffDict_,
            0.09
	)
    ),
    a0_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "a0",
	    this->coeffDict_,
            1.0
        )
    ),
    a1_
    (
	dimensioned<scalar>::lookupOrAddToDict
	(
	    "a1",
	    this->coeffDict_,
	    18.57
	)
    ),
    a2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "a2",
	    this->coeffDict_,
            112.0
        )
    ),
    a3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "a3",
	    this->coeffDict_,
            331.5
        )
    ),
    a4_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "a4",
	    this->coeffDict_,
            437.8
        )
    ),
    a5_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "a5",
	    this->coeffDict_,
            147.5
        )
    ),
//     k1_
//     (
//         dimensioned<scalar>::lookupOrAddToDict
//         (
//             "k1",
// 	    this->coeffDict_,
//             1.12
//         )
//     ),
//     k2_
//     (
//         dimensioned<scalar>::lookupOrAddToDict
//         (
//             "k2",
// 	    this->coeffDict_,
//             3.32125
//         )
//     ),
//     k3_
//     (
//         dimensioned<scalar>::lookupOrAddToDict
//         (
//             "k3",
// 	    this->coeffDict_,
//             -1.9304
//         )
//     ),
//     k4_
//     (
//         dimensioned<scalar>::lookupOrAddToDict
//         (
//             "k4",
// 	    this->coeffDict_,
//             -5.2397
//         )
//     ),
//     k5_
//     (
//         dimensioned<scalar>::lookupOrAddToDict
//         (
//             "k5",
// 	    this->coeffDict_,
//             2.94
//         )
//     ),
//     k6_
//     (
//         dimensioned<scalar>::lookupOrAddToDict
//         (
//             "k6",
// 	    this->coeffDict_,
//             15.96
//         )
//     ),
//     k7_
//     (
//         dimensioned<scalar>::lookupOrAddToDict
//         (
//             "k7",
// 	    this->coeffDict_,
//             21.66
//         )
//     ),
//     k8_
//     (
//         dimensioned<scalar>::lookupOrAddToDict
//         (
//             "k8",
// 	    this->coeffDict_,
//             1.28
//         )
//     ),
    kt_
    (
        IOobject
        (
            IOobject::groupName("kt", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    kl_
    (
        IOobject
        (
            IOobject::groupName("kl", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    epsilon_
    (
        IOobject
        (
            "epsilon",
            this->runTime_.timeName(),
            this->mesh_
        ),
        kt_*omega_ + D(kl_)
	+ D(kt_)
    ),
    v2_
    (
        IOobject
        (
            IOobject::groupName("v2", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
//     //***********************************************************	
//     // Trying original formulation of CmuRC (Dhakal 2011, York 2009, Gatski 1993)
//     CmuRC_
//     (
//         IOobject
//         (
//             IOobject::groupName("CmuRC", U.group()),
//             this->runTime_.timeName(),
//             this->mesh_,
//             IOobject::MUST_READ,
//             IOobject::AUTO_WRITE
//         ),
//         this->mesh_
//     ),
//     //***********************************************************    
    // New output variables (for monitoring)
    xOut_
    (
      IOobject
      (
	  "x",
	  this->runTime_.timeName(),
	  this->mesh_,
	  IOobject::NO_READ,
	  IOobject::AUTO_WRITE
      ),
      this->mesh_,
      dimensionedScalar("xOut", dimless, 0.0)
    ),
    etaOut_
    (
      IOobject
      (
	    "eta",
	    this->runTime_.timeName(),
	    this->mesh_,
	    IOobject::NO_READ,
	    IOobject::AUTO_WRITE
      ),
      this->mesh_,
      dimensionedScalar("etaOut", dimless, 0.0)
    ),
    etaEffOut_
    (
      IOobject
      (
	    "etaEff",
	    this->runTime_.timeName(),
	    this->mesh_,
	    IOobject::NO_READ,
	    IOobject::AUTO_WRITE
      ),
      this->mesh_,
      dimensionedScalar("etaEffOut", dimless, 0.0)
    ),
    psiOut_
    (
      IOobject
      (
	    "psi",
	    this->runTime_.timeName(),
	    this->mesh_,
	    IOobject::NO_READ,
	    IOobject::AUTO_WRITE
      ),
      this->mesh_,
      dimensionedScalar("psiOut", inv(dimTime), 0.0)
    ),
    v2Min_(dimensionedScalar("v2Min", v2_.dimensions(), SMALL)),
    y_(wallDist::New(this->mesh_).y())
    {
	bound(kt_, this->kMin_);
	bound(kl_, this->kMin_);
	bound(v2_, v2Min_);
	bound(omega_, this->omegaMin_);
	//bound(epsilon_, this->epsilonMin_);

   	if (type == typeName)
	{
		this->nut_.correctBoundaryConditions();
		this->printCoeffs(type); 
	}
    }


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool kkLOmegav2<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel> >::read())
    {
		A0_.readIfPresent(this->coeffDict());
		AS_.readIfPresent(this->coeffDict());
		Anu_.readIfPresent(this->coeffDict());
		ABP_.readIfPresent(this->coeffDict());
		ANAT_.readIfPresent(this->coeffDict());
		ATS_.readIfPresent(this->coeffDict());
		CBPcrit_.readIfPresent(this->coeffDict());
		CNC_.readIfPresent(this->coeffDict());
		CNATcrit_.readIfPresent(this->coeffDict());
		CINT_.readIfPresent(this->coeffDict());
		CTScrit_.readIfPresent(this->coeffDict());
		CRNAT_.readIfPresent(this->coeffDict());
		Cl1_.readIfPresent(this->coeffDict());
		Cl2_.readIfPresent(this->coeffDict());
		CR_.readIfPresent(this->coeffDict());
		CalphaTheta_.readIfPresent(this->coeffDict());
		CSS_.readIfPresent(this->coeffDict());
		CtauL_.readIfPresent(this->coeffDict());
		Comega1_.readIfPresent(this->coeffDict());
		Comega2_.readIfPresent(this->coeffDict());
		Comega3_.readIfPresent(this->coeffDict());
		ComegaR_.readIfPresent(this->coeffDict());
		Clambda_.readIfPresent(this->coeffDict());
		Cmustd_.readIfPresent(this->coeffDict());
		Prtheta_.readIfPresent(this->coeffDict());
		Sigmak_.readIfPresent(this->coeffDict());
		Sigmaomega_.readIfPresent(this->coeffDict());

		CRpsi_.readIfPresent(this->coeffDict());
		BetaStar_.readIfPresent(this->coeffDict());
		a0_.readIfPresent(this->coeffDict());
		a1_.readIfPresent(this->coeffDict());
		a2_.readIfPresent(this->coeffDict());
		a3_.readIfPresent(this->coeffDict());
		a4_.readIfPresent(this->coeffDict());
		a5_.readIfPresent(this->coeffDict());
		
// 		k1_.readIfPresent(this->coeffDict());
// 		k2_.readIfPresent(this->coeffDict());
// 		k3_.readIfPresent(this->coeffDict());
// 		k4_.readIfPresent(this->coeffDict());
// 		k5_.readIfPresent(this->coeffDict());
// 		k6_.readIfPresent(this->coeffDict());
// 		k7_.readIfPresent(this->coeffDict());
// 		k8_.readIfPresent(this->coeffDict());


        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
void kkLOmegav2<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

	// Local references
	const alphaField& alpha = this->alpha_;
	const rhoField& rho = this->rho_;
	const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
	const volVectorField& U = this->U_;

	eddyViscosity<RASModel<BasicTurbulenceModel> >::correct();

	const volTensorField gradU(fvc::grad(U));
	const volScalarField Omega(sqrt(2.0)*mag(skew(gradU)));
	const volScalarField S2(2.0*magSqr(dev(symm(gradU))));
	
	// Modified/effective rotation-rate magnitude (refer to York (2009)) - for calc of x
	// To avoid the use of rotating frame
	const volScalarField W(mag(((9.0/4.0)*Omega)-((5.0/4.0)*sqrt(S2))));
	
	const volScalarField lambdaT(sqrt(kt_) / (omega_ + this->omegaMin_));

	const volScalarField lambdaEff(min((Clambda_*y_), lambdaT));

	const volScalarField fw
	(
	    pow
	    (
			(lambdaEff
			/(lambdaT + dimensionedScalar("ROOTVSMALL", dimLength, ROOTVSMALL))),
			(2.0 / 3.0)
	    )
	);

	const volScalarField ktS(fSS(Omega)*fw*v2_);

	const volScalarField nutS
	(
	    fnu(ReT(fw))
	    *fINT()
	    *Cmu(sqrt(S2))
	    *sqrt(ktS)
	    *lambdaEff
	);

	const volScalarField PkT(nutS*S2);

	const volScalarField ktL(v2_ - ktS);

	const volScalarField ReOmega((sqr(y_)*Omega) / this->nu());

	const volScalarField nutL
	(
	    min
	    (
		(
			((ftauL(lambdaEff, ktL, Omega)
			*Cl1_
			*Omega
			*sqr(lambdaEff)
			*sqrt(ktL)
			*lambdaEff)
			/ this->nu())
			+
			(sqrt(v2_/(kt_ + dimensionedScalar("ROOTVSMALL", dimLength*dimLength*inv(dimTime)*inv(dimTime), ROOTVSMALL)))
			*BetaTS(ReOmega)
			*Cl2_
			*ReOmega
			*Omega
			*sqr(y_))
		)
		,
		((0.5*(kl_ + ktL)) / (sqrt(S2) + dimensionedScalar("ROOTVSMALL", inv(dimTime), ROOTVSMALL)))
	    )
	);


	const volScalarField PkL(nutL*S2);

	const volScalarField alphaTEff
	(
	    alphaT(lambdaEff, fnu(ReT(fw)), ktS)
	);

	const volScalarField BetaBP(1.0 - exp(-phiBP(Omega) / ABP_));

	const volScalarField RBP
	(
	    (CR_*BetaBP*kl_*omega_)
	    / (fw + dimensionedScalar("SMALL", dimless, ROOTVSMALL))
	);


	const volScalarField fNATcrit(1.0 - exp(-CNC_*sqrt(kl_)*y_ / this->nu()));
	
	const volScalarField BetaNAT(1.0 - exp(-phiNAT(ReOmega,fNATcrit) / ANAT_));

	const volScalarField RNAT
	(
		CRNAT_*BetaNAT*kl_*Omega
	);

	
	const volScalarField Dl(D(kl_));

	const volScalarField Dt(D(kt_));

	const volScalarField Dv2(D(v2_));


	// Original equation for x = omegaStar/S (refer to Dhakal (2011))
	// omegaStar is the reference frame rotation rate
	// will be termed omegam here to avoid confusion with omegaStar used in FW function
	// refer to York (2009) for the formulation
	const volScalarField omegam(0.5*(sqrt(S2)-Omega));

	// Original equation for x is omegam/S (refer to Dhakal (2011))
	// const x(omegam / (sqrt(S2) + this->omegaMin_));

	// Alternative equation for x (refer to Dhakal (2011), from York (2009)), works well
	const volScalarField x((2.0/9.0)*(1-(W/(sqrt(S2) + dimensionedScalar("ROOTVSMALL", inv(dimTime), ROOTVSMALL)))));
	
		
	// Weak equilibrium ratio of rotating to non-rotating eddy viscosity (Cmu rot / Cmu non-rot)
	const volScalarField eta
	(

	  max
	  (
		(
		(a5_ * pow(x, 5.0))
		- (a4_ * pow(x, 4.0))
		+ (a3_ * pow(x, 3.0))
		- (a2_ * pow(x, 2.0))
		+ (a1_ * x)
		+ a0_ ),
	    // Needs to be bound to remain non-negative (refer to Dhakal (2011))
	    scalar(0)
	  )
	);
	
// 	//***********************************************************
// 	// Trying original formulation of CmuRC (Dhakal 2011, York 2009, Gatski 1993)
// 	const volScalarField eta
// 	(
// 	  max
// 	  (
// 	    (CmuRC_/Cmustd_)
// 	  // Needs to be bound to remain non-negative (refer to Dhakal (2011))
// 	  ,scalar(0))
// 	);
// 	//***********************************************************
	
	
	
	// Blending function (1 very close to wall, 0 far from wall) (refer to Chitta (2015))
	const volScalarField FW
	(
	    tanh
		(
		    pow4
		    (
			(scalar(200) * this->nu())
			/ ((omegaStar(Dt,Dl) + this->omegaMin_) * sqr(y_))
		    )
		)
	);

	// Applying near-wall limitation on eta (refer to Chitta (2015))
	const volScalarField etaEff
	(
	    (FW* min( 1.0 , eta))
	    + (( 1.0 - FW) * eta)
	);
	
	// psi in v2 equation (refer to Dhakal (2011) and Chitta (2015))
	const volScalarField psi
	(
	    CRpsi_
	    * BetaStar_
	    * omega_
	);
        
// 	//***********************************************************
// 	// Trying original formulation of CmuRC (Dhakal 2011, York 2009, Gatski 1993)
// 	const volScalarField phiS
// 	(
// 	    sqrt(S2) * (kt_ + kl_) / (epsilon_ + dimensionedScalar("ROOTVSMALL", sqr(dimLength) * pow3(inv(dimTime)), ROOTVSMALL))
// 	);
// 	
// 	const volScalarField phiW
// 	(
// 	    W * (kt_ + kl_) / (epsilon_ + dimensionedScalar("ROOTVSMALL", sqr(dimLength) * pow3(inv(dimTime)), ROOTVSMALL))
// 	);
// 	
// 	//***********************************************************
	
    // Writing out some variables for monitoring 
    xOut_ = x;
    etaOut_ = eta;
    etaEffOut_ = etaEff;
    psiOut_ = psi;
    
    //--------------------------------------------------------------
    //------------------TRANSPORT EQNS START HERE-------------------
    //--------------------------------------------------------------

//     //**************************************************************
//     // Trying original formulation of CmuRC (Dhakal 2011, York 2009, Gatski 1993)
//     tmp<fvScalarMatrix> CmuRCEqn
//     (
//        fvm::SuSp((k1_ + (k2_ * CmuRC_ * sqr(phiS)) + (k3_ * CmuRC_ * phiS) + (k4_ * sqr(CmuRC_) * pow3(phiS))) /CmuRC_,CmuRC_)
//        ==
//        fvm::SuSp(k5_ + (k6_ * CmuRC_ * sqr(phiS)) + (k7_ * sqr(CmuRC_) * pow4(phiS)) + (k8_ * sqr(phiW)), CmuRC_)
// 	 
//     );
//     //**************************************************************
       
//     //CmuRCEqn.ref().relax();
//     //CmuRCEqn.ref().boundaryManipulate(CmuRC_.boundaryFieldRef());
//     solve(CmuRCEqn);
//     //bound(CmuRC_, CmuRCMin_);
    
    // Update omega at the wall
    omega_.boundaryFieldRef().updateCoeffs(); // New framework from 4.x
    //omega_.boundaryField().updateCoeffs();

    // Turbulence specific dissipation rate equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(alpha, rho, omega_)
      + fvm::div(alphaRhoPhi, omega_)
      - fvm::laplacian(alpha * rho * DomegaEff(alphaTEff), omega_)
     ==
        (
	    (alpha * rho
	    *Comega1_
	    *PkT
	    *omega_
	    *sqrt(kt_ / (v2_ + v2Min_))
	    )
	    /(kt_ + this->kMin_)
	)
      - fvm::SuSp
        (
            (alpha* rho * (1.0 - (ComegaR_/(fw + dimensionedScalar("SMALL", dimless, ROOTVSMALL))))
	    *(RBP + RNAT))
	    / (kt_ + this->kMin_)
	    , omega_
        )
      - fvm::Sp
	(
	    alpha * rho
	    *Comega2_
	    *sqr(fw)
	    *omega_
	    , omega_
	)
      + (
            alpha * rho
            *Comega3_
	    *fomega(lambdaEff, lambdaT)
	    *alphaTEff
	    *sqr(fw)
	    *sqrt(kt_)
	    /pow3(y_)
	)
    );

    omegaEqn.ref().relax(); // New framework from 4.x
    //omegaEqn().relax();
    omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef()); // New framework from 4.x
    //omegaEqn().boundaryManipulate(omega_.boundaryField());
    solve(omegaEqn);
    bound(omega_, this->omegaMin_);


    // Laminar kinetic energy equation
    tmp<fvScalarMatrix> klEqn
    (
        fvm::ddt(alpha, rho, kl_)
      + fvm::div(alphaRhoPhi, kl_)
      - fvm::laplacian(alpha * rho * this->nu(), kl_)
     ==
        alpha * rho * PkL
      - fvm::Sp(alpha * rho * ((RBP + RNAT + Dl)/(kl_ + this->kMin_)), kl_)
    );
    
    klEqn.ref().relax();
    klEqn.ref().boundaryManipulate(kl_.boundaryFieldRef());
    solve(klEqn);
    bound(kl_, this->kMin_);


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> ktEqn
    (
        fvm::ddt(alpha, rho, kt_)
      + fvm::div(alphaRhoPhi, kt_)
      - fvm::laplacian(alpha * rho * DkEff(alphaTEff), kt_)
     ==
        alpha * rho * (PkT + RBP + RNAT)
      - fvm::Sp(alpha * rho * (omega_ + (Dt/(kt_ + this-> kMin_))), kt_)
    );

    ktEqn.ref().relax();
    ktEqn.ref().boundaryManipulate(kt_.boundaryFieldRef());
    solve(ktEqn);
    bound(kt_, this->kMin_);


    // Turbulence stress normal to streamlines equation
    tmp<fvScalarMatrix> v2Eqn
    (
	fvm::ddt(alpha, rho, v2_)
      + fvm::div(alphaRhoPhi, v2_)
      - fvm::laplacian(alpha * rho * Dv2Eff(alphaTEff), v2_)
     ==
        fvm::Sp(alpha * rho * ((PkT + RBP + RNAT)/(kt_+ this->kMin_)), v2_)
      - fvm::Sp(alpha * rho * (omega_ + (Dv2/(v2_ + v2Min_))), v2_)
    // ************************** TEST ****************************  
      + fvm::Sp(alpha * rho * psi * (((sqr(etaEff) * kt_ )/(v2_+ v2Min_))-1), v2_)
    // *************************************************************
    );

    v2Eqn.ref().relax();
    v2Eqn.ref().boundaryManipulate(v2_.boundaryFieldRef());
    solve(v2Eqn);
    bound(v2_, v2Min_);

    //-----------------------------------------------------------

    // Update total fluctuation kinetic energy dissipation rate
    epsilon_ = (kt_*omega_) + Dl + Dt;
    bound(epsilon_, this->epsilonMin_);
    
    this->nut_ = nutS + nutL;
    this->nut_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
