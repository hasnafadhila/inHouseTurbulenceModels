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

#include "kOmegav2.H"
#include "bound.H"
#include "wallDist.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

	// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //
	
	template<class BasicTurbulenceModel>
	tmp<volScalarField> kOmegav2<BasicTurbulenceModel>::kOmegav2::F1(const volScalarField& CDkomega) const
	{
	    tmp<volScalarField> CDkomegaPlus = max
	    (
		CDkomega,
		dimensionedScalar("1.0e-10", dimless/sqr(dimTime), 1.0e-10)
	    );

	    tmp<volScalarField> arg1 = min
	    (
		    max
		    (
			(sqrt(v2_)/(omega_*y_)),
			(scalar(500)*this->nu()/(sqr(y_)*omega_))
		    ),
		    ((4*Sigmaomega2_*k_)/(CDkomegaPlus*sqr(y_)))
	     );

	    return tanh(pow4(arg1));
	}
	
	template<class BasicTurbulenceModel>
	tmp<volScalarField> kOmegav2<BasicTurbulenceModel>::ReT(const volScalarField& fw) const
	{
		return ((sqr(fw)* v2_) / this->nu() / (omega_ + this->omegaMin_));
	}

	template<class BasicTurbulenceModel>
	tmp<volScalarField> kOmegav2<BasicTurbulenceModel>::fnu(const volScalarField& ReT) const
	{
		return(1.0 - exp(-sqrt(ReT) / Anu_));
	}

	template<class BasicTurbulenceModel>
	tmp<volScalarField> kOmegav2<BasicTurbulenceModel>::fINT() const
	{
		return
			(
				min
				(
					(v2_ / CINT_/ (k_ + this->kMin_)),
					dimensionedScalar("1.0", dimless, 1.0)
				)
			);
	}

	template<class BasicTurbulenceModel>
	tmp<volScalarField> kOmegav2<BasicTurbulenceModel>::fSS(const volScalarField& Omega) const
	{
		return(exp(-sqr((CSS_*this->nu()*Omega) / (v2_ + this->kMin_))));
	}

	template<class BasicTurbulenceModel>
	tmp<volScalarField> kOmegav2<BasicTurbulenceModel>::Cmu(const volScalarField& S) const
	{
		return(1.0 / (A0_ + (AS_*(S / (omega_ + this->omegaMin_)))));
	}

	template<class BasicTurbulenceModel>
	tmp<volScalarField> kOmegav2<BasicTurbulenceModel>::BetaTS(const volScalarField& ReOmega) const
	{
		return(scalar(1) - exp(-sqr(max((ReOmega - CTScrit_), scalar(0))) / ATS_));
	}

	template<class BasicTurbulenceModel>
	tmp<volScalarField> kOmegav2<BasicTurbulenceModel>::ftauL
	(
		const volScalarField& lambdaEff,
		const volScalarField& v2L,
		const volScalarField& Omega
	) const
	{
		return
			(
				scalar(1)
				- exp
				(
					-CtauL_*v2L
					/
					(
						sqr(lambdaEff*Omega)
						+ dimensionedScalar("ROOTVSMALL", dimLength*dimLength*inv(dimTime)*inv(dimTime), ROOTVSMALL)
					)
				)
			);
	}

	template<class BasicTurbulenceModel>
	tmp<volScalarField> kOmegav2<BasicTurbulenceModel>::alphaT
	(
		const volScalarField& lambdaEff,
		const volScalarField& fnu,
		const volScalarField& v2S
	) const
	{
		return(fnu*BetaStar_*sqrt(v2S)*lambdaEff);
	}

	template<class BasicTurbulenceModel>
	tmp<volScalarField> kOmegav2<BasicTurbulenceModel>::phiBP(const volScalarField& Omega) const
	{
		return
		(
		    //min
		    //(
				max
				(
					(v2_ / this->nu() 
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
	tmp<volScalarField> kOmegav2<BasicTurbulenceModel>::phiNAT
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
	tmp<volScalarField> kOmegav2<BasicTurbulenceModel>::D(const volScalarField& k) const
	{
		return 2*this->nu()*magSqr(fvc::grad(sqrt(k)));
	}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void kOmegav2<BasicTurbulenceModel>::correctNut()
{
    // Currently this function is not implemented due to the complexity of
    // evaluating nut.  Better calculate nut at the end of correct()
    // notImplemented("kOmegav2::correctNut()"); //Needs to be stopped in OF4x it's being called from function validate()
    // Someone is going to have to implement.... volunteers? Dhila? :-)
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kOmegav2<BasicTurbulenceModel>::kOmegav2
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
            3.8
        )
    ),
    ABP_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ABP",
            this->coeffDict_,
            0.2
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
            1.5
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
            1450
        )
    ),
    CINT_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CINT",
	    this->coeffDict_,
            0.95
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
    C11_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C11",
	    this->coeffDict_,
            3.4e-6
        )
    ),
    C12_
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
            0.32
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
            3.0
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
    ComegaR_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ComegaR",
	    this->coeffDict_,
            1.15
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
    BetaStar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "BetaStar",
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
    Sigmaomega2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Sigmaomega2",
	     this->coeffDict_,
            1.856
        )
    ),
    k_
    (
        IOobject
        (
            IOobject::groupName("k", U.group()),
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
        this->mesh_// 
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
    epsilon_
    (
        IOobject
        (
            "epsilon",
            this->runTime_.timeName(),
            this->mesh_
        ),
        k_*omega_ + D(k_)
    ),
    /*
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
    ),*/
    y_(wallDist::New(this->mesh_).y())
    {
	bound(k_, this->kMin_);
	bound(v2_, this->kMin_);
	bound(omega_, this->omegaMin_);
	bound(epsilon_, this->epsilonMin_);

   	if (type == typeName)
	{
	        // Evaluating nut_ is complex so start from the field read from file
		this->nut_.correctBoundaryConditions();
		this->printCoeffs(type); 
	}
    }


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool kOmegav2<BasicTurbulenceModel>::read()
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
		C11_.readIfPresent(this->coeffDict());
		C12_.readIfPresent(this->coeffDict());
		CR_.readIfPresent(this->coeffDict());
		CalphaTheta_.readIfPresent(this->coeffDict());
		CSS_.readIfPresent(this->coeffDict());
		CtauL_.readIfPresent(this->coeffDict());
		Comega1_.readIfPresent(this->coeffDict());
		Comega2_.readIfPresent(this->coeffDict());
		ComegaR_.readIfPresent(this->coeffDict());
		Clambda_.readIfPresent(this->coeffDict());
		BetaStar_.readIfPresent(this->coeffDict());
		Prtheta_.readIfPresent(this->coeffDict());
		Sigmak_.readIfPresent(this->coeffDict());
		Sigmaomega_.readIfPresent(this->coeffDict());
		Sigmaomega2_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
void kOmegav2<BasicTurbulenceModel>::correct()
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
	//volScalarField& nut = this->nut_;

	eddyViscosity<RASModel<BasicTurbulenceModel> >::correct();

	const volTensorField gradU(fvc::grad(U));
	const volScalarField Omega(sqrt(2.0)*mag(skew(gradU)));
	const volScalarField S2(2.0*magSqr(dev(symm(gradU))));
	
	
	const volScalarField lambdaT(sqrt(v2_) / (omega_ + this->omegaMin_));

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

	
	const volScalarField v2S(fSS(Omega)*fw*v2_);
	
	const volScalarField nutS
	(
	    fnu(ReT(fw))
	    *fINT()
	    *Cmu(sqrt(S2))
	    *sqrt(v2S)
	    *lambdaEff
	);

	const volScalarField Pk(this->nut_*S2);

	const volScalarField v2L(v2_ - v2S);

	const volScalarField ReOmega((sqr(y_)*Omega) / this->nu());
	
	const volScalarField deff(lambdaEff/Clambda_);
	
	const volScalarField nutL
	(
	    min
	    (
		(
			((ftauL(lambdaEff, v2L, Omega)
			*C11_
			*Omega
			*sqr(lambdaEff)
			*sqrt(v2L)
			*lambdaEff)
			/ this->nu())
			+
			(BetaTS(ReOmega)
			*C12_
			*(sqr(deff)*Omega/this->nu())
			*sqr(deff)
			*Omega)
			
		)
		,
		((0.5*(k_-v2S)) / (sqrt(S2) + dimensionedScalar("ROOTVSMALL", inv(dimTime), ROOTVSMALL)))
	    )
	);

	const volScalarField Pv2(nutS*S2);
	
	const volScalarField Pomega(Comega1_*omega_*nutS*S2/(v2_+this->kMin_));

	const volScalarField alphaTEff
	(
	    alphaT(lambdaEff, fnu(ReT(fw)), v2S)
	);

	const volScalarField BetaBP(1.0 - exp(-phiBP(Omega) / ABP_));

	const volScalarField RBP
	(
	    (CR_*BetaBP*(k_-v2_)*omega_)
	    / (fw + dimensionedScalar ("SMALL", dimless, ROOTVSMALL))
	);

	const volScalarField fNATcrit(1.0 - exp(-CNC_*sqrt(k_)*y_ / this->nu()));

	const volScalarField BetaNAT(1.0 - exp(-phiNAT(ReOmega,fNATcrit) / ANAT_));

	const volScalarField RNAT
	(
		CRNAT_*BetaNAT*(k_-v2_)*Omega
	);
	
	const volScalarField Dk(D(k_));

	const volScalarField Dv2(D(v2_));
    
    //--------------------------------------------------------------
    //------------------TRANSPORT EQNS START HERE-------------------
    //--------------------------------------------------------------

    // Update omega and G at the wall
    //omega_.boundaryField().updateCoeffs();
    omega_.boundaryFieldRef().updateCoeffs();

    const volScalarField CDkomega(2*Sigmaomega2_*(fvc::grad(k_) & fvc::grad(omega_))/(omega_+this->omegaMin_));
    
    // Turbulence specific dissipation rate equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(alpha, rho, omega_)
      + fvm::div(alphaRhoPhi, omega_)
      - fvm::laplacian(alpha * rho * DomegaEff(alphaTEff), omega_)
     ==
	fvm::SuSp
	(
	    alpha * rho
	    *Pomega
	    /(omega_+this->omegaMin_)
	    , omega_
	)
      + fvm::SuSp
        (
            (alpha* rho * ((ComegaR_/(fw + dimensionedScalar ("SMALL", dimless, ROOTVSMALL)))-1.0)
	    *(RBP + RNAT))
	    / (v2_ + this->kMin_)
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
      + fvm::SuSp
	(
	alpha * rho
	* BetaStar_
	* (1- F1(CDkomega)) * fSS(Omega)
	* CDkomega
	/ (omega_ + this->omegaMin_)
	, omega_
	)
    );

    //omegaEqn().relax();
    omegaEqn.ref().relax(); // New framework from 4.x
    //fvOptions.constrain(omegaEqn.ref()); // New framework from 4.x
    //omegaEqn().boundaryManipulate(omega_.boundaryField());
    omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef()); // New framework from 4.x
    solve(omegaEqn);
    //fvOptions.correct(omega_); // New framework from 4.x
    bound(omega_, this->omegaMin_);


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha * rho * DkEff(alphaTEff), k_)
     ==
        ( alpha * rho * Pk )
	- fvm::Sp
	  (
	    alpha * rho
	    * min(omega_*k_,omega_*v2_) / (k_ + this->kMin_)
	    , k_
	  )
      - fvm::Sp(alpha * rho * Dk/(k_ + this-> kMin_), k_)
    );

    //kEqn().relax();
    kEqn.ref().relax(); // New framework from 4.x
    //fvOptions.constrain(klEqn.ref()); // New framework from 4.x
    //kEqn().boundaryManipulate(k_.boundaryField());
    kEqn.ref().boundaryManipulate(k_.boundaryFieldRef()); // New framework from 4.x
    solve(kEqn);
    //fvOptions.correct(k_); // New framework from 4.x
    bound(k_, this->kMin_);


    // Turbulence stress normal to streamlines equation
    tmp<fvScalarMatrix> v2Eqn
    (
	fvm::ddt(alpha, rho, v2_)
      + fvm::div(alphaRhoPhi, v2_)
      - fvm::laplacian(alpha * rho * Dv2Eff(alphaTEff), v2_)
     ==
        fvm::Sp(alpha * rho * (Pv2 + RBP + RNAT)/(v2_+ this->kMin_), v2_)
      - fvm::Sp(alpha * rho * (omega_ + (Dv2/(v2_ + this->kMin_))), v2_)
    );

    v2Eqn.ref().relax(); // New framework from 4.x
    //v2Eqn().relax();
    v2Eqn.ref().boundaryManipulate(v2_.boundaryFieldRef()); // New framework from 4.x
    //v2Eqn().boundaryManipulate(v2_.boundaryField());
    solve(v2Eqn);
    bound(v2_, this->kMin_);

    //-----------------------------------------------------------

    // Update total fluctuation kinetic energy dissipation rate
    epsilon_ = k_*omega_ + Dk;
    
    this->nut_ = nutS + nutL;
    this->nut_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
