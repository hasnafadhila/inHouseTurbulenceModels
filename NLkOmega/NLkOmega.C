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

\*---------------------------------------------------------------------------*/

#include "NLkOmega.H"
#include "fvOptions.H"
#include "bound.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> NLkOmega<BasicTurbulenceModel>::Ret() const
{
    return (k_/(this->nu()*(omega_ + this->omegaMin_)));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> NLkOmega<BasicTurbulenceModel>::f1(const volScalarField& Ret) const
{
    return (1-exp(-pow(Ret,C1_)/C2_));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> NLkOmega<BasicTurbulenceModel>::f2(const volScalarField& Ret) const
{
    return (exp(-pow(Ret,C3_)/C4_));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> NLkOmega<BasicTurbulenceModel>::f3(const volScalarField& Ret) const
{
    return  (1-tanh(pow(Ret,C5_)/C6_));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> NLkOmega<BasicTurbulenceModel>::Cmu
(
		const volSymmTensorField& S,
		const volTensorField& W
) const
{
		volScalarField STilde((1/((omega_ + this->omegaMin_)))*sqrt(2.0)*mag(S));
		volScalarField OmegaTilde((1/((omega_ + this->omegaMin_)))*sqrt(2.0)*mag(W));
		volScalarField M(max(STilde,OmegaTilde));

		return (min(1.00,(1/(1+(0.01*sqr(M))))));

}

template<class BasicTurbulenceModel>
tmp<volScalarField> NLkOmega<BasicTurbulenceModel>::scalingTerm() const
{
		const volScalarField S2(2.0*magSqr(dev(symm(fvc::grad(this->U_)))));

		const volScalarField scalingTerm1 (k_/sqr((max(omega_, sqrt(S2)*kappa_) + this->omegaMin_)));

		return (max(scalingTerm1,dimensionedScalar("ZERO", sqr(dimLength), 0.0)));

}


template<class BasicTurbulenceModel>
void NLkOmega<BasicTurbulenceModel>::correctNut()
{
      correctNonlinearStress(fvc::grad(this->U_));
}


template<class BasicTurbulenceModel>
void NLkOmega<BasicTurbulenceModel>::correctNonlinearStress(const volTensorField& gradU)
{

		const volSymmTensorField S(symm(gradU));
    const volTensorField W(skew(gradU));

		const volScalarField CmuEff(Cmu(S,W));

		this->nut_ =  k_/(omega_ + this->omegaMin_);
		this->nut_.correctBoundaryConditions();
		fv::options::New(this->mesh_).correct(this->nut_);


    const volScalarField Cbeta1Eff ( min (
						min ( (CV1_*f1(Ret())*f2(Ret())) , CV1_ )
	      +   min ( (CB1_*f1(Ret())*f3(Ret())) , CB1_ )
				+   min ( (CL1_*(1-f3(Ret()))) , CL1_)
                                  ,  CV1_ )
			);

		const volScalarField Cbeta2Eff ( min (
					  min ( (CV2_*f1(Ret())*f2(Ret())) , CV2_ )
	      +   min ( (CB2_*f1(Ret())*f3(Ret())) , CB2_ )
				+   min ( (CL2_*(1-f3(Ret()))), CL2_ )
                                  ,  CV2_ )
			);

    this->nonlinearStress_ =
    (
    CmuEff * scalingTerm()
         *(
             Cbeta1Eff*dev(innerSqr(S))
           + Cbeta2Eff*twoSymm(S&W)
          )
    );

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
NLkOmega<BasicTurbulenceModel>::NLkOmega
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
    nonlinearEddyViscosity<RASModel<BasicTurbulenceModel>>
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

    betaStar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaStar",
            this->coeffDict_,
            0.09
        )
    ),
		Comega1_
		(
				dimensioned<scalar>::lookupOrAddToDict
				(
						"Comega1",
						this->coeffDict_,
						0.52
				)
		),
    Comega2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Comega2",
            this->coeffDict_,
            0.072
        )
    ),
    sigmaK_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaK",
            this->coeffDict_,
            0.5
        )
    ),
    sigmaOmega_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaOmega",
            this->coeffDict_,
            0.5
        )
    ),
		kappa_
		(
				dimensioned<scalar>::lookupOrAddToDict
				(
						"kappa",
						this->coeffDict_,
						2.50
				)
		),
		CV1_
		(
				dimensioned<scalar>::lookupOrAddToDict
				(
						"CV1",
						this->coeffDict_,
						160.0
				)
		),
		CB1_
		(
				dimensioned<scalar>::lookupOrAddToDict
				(
						"CB1",
						this->coeffDict_,
						25.0
				)
		),
		CL1_
		(
				dimensioned<scalar>::lookupOrAddToDict
				(
						"CL1",
						this->coeffDict_,
						10.2
				)
		),
		CV2_
		(
				dimensioned<scalar>::lookupOrAddToDict
				(
						"CV2",
						this->coeffDict_,
						122.0
				)
		),
		CB2_
		(
				dimensioned<scalar>::lookupOrAddToDict
				(
						"CB2",
						this->coeffDict_,
					  15.0
				)
		),
		CL2_
		(
				dimensioned<scalar>::lookupOrAddToDict
				(
						"CL2",
						this->coeffDict_,
						8.0
				)
		),
		C1_
		(
				dimensioned<scalar>::lookupOrAddToDict
				(
						"C1",
						this->coeffDict_,
						0.92				)
		),
		C2_
		(
				dimensioned<scalar>::lookupOrAddToDict
				(
						"C2",
						this->coeffDict_,
						0.01
				)
		),
		C3_
		(
				dimensioned<scalar>::lookupOrAddToDict
				(
						"C3",
						this->coeffDict_,
						0.40
				)
		),
		C4_
		(
				dimensioned<scalar>::lookupOrAddToDict
				(
						"C4",
						this->coeffDict_,
						0.18
				)
		),
		C5_
		(
				dimensioned<scalar>::lookupOrAddToDict
				(
						"C5",
						this->coeffDict_,
						1.90
				)
		),
		C6_
		(
				dimensioned<scalar>::lookupOrAddToDict
				(
						"C6",
						this->coeffDict_,
						70.00
				)
		),
    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
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
            IOobject::groupName("omega", alphaRhoPhi.group()),
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
            this->mesh_,
	          IOobject::NO_READ,
	          IOobject::AUTO_WRITE
        ),
       betaStar_*k_*omega_
    ),
		y_(wallDist::New(this->mesh_).y())
{
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool NLkOmega<BasicTurbulenceModel>::read()
{
    if (nonlinearEddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        betaStar_.readIfPresent(this->coeffDict());
        Comega1_.readIfPresent(this->coeffDict());
        Comega2_.readIfPresent(this->coeffDict());
        sigmaK_.readIfPresent(this->coeffDict());
        sigmaOmega_.readIfPresent(this->coeffDict());
				kappa_.readIfPresent(this->coeffDict());
				CV1_.readIfPresent(this->coeffDict());
				CB1_.readIfPresent(this->coeffDict());
				CL1_.readIfPresent(this->coeffDict());
				CV2_.readIfPresent(this->coeffDict());
				CB2_.readIfPresent(this->coeffDict());
				CL2_.readIfPresent(this->coeffDict());
				C1_.readIfPresent(this->coeffDict());
				C2_.readIfPresent(this->coeffDict());
				C3_.readIfPresent(this->coeffDict());
				C4_.readIfPresent(this->coeffDict());
				C5_.readIfPresent(this->coeffDict());
				C6_.readIfPresent(this->coeffDict());


        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
void NLkOmega<BasicTurbulenceModel>::correct()
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
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    nonlinearEddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    tmp<volTensorField> tgradU = fvc::grad(U);

		const volTensorField& gradU = tgradU();

		const volSymmTensorField S(symm(gradU));
		const volTensorField W(skew(gradU));

		const volScalarField Omega(sqrt(2.0)*mag(skew(gradU)));
		const volScalarField S2(2.0*magSqr(dev(symm(gradU))));

		volScalarField STilde((1/((omega_ + this->omegaMin_)))*sqrt(2.0)*mag(S));
		volScalarField OmegaTilde((1/((omega_ + this->omegaMin_)))*sqrt(2.0)*mag(W));

    volScalarField G
    (
        this->GName(),
        (nut*(gradU && dev(twoSymm(gradU))))
    );

		const volScalarField Pk( (nut*(gradU && dev(twoSymm(gradU)))) -  (this->nonlinearStress_ && gradU) );


    // Update omega and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();

    // Turbulence specific dissipation rate equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(alpha, rho, omega_)
      + fvm::div(alphaRhoPhi, omega_)
      - fvm::laplacian(alpha*rho*(sigmaOmega_*k_/omega_ + this->nu()), omega_)
     ==
        Comega1_*alpha*rho*G*omega_/k_
      - fvm::SuSp(((2.0/3.0)*Comega1_)*alpha*rho*divU, omega_)
      - fvm::Sp(Comega2_*alpha*rho*omega_, omega_)
      + fvOptions(alpha, rho, omega_)
    );

    omegaEqn.ref().relax();
    fvOptions.constrain(omegaEqn.ref());
    omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
    solve(omegaEqn);
    fvOptions.correct(omega_);
    bound(omega_, this->omegaMin_);


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*(sigmaK_*k_/omega_ + this->nu()), k_)
     ==
				min(alpha*rho*Pk,alpha*rho*20*betaStar_*omega_*k_)
      - fvm::SuSp((2.0/3.0)*alpha*rho*divU, k_)
      - fvm::Sp(betaStar_*alpha*rho*omega_, k_)
      + fvOptions(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);


		// Update total fluctuation kinetic energy dissipation rate
    epsilon_ = betaStar_*k_*omega_;
    bound(epsilon_, this->epsilonMin_);

    correctNonlinearStress(gradU);

    tgradU.clear();


}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
