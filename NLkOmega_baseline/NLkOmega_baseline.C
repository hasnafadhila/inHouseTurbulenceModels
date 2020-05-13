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

\*---------------------------------------------------------------------------*/

#include "NLkOmega_baseline.H"
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
tmp<volScalarField> NLkOmega_baseline<BasicTurbulenceModel>::Cmu
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
tmp<volScalarField> NLkOmega_baseline<BasicTurbulenceModel>::scalingTerm() const
{
		const volScalarField S2(2.0*magSqr(dev(symm(fvc::grad(this->U_)))));

		const volScalarField scalingTerm1 (k_/sqr((max(omega_, sqrt(S2)*kappa_) + this->omegaMin_)));

		return (max(scalingTerm1,dimensionedScalar("ZERO", sqr(dimLength), 0.0)));

}

template<class BasicTurbulenceModel>
void NLkOmega_baseline<BasicTurbulenceModel>::correctNut()
{
      correctNonlinearStress(fvc::grad(this->U_));
}

template<class BasicTurbulenceModel>
void NLkOmega_baseline<BasicTurbulenceModel>::correctNonlinearStress(const volTensorField& gradU)
{
	const volSymmTensorField S(symm(gradU));
	const volTensorField W(skew(gradU));

	const volScalarField CmuEff(Cmu(S,W));

	this->nut_ =  k_/(omega_ + this->omegaMin_);
	this->nut_.correctBoundaryConditions();
	fv::options::New(this->mesh_).correct(this->nut_);

    this->nonlinearStress_ =
    (
    CmuEff * scalingTerm()
         *(
             Cbeta1_*dev(innerSqr(S))
           + Cbeta2_*twoSymm(S&W)
          )
    );

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
NLkOmega_baseline<BasicTurbulenceModel>::NLkOmega_baseline
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
    Cbeta1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cbeta1",
            this->coeffDict_,
            10.20
        )
    ),
    Cbeta2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cbeta2",
            this->coeffDict_,
            8.00
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
bool NLkOmega_baseline<BasicTurbulenceModel>::read()
{
    if (nonlinearEddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        betaStar_.readIfPresent(this->coeffDict());
        Comega1_.readIfPresent(this->coeffDict());
        Comega2_.readIfPresent(this->coeffDict());
        sigmaK_.readIfPresent(this->coeffDict());
        sigmaOmega_.readIfPresent(this->coeffDict());
				kappa_.readIfPresent(this->coeffDict());
        Cbeta1_.readIfPresent(this->coeffDict());
        Cbeta2_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
void NLkOmega_baseline<BasicTurbulenceModel>::correct()
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
