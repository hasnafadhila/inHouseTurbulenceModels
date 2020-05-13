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

#include "kOmegaLKE.H"
#include "fvOptions.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void kOmegaLKE<BasicTurbulenceModel>::correctNut
(
    const volScalarField& S2
    //const volScalarField& Omega
)
{
    /*
    this->nut_ = to be implemented
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
    */
}

template<class BasicTurbulenceModel>
void kOmegaLKE<BasicTurbulenceModel>::correctNut()
{
    /*
     * TO BE IMPLEMENTED *
    tmp<volTensorField> tgradU = fvc::grad(this->U_);
    volScalarField Omega(sqrt(2.0)*mag(skew(tgradU())));
    volScalarField S2(2*magSqr(dev(symm(tgradU()))));
    tgradU.clear();
    //correctNut(S2, Omega);
    correctNut(S2);
    */
}

// * * * * * * * * * * * * New Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaLKE<BasicTurbulenceModel>::ReLambda() const
{
    return max(mag(this->U_), sqrt(this->kMin_))*y_/this->nu();
} 

template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaLKE<BasicTurbulenceModel>::ReUpsilon() const
{
    return pow(2.0*this->nu()*kl_/sqr(y_)*this->nu(), 0.25)*y_/this->nu();
} 

template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaLKE<BasicTurbulenceModel>::fSS(const volScalarField& Omega) const
{
    dimensionedScalar yMin_ ("yMin_", dimLength, ROOTVSMALL);
    return exp(-pow(cSS_*this->nu()*Omega/max(k_, this->kMin_), 2));
} 

template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaLKE<BasicTurbulenceModel>::fv() const
{
    const volScalarField ReT (k_/(this->nu()*max(omega_, this->omegaMin_))); 
    //const volScalarField fv (1 - exp(-pow(ReT, 0.5)/cV_));
    return (1 - exp(-pow(ReT, 0.5)/cV_));
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kOmegaLKE<BasicTurbulenceModel>::kOmegaLKE
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
    eddyViscosity<RASModel<BasicTurbulenceModel>>
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

    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.09
        )
    ),
    cOmega1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cOmega1_",
            this->coeffDict_,
            0.52
        )
    ),
    cOmega2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cOmega2_",
            this->coeffDict_,
            0.0708
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
    sigmaKl_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaKl",
            this->coeffDict_,
            0.0125
        )
    ),
    tuInf_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "tuInf",
            this->coeffDict_,
            0
        )
    ),
    cCrit_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cCrit",
            this->coeffDict_,
            76500
        )
    ),
    cSS_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cSS",
            this->coeffDict_,
            1.45
        )
    ),
    cV_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cV",
            this->coeffDict_,
            0.43
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
    
    y_(wallDist::New(this->mesh_).y())
    
{
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
    
    if (tuInf_.value() == 0)
    {
        FatalErrorInFunction
                        << "A value for the freestream turbulence level (tuInf) has not been provided" << endl
                        << "tuInf must be defined in the coefficients dictionary:" << endl
                        << this->coeffDict_.name() << endl
                        << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool kOmegaLKE<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        cOmega1_.readIfPresent(this->coeffDict());
        cOmega2_.readIfPresent(this->coeffDict());
        sigmaK_.readIfPresent(this->coeffDict());
        sigmaOmega_.readIfPresent(this->coeffDict());
        sigmaKl_.readIfPresent(this->coeffDict());
        tuInf_.readIfPresent(this->coeffDict());
        cCrit_.readIfPresent(this->coeffDict());
        cSS_.readIfPresent(this->coeffDict());
        cV_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
void kOmegaLKE<BasicTurbulenceModel>::correct()
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
    dimensionedScalar& omegaMin_ = this->omegaMin_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    
    tmp<volTensorField> tgradU = fvc::grad(U);
    const volTensorField& gradU = tgradU();
    const volScalarField Omega(sqrt(2.0)*mag(skew(gradU)));
    const volScalarField S2(2.0*magSqr(dev(symm(gradU))));
    
    volScalarField G
    (
        this->GName(),
        nut*(tgradU() && dev(twoSymm(tgradU())))
    );
    
    
    //const volScalarField Pk( k_/omega_*(tgradU() && dev(twoSymm(tgradU()))) );
    const volScalarField Pk( (tgradU() && dev(twoSymm(tgradU()))) );
    
    ///////////////// TRANSITION FUNCTIONS AND MODIFICATIONS //////////////////////////////
    
    const dimensionedScalar eta_ (0.02974*tanh(59.79*pow(tuInf_, 1.191) + 1.65e-13));
    
    const volScalarField nutL1
    (
        eta_*kl_*sqrt(S2)*pow( ReUpsilon() , -1.30)*pow( ReLambda() , 0.5)/max(S2, sqr(mag(U)/y_))
    );
    
    const volScalarField ReL 
    (
        kl_/min(this->nu(), nutL1)/max(Omega, omegaMin_)
    );
    
    const volScalarField gamma 
    (
        min (pow( ReL, 2), cCrit_ )/cCrit_
    ); 
    

    //const volScalarField ReT (k_/(this->nu()*max(omega_, omegaMin_))); 
    //const volScalarField fv (1 - exp(-pow(ReT, 0.5)/cV_));
    
    ////////////////// Laminar model support functions ///////////////////////////////////
    
    const volScalarField gammaL
    (
        sqrt(kl_)*y_
    );
    
    const volScalarField Pkl
    (      
        eta_*kl_*sqrt(S2)*pow( ReUpsilon() , -1.30)*pow( ReLambda() , 0.5)
    );
    
    
    ///////////////// OUTPUT VARIABLES ///////////////////////////////////////////////////
    
    // PLACE HOLDER FOR OUTPUT VARIABLES TO BE ADDED IF NEEDED
    
    ///////////////////////////////////////////////////////////////////////////////////////
    
    // Laminar kinetic energy equation
    tmp<fvScalarMatrix> klEqn
    (
        fvm::ddt(alpha, rho, kl_)
      + fvm::div(alphaRhoPhi, kl_)
      - fvm::laplacian(alpha*rho*( this->nu() + sigmaKl_*gammaL ), kl_)
     ==
        alpha*rho*Pkl
      - fvm::SuSp((2.0/3.0)*alpha*rho*divU, kl_)
      - fvm::Sp(alpha*rho*2*this->nu()/sqr(y_), kl_)
      + fvOptions(alpha, rho, kl_)
    );

    klEqn.ref().relax();
    fvOptions.constrain(klEqn.ref());
    solve(klEqn);
    fvOptions.correct(kl_);
    bound(kl_, this->kMin_);
    
    // Update omega and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();

    // Turbulence specific dissipation rate equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(alpha, rho, omega_)
      + fvm::div(alphaRhoPhi, omega_)
      - fvm::laplacian(alpha*rho*DomegaEff(sigmaOmega_*gamma*k_/max(omega_,this->omegaMin_)), omega_)
     ==
        cOmega1_*alpha*rho*k_/max(omega_,this->omegaMin_)*Pk*omega_/max(k_,this->kMin_)
      - fvm::SuSp(((2.0/3.0)*cOmega1_)*alpha*rho*divU, omega_)
      - fvm::Sp(cOmega2_*alpha*rho*omega_, omega_) 
      + fvOptions(alpha, rho, omega_)
      + fvm::Sp
        (
            alpha*rho
            *max(
                    (0.125)*(fvc::grad(k_) & fvc::grad(omega_))
                , 
                    dimensionedScalar("0", dimless/pow3(dimTime), 0)
                )
            /max(omega_,this->omegaMin_)/max(omega_,this->omegaMin_), omega_
        )
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
      - fvm::laplacian(alpha*rho*DkEff(sigmaK_*gamma*k_/max(omega_,this->omegaMin_)), k_)
     ==
        fv()*alpha*rho*k_/max(omega_,this->omegaMin_)*Pk*gamma
      - fvm::SuSp((2.0/3.0)*alpha*rho*divU, k_)
      //- Cmu_*alpha*rho*omega_*k_*gamma
      - fvm::Sp(Cmu_*alpha*rho*gamma*omega_, k_)
      + fvOptions(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    tgradU.clear();
    
    const volScalarField nutL
    (
        eta_*kl_*sqrt(S2)*pow( ReUpsilon() , -1.30)*pow( ReLambda() , 0.5)/max(S2, sqr(mag(U)/y_))
    );
    
    this->nut_ = fSS(Omega)*k_/max(omega_,this->omegaMin_) + nutL;
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
