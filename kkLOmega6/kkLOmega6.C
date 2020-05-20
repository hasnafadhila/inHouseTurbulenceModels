/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2015 OpenFOAM Foundation
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

#include "kkLOmega6.H"
//#include "fvOptions.H"
#include "bound.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
template<class BasicTurbulenceModel>
tmp<volScalarField> kkLOmega6<BasicTurbulenceModel>::lambdaT() const
{
    return (sqrt(kt_)/(omega_ + this->omegaMin_));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> kkLOmega6<BasicTurbulenceModel>::lambdaEff(const volScalarField& lambdaT) const
{
    return (min((Clambda_*y_), lambdaT));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> kkLOmega6<BasicTurbulenceModel>::Ret(const volScalarField& fw) const
{
    return (sqr(fw)*kt_/this->nu()/(omega_ + this->omegaMin_));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> kkLOmega6<BasicTurbulenceModel>::fnu(const volScalarField& Ret) const
{
    return (1.0 - exp(-sqrt(Ret)/Anu_));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kkLOmega6<BasicTurbulenceModel>::fINT() const
{
    return
    (
            min
            (
                kt_/(Cint_*(kl_ + kt_ + this->kMin_)),
                dimensionedScalar("1.0", dimless, 1.0)
            )
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kkLOmega6<BasicTurbulenceModel>::fSS(const volScalarField& Omega) const
{
    return (exp(-sqr(Css_*this->nu()*Omega/(kt_ + this->kMin_))));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kkLOmega6<BasicTurbulenceModel>::Cmu(const volScalarField& S) const
{
    return (1.0/(A0_ + As_*(S/(omega_ + this->omegaMin_))));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kkLOmega6<BasicTurbulenceModel>::BetaTS(const volScalarField& ReOmega) const
{
    return (scalar(1) - exp(-sqr(max(ReOmega - CtsCrit_, scalar(0)))/Ats_));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kkLOmega6<BasicTurbulenceModel>::fTaul
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
            -CTaul_*ktL
            /
            (
                sqr
                (
                    lambdaEff*Omega
                    + dimensionedScalar
                    (
                        "ROOTVSMALL",
                        dimLength*inv(dimTime),
                        ROOTVSMALL
                    )
                )
            )
        )
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kkLOmega6<BasicTurbulenceModel>::alphaT
(
    const volScalarField& lambdaEff,
    const volScalarField& fnu,
    const volScalarField& ktS
) const
{
    return (fnu*CmuStd_*sqrt(ktS)*lambdaEff);
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kkLOmega6<BasicTurbulenceModel>::fOmega
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
                    lambdaT
                    + dimensionedScalar
                    (
                        "ROOTVSMALL",
                        lambdaT.dimensions(),
                        ROOTVSMALL
                    )
                )
            )
        )
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kkLOmega6<BasicTurbulenceModel>::phiBP(const volScalarField& Omega) const
{
    return
    (
        //min
        //(
            max
            (
                kt_/this->nu()
             / (
                    Omega
                  + dimensionedScalar
                    (
                        "ROOTVSMALL",
                        Omega.dimensions(),
                        ROOTVSMALL
                    )
                )
              - CbpCrit_,
                scalar(0)
            )
            //scalar(50.0)
        //)
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kkLOmega6<BasicTurbulenceModel>::phiNAT
(
    const volScalarField& ReOmega,
    const volScalarField& fNatCrit
) const
{
    return
    (
        max
        (
            ReOmega
          - CnatCrit_
            / (
                fNatCrit + dimensionedScalar("ROOTVSMALL", dimless, ROOTVSMALL)
            ),
            scalar(0)
        )
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kkLOmega6<BasicTurbulenceModel>::D(const volScalarField& k) const
{
    return this->nu()*magSqr(fvc::grad(sqrt(k)));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> kkLOmega6<BasicTurbulenceModel>::fw
(
    const volScalarField& lambdaEff,
    const volScalarField& lambdaT
) const
{
    return
    (
        pow
        (
            lambdaEff
           /(lambdaT + dimensionedScalar("SMALL", dimLength, ROOTVSMALL)),
            2.0/3.0
        )
    );
}

template<class BasicTurbulenceModel>
tmp<volScalarField> kkLOmega6<BasicTurbulenceModel>::ktS
(
    const volScalarField& fSS,
    const volScalarField& fw
) const
{
    return(fSS*fw*kt_);
}

template<class BasicTurbulenceModel>
tmp<volScalarField> kkLOmega6<BasicTurbulenceModel>::nutS
(
    const volScalarField& fnu,
    const volScalarField& fINT,
    const volScalarField& Cmu,
    const volScalarField& ktS,
    const volScalarField& lambdaEff
) const
{
    return
    (
      fnu*fINT*Cmu*sqrt(ktS)*lambdaEff
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void kkLOmega6<BasicTurbulenceModel>::correctNut()
{
    // Currently this function is not implemented due to the complexity of
    // evaluating nut.  Better calculate nut at the end of correct()
    //notImplemented("kkLOmega6::correctNut()"); //Needs to be stopped in OF4x it's being called from function validate()
    // Someone is going to have to implement.... volunteers? Dhila? :-)
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kkLOmega6<BasicTurbulenceModel>::kkLOmega6
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

    A0_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "A0",
            this->coeffDict_,
            4.04
        )
    ),
    As_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "As",
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
    Abp_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Abp",
            this->coeffDict_,
            0.6
        )
    ),
    Anat_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Anat",
            this->coeffDict_,
            200
        )
    ),
    Ats_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ats",
            this->coeffDict_,
            200
        )
    ),
    CbpCrit_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CbpCrit",
            this->coeffDict_,
            1.2
        )
    ),
    Cnc_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cnc",
            this->coeffDict_,
            0.1
        )
    ),
    CnatCrit_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CnatCrit",
            this->coeffDict_,
            1250
        )
    ),
    Cint_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cint",
            this->coeffDict_,
            0.75
        )
    ),
    CtsCrit_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CtsCrit",
            this->coeffDict_,
            1000
        )
    ),
    CrNat_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CrNat",
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
            "C12",
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
    Css_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Css",
            this->coeffDict_,
            1.5
        )
    ),
    CTaul_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CTaul",
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
    CmuStd_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CmuStd",
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
    sigmaK_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaK",
            this->coeffDict_,
            1
        )
    ),
    sigmaOmega_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaOmega",
            this->coeffDict_,
            1.17
        )
    ),
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
            this->mesh_,
	    IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        kt_*omega_ + D(kl_) + D(kt_)
    ),
    SOut_
    (
      IOobject
      (
	  "S",
	  this->runTime_.timeName(),
	  this->mesh_,
	  IOobject::NO_READ,
	  IOobject::AUTO_WRITE
      ),
      this->mesh_,
      dimensionedSymmTensor("SOut", dimensionSet(0,0,-1,0,0,0,0), symmTensor::zero)
    ),
    OmegaOut_
    (
      IOobject
      (
	  "Omega",
	  this->runTime_.timeName(),
	  this->mesh_,
	  IOobject::NO_READ,
	  IOobject::AUTO_WRITE
      ),
      this->mesh_,
      dimensionedTensor("OmegaOut", dimensionSet(0,0,-1,0,0,0,0), tensor::zero)
    ),
    nutSOut_
    (
      IOobject
      (
	  "nuts",
	  this->runTime_.timeName(),
	  this->mesh_,
	  IOobject::NO_READ,
	  IOobject::AUTO_WRITE
      ),
      this->mesh_,
      dimensionedScalar("nutSOut", dimensionSet(0,2,-1,0,0,0,0), 0.0)
    ),
    epsilontsOut_
    (
      IOobject
      (
	  "epsilonts",
	  this->runTime_.timeName(),
	  this->mesh_,
	  IOobject::NO_READ,
	  IOobject::AUTO_WRITE
      ),
      this->mesh_,
      dimensionedScalar("epsilontsOut", dimensionSet(0,2,-3,0,0,0,0), 0.0)
    ),
    ktSOut_
    (
      IOobject
      (
	  "kts",
	  this->runTime_.timeName(),
	  this->mesh_,
	  IOobject::NO_READ,
	  IOobject::AUTO_WRITE
      ),
      this->mesh_,
      dimensionedScalar("ktSOut", dimensionSet(0,2,-2,0,0,0,0), 0.0)
    ),
    PklOut_
    (
      IOobject
      (
	  "Pkl",
	  this->runTime_.timeName(),
	  this->mesh_,
	  IOobject::NO_READ,
	  IOobject::AUTO_WRITE
      ),
      this->mesh_,
      dimensionedScalar("PklOut", dimensionSet(0,2,-3,0,0,0,0), 0.0)
    ),
    PktOut_
    (
      IOobject
      (
	  "Pkt",
	  this->runTime_.timeName(),
	  this->mesh_,
	  IOobject::NO_READ,
	  IOobject::AUTO_WRITE
      ),
      this->mesh_,
      dimensionedScalar("PktOut", dimensionSet(0,2,-3,0,0,0,0), 0.0)
    ),
    RetOut_
    (
      IOobject
      (
	  "Ret",
	  this->runTime_.timeName(),
	  this->mesh_,
	  IOobject::NO_READ,
	  IOobject::AUTO_WRITE
      ),
      this->mesh_,
      dimensionedScalar("Ret", dimensionSet(0,0,0,0,0,0,0), 0.0)
    ),
    fvOut_
    (
      IOobject
      (
	  "fv",
	  this->runTime_.timeName(),
	  this->mesh_,
	  IOobject::NO_READ,
	  IOobject::AUTO_WRITE
      ),
      this->mesh_,
      dimensionedScalar("fv", dimensionSet(0,0,0,0,0,0,0), 0.0)
    ),
    y_(wallDist::New(this->mesh_).y())
{
    bound(kt_, this->kMin_);
    bound(kl_, this->kMin_);
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
bool kkLOmega6<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel> >::read())
    {
        A0_.readIfPresent(this->coeffDict());
        As_.readIfPresent(this->coeffDict());
        Anu_.readIfPresent(this->coeffDict());
        Abp_.readIfPresent(this->coeffDict());
        Anat_.readIfPresent(this->coeffDict());
        Abp_.readIfPresent(this->coeffDict());
        Ats_.readIfPresent(this->coeffDict());
        CbpCrit_.readIfPresent(this->coeffDict());
        Cnc_.readIfPresent(this->coeffDict());
        CnatCrit_.readIfPresent(this->coeffDict());
        Cint_.readIfPresent(this->coeffDict());
        CtsCrit_.readIfPresent(this->coeffDict());
        CrNat_.readIfPresent(this->coeffDict());
        C11_.readIfPresent(this->coeffDict());
        C12_.readIfPresent(this->coeffDict());
        CR_.readIfPresent(this->coeffDict());
        CalphaTheta_.readIfPresent(this->coeffDict());
        Css_.readIfPresent(this->coeffDict());
        CTaul_.readIfPresent(this->coeffDict());
        Comega1_.readIfPresent(this->coeffDict());
        Comega2_.readIfPresent(this->coeffDict());
        Comega3_.readIfPresent(this->coeffDict());
        ComegaR_.readIfPresent(this->coeffDict());
        Clambda_.readIfPresent(this->coeffDict());
        CmuStd_.readIfPresent(this->coeffDict());
        Prtheta_.readIfPresent(this->coeffDict());
        sigmaK_.readIfPresent(this->coeffDict());
        sigmaOmega_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
void kkLOmega6<BasicTurbulenceModel>::correct()
{

    if (!this->turbulence_)
    {
        return;
    }

    // Local references//
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U_ = this->U_;
    const dimensionedScalar& kMin_ = this->kMin_;
    const dimensionedScalar& omegaMin_ = this->omegaMin_;
    const dimensionedScalar& epsilonMin_ = this->epsilonMin_;
    volScalarField& nut_ = this->nut_;
    volScalarField& omega_ = this->omega_;
    volScalarField& kt_ = this->kt_;
    volScalarField& kl_ = this->kl_;


    eddyViscosity<RASModel<BasicTurbulenceModel> >::correct();


    tmp<volTensorField> tgradU(fvc::grad(U_));
    const volTensorField& gradU = tgradU();

    const volSymmTensorField S(symm(gradU));
    const volTensorField W(skew(gradU));

    const volScalarField Omega(sqrt(2.0)*mag(skew(gradU)));
    const volScalarField S2(2.0*magSqr(dev(symm(gradU))));

    const volScalarField fwEff(fw(lambdaEff(lambdaT()),lambdaT()));
    const volScalarField ktSEff(ktS(fSS(Omega),fwEff));
    const volScalarField RetEff(Ret(fwEff));
    const volScalarField nutSEff(nutS(fnu(RetEff),fINT(),Cmu(sqrt(S2)),ktSEff,lambdaEff(lambdaT())));
    const volScalarField Pkt(nutSEff*S2);

    const volScalarField ktL(kt_ - ktSEff);
    const volScalarField ReOmega(sqr(y_)*Omega/this->nu());
    const volScalarField yEff(lambdaEff(lambdaT())/Clambda_);
    const volScalarField ReOmegaEff(sqr(yEff)*Omega/this->nu());
    const volScalarField nutl
    (
        min
        (
            ((fTaul(lambdaEff(lambdaT()), ktL, Omega)*C11_*Omega*sqr(lambdaEff(lambdaT()))
           *sqrt(ktL)*lambdaEff(lambdaT())/this->nu())
          //+ (BetaTS(ReOmega)*C12_*ReOmega*sqr(y_)*Omega))
          + (BetaTS(ReOmega)*C12_*ReOmegaEff*sqr(yEff)*Omega))
        ,
            (0.5*(kl_ + ktL)/(sqrt(S2) + omegaMin_))
        )
    );

    const volScalarField Pkl(nutl*S2);

    const volScalarField alphaTEff
    (
        alphaT(lambdaEff(lambdaT()), fnu(RetEff), ktSEff)
    );

    const volScalarField Rbp
    (
        CR_
        *(1.0 - exp(-phiBP(Omega)/Abp_)) //betaBP
	*kl_
	*omega_
       /(fwEff + ROOTVSMALL)
    );

    const volScalarField fNatCrit(1.0 - exp(-Cnc_*sqrt(kl_)*y_/this->nu()));

    const volScalarField Rnat
    (
        CrNat_
        *(1.0 - exp(-phiNAT(ReOmega, fNatCrit)/Anat_)) //betaNAT
	*kl_
	*Omega
    );


    const volScalarField epsilonts((this->omega_ * ktSEff) + D(ktSEff)); //THIS IS THE CURRENT ONE


    // ***************************** OUTPUT TERMS ********************************************
    // ***************************************************************************************
    // ***************************************************************************************

    SOut_ = S;
    OmegaOut_ = W;

    nutSOut_ = nutSEff;
    epsilontsOut_ = epsilonts;
    ktSOut_ = ktSEff;
    PklOut_ = Pkl;
    PktOut_ = Pkt;

    RetOut_ = RetEff;

    const volScalarField fvEff(fnu(RetEff));

    fvOut_ = fvEff;
    // ******************************************************************************************
    // ******************************************************************************************
    // ******************************************************************************************


    //TRANSPORT EQUATIONS

    omega_.boundaryFieldRef().updateCoeffs();

    tmp<fvScalarMatrix> omegaEqn
    (
     fvm::ddt(alpha, rho, omega_)
      + fvm::div(alphaRhoPhi, omega_)
      - fvm::laplacian(alpha*rho*DomegaEff(alphaTEff), omega_)
     ==
        alpha*rho*Comega1_*Pkt*omega_/(kt_ + kMin_)
      - fvm::SuSp
        (
            alpha*rho*(1.0 - ComegaR_/(fwEff + ROOTVSMALL))*(Rbp + Rnat)/(kt_ + kMin_)
          , omega_
        )
      - fvm::Sp(alpha*rho*Comega2_*sqr(fwEff)*omega_, omega_)
      + alpha*rho*(Comega3_*fOmega(lambdaEff(lambdaT()), lambdaT())*alphaTEff*sqr(fwEff)*sqrt(kt_))()()/pow3(y_())
    );

    omegaEqn.ref().relax();
    omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
    solve(omegaEqn);
    bound(omega_, this->omegaMin_);


    const volScalarField Dl(D(kl_));

    // Laminar kinetic energy equation
    tmp<fvScalarMatrix> klEqn
    (
     fvm::ddt(alpha, rho, kl_)
      + fvm::div(alphaRhoPhi, kl_)
      - fvm::laplacian(alpha*rho*this->nu(), kl_)
     ==
        alpha*rho*Pkl
      - fvm::Sp(alpha*rho*((Rbp + Rnat + Dl)/(kl_ + kMin_)), kl_)
    );

    klEqn.ref().relax();
    klEqn.ref().boundaryManipulate(kl_.boundaryFieldRef());
    solve(klEqn);
    bound(kl_, this->kMin_);


    const volScalarField Dt(D(kt_));

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> ktEqn
    (
     fvm::ddt(alpha, rho, kt_)
      + fvm::div(alphaRhoPhi, kt_)
      - fvm::laplacian(alpha*rho*DkEff(alphaTEff), kt_)
     ==
        alpha*rho*Pkt
      + alpha*rho*(Rbp + Rnat)
      - fvm::Sp(alpha*rho*(omega_ + Dt/(kt_+ kMin_)), kt_)
    );

    ktEqn.ref().relax();
    ktEqn.ref().boundaryManipulate(kt_.boundaryFieldRef());
    solve(ktEqn);
    bound(kt_, this->kMin_);


    // Update total fluctuation kinetic energy dissipation rate
    epsilon_ = kt_*omega_ + Dl + Dt;
    bound(epsilon_, epsilonMin_);


    // Re-calculate turbulent viscosity
    nut_ = nutSEff + nutl;
    nut_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
