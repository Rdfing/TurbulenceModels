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

#include "kOmega20062D.H"
#include "fvOptions.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void kOmega20062D<BasicTurbulenceModel>::correctNut()
{
    //this->nut_ = k_/omega_; //orgional
    this->nut_ = k_/ max(omega_, Clim_*sqrt(2.0/0.09*magSqr(symm(fvc::grad(this->U_)))));
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kOmega20062D<BasicTurbulenceModel>::kOmega20062D
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

    Cmu_ //This is beta_star
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaStar",
            this->coeffDict_,
            0.09
        )
    ),
    // beta_
    // (
    //     dimensioned<scalar>::lookupOrAddToDict
    //     (
    //         "beta",
    //         this->coeffDict_,
    //         0.072
    //     )
    // ),
    gamma_ // not use alpha since alphas has been used by the volume fraction
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma",
            this->coeffDict_,
            0.52 // 13/25
        )
    ),

    sigmaK_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK",
            this->coeffDict_,
            0.6
        )
    ),

    sigmaOmega_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega",
            this->coeffDict_,
            0.5
        )
    ),

    Clim_ //added
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Clim",
            this->coeffDict_,
            0.875 //Clim=7/8
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

    // added
    // nut_
    // (
    //     IOobject
    //     (
    //         "nut",
    //         runTime_.timeName(),
    //         mesh_,
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     autoCreateNut("nut", mesh_)
    // )

    fBeta_
    (   
        IOobject
        (
            "fBeta",
            this->runTime_.timeName(),
            this->mesh_
        ),
        this->mesh_,
        dimless
    ),

    Chi_
    (   
        IOobject
        (
            "Chi",
            this->runTime_.timeName(),
            this->mesh_
        ),
        this->mesh_,
        dimless
    ),

    absChi_
    (   
        IOobject
        (
            "absChi",
            this->runTime_.timeName(),
            this->mesh_
        ),
        this->mesh_,
        dimless
    ),
    beta_
    (   
        IOobject
        (
            "beta",
            this->runTime_.timeName(),
            this->mesh_
        ),
        this->mesh_,
        dimless
    ),
    sigmaD_
    (   
        IOobject
        (
            "alphad",
            this->runTime_.timeName(),
            this->mesh_
        ),
        this->mesh_,
        dimensionedScalar("alphad", dimless, 0.125)
    )
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
bool kOmega20062D<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        // beta_.readIfPresent(this->coeffDict());
        gamma_.readIfPresent(this->coeffDict());
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
void kOmega20062D<BasicTurbulenceModel>::correct()
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

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

/**********************************************************************/
    // added
    volSymmTensorField Sij(symm(fvc::grad(U)));
    volTensorField Omij(-skew(fvc::grad(U)));
/**********************************************************************/
    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField G
    (
        this->GName(),
        nut*(tgradU() && dev(twoSymm(tgradU())))
    );
    tgradU.clear();

    // Update omega and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();
/**********************************************************************/

    // sigmaD switch
    volScalarField sigmaDCheck(fvc::grad(k_) & fvc::grad(omega_)); 

    forAll(sigmaD_,celli)
    {
        if (sigmaDCheck[celli] <= 0.0) //if (sigmaDCheck[celli] <= 0.0001)
        {
            sigmaD_[celli]=scalar(0);
        }
        else
        {
            sigmaD_[celli]=scalar(0.125); //1/8
        }
    }

    // coef in the cross diffusion
    volScalarField CDkOmega
    (
        (1*sigmaD_)*(fvc::grad(k_) & fvc::grad(omega_))/omega_
    );

    // for 3D model
    Chi_ = (Omij & Omij) && Sij /pow((Cmu_*omega_),3); // does not do anything for 2D
    absChi_ = mag(Chi_); // does not do anything for 2D
    //fBeta_ = (1.0+85.0*absChi_)/(1.0+100.0*absChi_);

    // for 2D
    fBeta_ = 1.0; 

    beta_ = 0.0708*fBeta_;


    // Turbulence specific dissipation rate equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(alpha, rho, omega_)
      + fvm::div(alphaRhoPhi, omega_)
      - fvm::laplacian(alpha*rho*DomegaEff(), omega_)
     ==
        gamma_*alpha*rho*G*omega_/k_ // see note
      - fvm::SuSp(((2.0/3.0)*gamma_)*alpha*rho*divU, omega_) // see note
      - fvm::Sp(beta_*alpha*rho*omega_, omega_) // implicit treatment of the negtive term will enhace the diagnal dominance of the matrix
      - fvm::SuSp(alpha*rho*(- scalar(1))*CDkOmega()/omega_(),omega_) // cross diffusion term in 2006
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
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha*rho*G // wiithout production limiter
        //alpha*rho*min(G, (20*Cmu_)*k_*omega_) //with production limiter
      - fvm::SuSp((2.0/3.0)*alpha*rho*divU, k_)
      - fvm::Sp(Cmu_*alpha*rho*omega_, k_)
      + fvOptions(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    // update viscosity
    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
