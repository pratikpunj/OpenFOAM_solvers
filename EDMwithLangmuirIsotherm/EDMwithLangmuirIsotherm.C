/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

Application
    EDMwithLangmuirIsotherm

Group
    custom

Description
    Binding condition solver for the application of Liquid chromatoraphy with Langmuir Isotherm and first order adsorption kinetics.


    \vartable
        C       | Concentration
        q       | Adsorbed concentration
        U       | Superficial/Geometrical Velocity
	D	| Dispersion Tensor
	eps	| Porous zone total porosity
	q_s	| Maximum concentration of solute adsorbed
	b	| Langmuri equilibrium constant
	k	| First order adsorption rate constant

    \endvartable

    \heading Required fields
    \plaintable
        C       | Concentration
        q       | Adsorbed concentration
    \endplaintable

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "simpleControl.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Non Binding contration distribution solver with analytical velocity field."
    );

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating Concentration Field\n" << endl;

    

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
        
        #include "CourantNo.H"

        while (simple.correctNonOrthogonal())
        {
            tmp<volScalarField> q_e = (q_s * b * C) / (1.0 + (b*C));
	       
	    fvScalarMatrix qEqn
	    (
    		fvm::ddt(q)
  		+ fvm::Sp(k,q)
	        == k*q_e 
	    );

	    qEqn.relax();
	    qEqn.solve();
	    q.correctBoundaryConditions();

	    fvScalarMatrix CEqn
            (
                fvm::ddt(C)
              + ((1/eps)*fvm::div(phi, C))
              - fvm::laplacian(D, C)
	      + ((1-eps)/eps)*fvc::ddt(q)
	      ==
	        fvOptions(C)
            );

            CEqn.relax();
	    fvOptions.constrain(CEqn);
            CEqn.solve();
	    fvOptions.correct(C);
            C.correctBoundaryConditions();
        }

        runTime.write();
        
        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}




// ************************************************************************* //
