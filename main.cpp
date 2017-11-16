#include <iostream>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <string>

// Change this next line to your PATH:
#include "/home/asdf/gems/gmml/includes/MolecularModeling/assembly.hpp"
//#include "/home/oliver/Programs/gems/gmml/includes/gmml.hpp"

using namespace MolecularModeling;
using namespace gmml;
typedef std::vector<Atom*> AtomVector;

int main(){

    Assembly example1("example1.pdb", PDB);
    example1.BuildStructureByDistance(); // Sets the bonding information in AtomNode based on distance

    Atom *atom1, *atom2, *atom3, *atom4; //Pointers to the atoms we will rotate
    //Go through each atom, and select the four that define the dihedral we want to rotate.
    AtomVector atoms = example1.GetAllAtomsOfAssembly();
    for(AtomVector::iterator it1 = atoms.begin(); it1 != atoms.end(); ++it1)
    {
        Atom *atom = *it1; // Makes readability and my life easier, as the IDE will auto-suggest Atom methods.

        if ( (atom->GetName().compare("C1")==0) && (atom->GetResidue()->GetName().compare("0SA")==0) )
        {
            atom1 = atom;
        }
        if ( (atom->GetName().compare("C2")==0) && (atom->GetResidue()->GetName().compare("0SA")==0) )
        {
            atom2 = atom;
        }
        if ( (atom->GetName().compare("O3")==0) && (atom->GetResidue()->GetName().compare("3LB")==0) )
        {
            atom3 = atom;
        }
        if ( (atom->GetName().compare("C3")==0) && (atom->GetResidue()->GetName().compare("3LB")==0) )
        {
            atom4 = atom;
        }
    }

    //I should check if they are set, but I won't
    std::cout << atom1->GetId() << ", " << atom2->GetId() << ", " << atom3->GetId() << ", " << atom4->GetId() << std::endl;
    example1.SetDihedral(atom1, atom2, atom3, atom4, 60.0);

    // Write out a pdb file with the rotated dihedral. Can view this file in a program like VMD.
    PdbFileSpace::PdbFile *outputPdbFile = example1.BuildPdbFileStructureFromAssembly(-1,0);
    outputPdbFile->Write("example1_rotated.pdb");
}
