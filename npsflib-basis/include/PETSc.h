#ifndef PETSc_H
#define PETSc_H

#include "petscts.h"
#include "petscksp.h"
#include <vector>
#include <iostream>

using namespace std; // to save the std:: in front of vector and cout

namespace PETSc_functions{

  void Create_Vec (Vec &vec, int nsize);
  void Get_Local_Info ( Vec vec, vector<int> &nlocal);
  void Create_Mat ( Mat &matrix, vector<int> nlocal, PetscInt n, PetscInt m );
  void Assembly_vector ( Vec &vec );
  void Assembly_matrix_final ( Mat &matrix );
  void Assembly_matrix_flush ( Mat &matrix );
  void Create_Viewer (PetscViewer &viewer, string filename );

}

#endif

