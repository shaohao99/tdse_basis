#include "PETSc.h"
#include "petscts.h"
#include "petscksp.h"

//using namespace std;

namespace PETSc_functions{

  void Create_Vec (Vec &vec, int nsize)
  {
   VecCreate(PETSC_COMM_WORLD, &vec);
   VecSetSizes(vec, PETSC_DECIDE, nsize);
   VecSetFromOptions(vec);   // used by setting "-vec_type mpi" at run time
  }

  void Get_Local_Info ( Vec vec, vector<int> &nlocal)
  {
      PetscInt rstart, rend;
      VecGetOwnershipRange(vec, &rstart, &rend);

      PetscInt local_number;
      VecGetLocalSize(vec, &local_number);

      nlocal.resize(3);
      nlocal[1] = rstart;   // start index
      nlocal[2] = rend;     // end index
      nlocal[0] = local_number;  // number of local data

        // cout<<"nlocal="<<nlocal[0]<<'\t'<<"rstart="<<nlocal[1]<<'\t'<<"rend="<<nlocal[2]<<endl;
  }


  void Create_Mat ( Mat &matrix, vector<int> nlocal, PetscInt n, PetscInt m )
  {
      MatCreate(PETSC_COMM_WORLD, &matrix);
      PetscInt local_number=nlocal[0];
      MatSetSizes(matrix, local_number, local_number, n, m);
      MatSetFromOptions(matrix);   // used by set "-mat_type mpiaij" at run time, default mpiaij
      // MatSetType(matrix, MATMPIDENSE);  // MATMPIBAIJ, MATMPIDENSE
      // MatSetBlockSize(matrix, wf.n1);
  }


  void Assembly_vector ( Vec &vec )
  {
     VecAssemblyBegin(vec);
     VecAssemblyEnd(vec);
  }

  void Assembly_matrix_final ( Mat &matrix )
  {
     MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY);
     MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY);
  }

  void Assembly_matrix_flush ( Mat &matrix )
  {
     MatAssemblyBegin(matrix, MAT_FLUSH_ASSEMBLY);
     MatAssemblyEnd(matrix, MAT_FLUSH_ASSEMBLY);
  }


  void Create_Viewer (PetscViewer &viewer, string filename )
  {
      // PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer);   
      PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
      PetscViewerSetType(viewer, PETSC_VIEWER_ASCII);
      PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_SYMMODU);   // PETSC_VIEWER_ + DEFAULT, ASCII_MATLAB, ASCII_COMMON, ASCII_INDEX, ASCII_SYMMODU

      char name[filename.size()];   // convert string to char
      for(int i=0; i<filename.size(); i++){
         name[i] = filename[i];
      }
        // cout<<name<<endl;
      PetscViewerFileSetName(viewer, name);
  }


}  // end namespace
