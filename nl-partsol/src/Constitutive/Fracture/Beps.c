/**
 * @file Beps.c
 * @author Miguel Molinos (@migmolper)
 * @brief
 * @version 0.1
 * @date 2022-05-25
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "Constitutive/Fracture/Beps.h"
#include "Globals.h"

void compute_Beps__Constitutive__(Particle MPM_Mesh, Mesh FEM_Mesh,
                                  bool Initialize_Beps) {

  unsigned Np = MPM_Mesh.NumGP;
  unsigned p;
  double DeltaX = FEM_Mesh.DeltaX;

/*
  Loop in the material point set
*/
#pragma omp parallel
  {
#pragma omp for private(p)
    for (p = 0; p < Np; p++) {

      //! Decide if compute Beps or not
      const double *dis_p = MPM_Mesh.Phi.dis.nM[p];
      if ((euclidean_norm__TensorLib__(dis_p) > 0.000001) ||
          (Initialize_Beps == true)) {
        free__SetLib__(&MPM_Mesh.Beps[p]);
      } else {
        continue;
      }

      //! Compute the normalizing distance
      unsigned MatIndx_p = MPM_Mesh.MatIdx[p];
      Material MatProp_p = MPM_Mesh.Mat[MatIndx_p];
      double Ceps_p = MatProp_p.Ceps;
      double eps_distance_p = Ceps_p * DeltaX;

      //! Position of the particle p
      const double *xp = MPM_Mesh.Phi.x_GC.nM[p];

      //! Loop over the list of nodes close to the particle
      unsigned I0_p = MPM_Mesh.I0[p];
      ChainPtr Nodal_Locality_p = FEM_Mesh.NodalLocality_0[I0_p];

      while (Nodal_Locality_p != NULL) {

        unsigned A = Nodal_Locality_p->Idx;

        ChainPtr Particle_Locality_I = FEM_Mesh.List_Particles_Node[A];

        //! Loop over the list of particles close to each node
        while (Particle_Locality_I != NULL) {
          unsigned q = Particle_Locality_I->Idx;

          //! Position of the particle q
          const double *xq = MPM_Mesh.Phi.x_GC.nM[q];

          //! Update list of particles close to p
          if (euclidean_distance__TensorLib__(xp, xq) <= eps_distance_p) {
            push__SetLib__(&MPM_Mesh.Beps[p], q);
          }

          Particle_Locality_I = Particle_Locality_I->next;

        } //  while Particle_Locality_I

        Nodal_Locality_p = Nodal_Locality_p->next;

      } // while Nodal_Locality_p

    } // for p
  }   // #pragma omp parallel
}
