#include "nl-partsol.h"

/*******************************************************/

/*! 
  \fn static void Eigenerosion(int p,Fields Phi,Material MatPro,
  ChainPtr * Beps,double DeltaX)

  \brief Function to compute is a material point is or not eroded. 
  Here the notation is the same as in \cite Pandolfi_2012

  \param p : Index of the particle
  \param Phi : Particles fields
  \param Properties : Define the material properties of the particle
  \param B_eps : Define the particles close to each particle
  \param DeltaX : Mesh size
*/
static void Eigenerosion(
  int p, 
  Fields Phi, 
  Material MatPro,
  ChainPtr * Beps, 
  double DeltaX)
{

  /* Read the required fields */
  Matrix chi = Phi.chi;
  Matrix W = Phi.W;
  Matrix Mass = Phi.mass;
  Matrix Rho = Phi.rho;
  Matrix Stress = Phi.Stress;
  
  /*!
    Define auxiliar variable 
  */
  Tensor Stress_p, EV_Stress_p; /* Stress tensor */
  double Ceps_p, m_p, rho_p, V_p, sum_p, G_p, Gf_p;
  double m_q, rho_q, V_q, W_q;
  int * Beps_p;
  int NumBeps_p;
  int q;
  
  /*!
    For the current particle get first principal stress 
  */	
  Stress_p = memory_to_tensor__TensorLib__(Stress.nM[p], 2);
  EV_Stress_p = Eigenvalues__TensorLib__(Stress_p);

  /*!
    Non broken GP Only traction 
  */ 
  if((chi.nV[p] < 1) && 
    (EV_Stress_p.n[0]>0))
  { 

      /*!
	Normalizing constant 
      */
      Ceps_p = MatPro.Ceps;
    
      /*!
	Normalizing constant
      */
      Gf_p = MatPro.Gf;

      /*!
	Include the main particle in the calculus
      */
      m_p = Mass.nV[p];
      rho_p = Rho.nV[p];
      V_p = m_p/rho_p;
      
      sum_p = V_p*W.nV[p];

      /*!
	Get a pointer with the list of neighbours
      */
      NumBeps_p = lenght__SetLib__(Beps[p]);
      Beps_p = set_to_memory__SetLib__(Beps[p],NumBeps_p);

      /*!
	Loop over the neighbours 
      */
      for(int j = 0; j < NumBeps_p ; j++)
	{
	  /*!
	    Index 
	  */
	  q = Beps_p[j];

	  /*!
	    Get volume of particle q
	  */
	  m_q = Mass.nV[q];
	  rho_q = Rho.nV[q];
	  V_q = m_q/rho_q;

	  /*!
	    Get m_p 
	  */
	  V_p += V_q;

	  if(chi.nV[q] < 1)
	    {
	      /*!
		Get internal work of particle q 
	      */
	      W_q = W.nV[q]; 	
	      /*!
		Add to sum_p 
	      */
	      sum_p += V_q*W_q;
	    }
	
	}

      /*!
	Free memory 
      */
      free(Beps_p);

      /*!
	Compute energy-release rate for the particle 
      */
      G_p = (Ceps_p*DeltaX/V_p)*sum_p;

      /*!
	Fracture criterium
      */
      if(G_p > Gf_p)
	{
	  chi.nV[p] = 1.0;
	}
      
    }

  /*!
    Free eigenvalues 
  */
  free__TensorLib__(EV_Stress_p);
    
  
}

/*******************************************************/

/*!
  \fn static void Eigensoftening(int p,
  Fields Phi,
  Material Properties,
  ChainPtr * Beps)

  \brief Function to compute is a material point is or not eroded. 
  Here the notation is the same as in \cite Navas_2017_ES

  \param p : Index of the particle
  \param Phi : Particles fields
  \param Properties : Define the material properties of the particle
  \param Beps : Table with the list of neighbours per particle.
*/
static void Eigensoftening(
  int p, 
  Fields Phi, 
  Material MatPro, 
  ChainPtr * Beps)
{

  int Ndim = NumberDimensions;
  /* Read the required fields */
  Matrix chi = Phi.chi;
  Matrix Mass = Phi.mass;
  Matrix Rho = Phi.rho;
  Matrix Stress = Phi.Stress;
  Matrix Strain = Phi.Strain;
  Matrix Strain_If = Phi.Strain_If;
  
  /* Define auxiliar variable */
  Tensor Stress_p, Stress_q, EV_Stress_p, EV_Stress_q; 
  Tensor Strain_p, EV_Strain_p; 
  
  /* Material properties of the eigensoftening algorithm */
  double ft_p, Wc_p, heps_p;
    
  double m_p, sum_p, Stress_eps_p, chi_p, Strain_eps_f_p;
  double m_q;
  int * Beps_p;
  int NumBeps_p;
  int q;

  /* Get the tensile strengt of the material */
  ft_p = MatPro.ft;
  
  /* Get the bandwidth of the cohesive fracture (Bazant) */
  heps_p = MatPro.heps;
  
  /* Get the critical opening displacement */
  Wc_p = MatPro.Wc;
    
  /* Only for intact particles */
  if((chi.nV[p] == 0.0) && (Strain_If.nV[p] == 0.0))
    {

      /*!
	Get the number of neighbours 
      */
      NumBeps_p = lenght__SetLib__(Beps[p]);
      
      if(NumBeps_p > 0)
	{

	  /*!
	    Get the neighbours 
	  */      
	  Beps_p = set_to_memory__SetLib__(Beps[p],NumBeps_p);    
      
	  /*!
	    For the current particle get the volume 
	  */
	  m_p = Mass.nV[p];

	  /*!
	    For the current particle get first principal stress 
	  */	
	  Stress_p = memory_to_tensor__TensorLib__(Stress.nM[p], 2);
	  EV_Stress_p = Eigenvalues__TensorLib__(Stress_p);
      
	  /*!
	    Add the first term to the sumation 
	  */
	  sum_p = m_p*EV_Stress_p.n[0];

	  /*!
	    Free eigenvalues 
	  */
	  free__TensorLib__(EV_Stress_p);

	  /*!
	    Loop over the neighbours 
	  */
	  for(int j = 0; j < NumBeps_p ; j++)
	    {

	      /*!
		Get the indedx of each particle q close to p 
	      */
	      q = Beps_p[j];

	      /*!
		For the current particle get the volume
		and first principal stress 
	      */
	      m_q = Mass.nV[q];
	      
	      if(chi.nV[q] != 1.0)
		{
		  Stress_q = memory_to_tensor__TensorLib__(Stress.nM[q], 2);
		  EV_Stress_q = Eigenvalues__TensorLib__(Stress_q);
	    
		  /*!
		    Get sum_p 
		  */
		  sum_p += m_q*EV_Stress_q.n[0];

		  /* Free eigenvalues */
		  free__TensorLib__(EV_Stress_q);	    
		}
	
	      /*!
		Add volumen contribution 
	      */
	      m_p += m_q;
	    }
    
	  /* Free memory */
	  free(Beps_p);

	  /*!
	    Get the equivalent critical stress 
	  */
	  Stress_eps_p = sum_p/m_p;

	  /*!
	    Store the principal strain when crack start 
	  */
	  if(Stress_eps_p>ft_p)
	    {

	      Strain_p = memory_to_tensor__TensorLib__(Strain.nM[p], 2);
	      EV_Strain_p = Eigenvalues__TensorLib__(Strain_p);
	  	  
	      /*!
		Strain during fracture 
	      */
	      Strain_If.nV[p] = EV_Strain_p.n[0];

	      /* Free eigenvalues */
	      free__TensorLib__(EV_Strain_p);	    	  
	    }
	
	}
      
    }
  
  /*!
    Compute the damage parameter if the particle is damaged 
  */
  else if(chi.nV[p] != 1.0)
    {

      Strain_p = memory_to_tensor__TensorLib__(Strain.nM[p], 2);
      EV_Strain_p = Eigenvalues__TensorLib__(Strain_p);
	  	  
      /*!
	Fracture criterium 
      */
      Strain_eps_f_p = EV_Strain_p.n[0]-Strain_If.nV[p];
      chi_p = Strain_eps_f_p*heps_p/Wc_p;
      chi.nV[p] = DMIN(1,DMAX(chi_p,chi.nV[p]));


      /*!
	Update stress
       */
      for(int i = 0 ; i<Ndim*Ndim ; i++)
	{
	  Stress.nM[p][i] = (1-chi.nV[p])*Stress.nM[p][i];
	}

      /* Free eigenvalues */
      free__TensorLib__(EV_Strain_p);	    	  
      
    }
  /*!
    The particle is broken chi == 1.0
  */
  else
    {
      for(int i = 0 ; i<Ndim*Ndim ; i++)
	{
	  Stress.nM[p][i] = 0.0;
	}
    }
 
}

/*******************************************************/

/*
  \fn static void ComputeBeps(int p, Particle Particles, Mesh Nodes)
  
  \param Generate the B$_{\epsilon}$ neibourhood of each particle

  \param p : Particle
  \param Particles : Information of the particle mesh
  \param Nodes : Informaction with the set of nodes
 */
static void ComputeBeps(
  int p, 
  Particle MPM_Mesh, 
  Mesh FEM_Mesh)
{

  int Ndim = NumberDimensions;
  int Mat_p = MPM_Mesh.MatIdx[p];
  int I0 = MPM_Mesh.I0[p];
  
  /* Search radious */
  double epsilon = MPM_Mesh.Mat[Mat_p].Ceps*FEM_Mesh.DeltaX;

  ChainPtr Set_NodesBeps = NULL;
  int * NodesBeps;
  int NumNodesBeps;

  /* Index of each node close to the particle */
  int I_Beps;
  /* Interator pointer in Beps */
  ChainPtr Particles_Beps = NULL; 
  /* Index of a particle close to the particle p */
  int q_Beps;
  
  /* Distance */
  Matrix x_GC = MPM_Mesh.Phi.x_GC;
  Matrix X_p = memory_to_matrix__MatrixLib__(Ndim,1,x_GC.nM[p]);
  Matrix X_q = memory_to_matrix__MatrixLib__(Ndim,1,NULL);
  Matrix Distance;

  /* Get nodes close to the particle */
  Set_NodesBeps = FEM_Mesh.NodalLocality[I0];
  NumNodesBeps = FEM_Mesh.SizeNodalLocality[I0];
  NodesBeps = set_to_memory__SetLib__(Set_NodesBeps,NumNodesBeps);

  /* Loop in the nodes close to the particle */
  for(int i = 0 ; i<NumNodesBeps ; i++){

    /* Get the index of nodes close to the particle */
    I_Beps = NodesBeps[i];

    /* List of particles close to the node */
    Particles_Beps = FEM_Mesh.I_particles[I_Beps];
    
    while(Particles_Beps != NULL){

      /* Get the index of each particle */
      q_Beps = Particles_Beps->I;

      /* In Beps only those particles of the same material */
      if(Mat_p == MPM_Mesh.MatIdx[q_Beps]){

	/* Get the vector with the coordinates of each particle */
	X_q.nV = x_GC.nM[q_Beps];

	/* Get a vector from the GP to the node */
	Distance = substraction__MatrixLib__(X_p,X_q);

	/* Asign to p only those particles in Beps */
	if (norm__MatrixLib__(Distance,2) < epsilon){
	  push__SetLib__(&MPM_Mesh.Beps[p],q_Beps);
	}

	/* Free distance vector */
	free__MatrixLib__(Distance);

      }

      /* Go to the next set of particles */
      Particles_Beps = Particles_Beps->next;
      
    }
    
  }

  /* Free pointer with the list of nodes close to the particle */
  free(NodesBeps);    
 
}

/*********************************************************************/

void compute_particle_Damage(
  int p, 
  Particle MPM_Mesh, 
  Mesh FEM_Mesh)
{

  int Ndim = NumberDimensions;
  int Mat_p = MPM_Mesh.MatIdx[p];  
  double DeltaX = FEM_Mesh.DeltaX;
  
  /* Get the material properties of the particle */
  Material MatProp = MPM_Mesh.Mat[Mat_p];

  /* Beps of all the particles */
  ChainPtr * Beps = MPM_Mesh.Beps;

  /* Select the eigenerosion algorithm */
  if(MatProp.Eigenerosion)
    {
      /* Update Beps of each particle p */
      ComputeBeps(p, MPM_Mesh, FEM_Mesh);
      
      /* Update the damage variable of the particle */
      Eigenerosion(p,MPM_Mesh.Phi,MatProp,Beps,DeltaX);

      /* Free the previous list and set to NULL */
      free__SetLib__(&MPM_Mesh.Beps[p]);
    }

  /* Select the eigensoftening algorithm */
  if(MatProp.Eigensoftening)
    {
      /* Update Beps of each particle p */
      ComputeBeps(p, MPM_Mesh, FEM_Mesh);

      /* Update the damage variable of the particle */
      Eigensoftening(p,MPM_Mesh.Phi,MatProp,Beps);

      /* Free the previous list and set to NULL */
      free__SetLib__(&MPM_Mesh.Beps[p]);   
    }
  
}

/*******************************************************/
