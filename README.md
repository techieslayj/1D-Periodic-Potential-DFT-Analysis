# 1D-Periodic-Potential-DFT-Analysis
Ongoing research/Master's Project involving quantum chemistry and DFT fundamentals

* Code done in MATLAB however models and comparison analysis were created/done using NWChem partnerred with parallel computing modules such as OpenMP/OpenMPI.


This is a multi-step Density Functional Theory (DFT) analysis that includes multiple code files and layers of mathematical/physical/chemical analysis. Basic quantum fundamentals are used such as Schrodinger, Hartree-Fock, Kohn-Sham and Bloch's theorem to name a few. I am continuously contributing to this project as it is a part of my Master's project I am putting together for graduate school.

I begin this repository by posting some of the more important code files for my analysis: 

# System background:
  1D system that is purely periodic in behavior. Preset box or cell size of 12 (NPTS) where all atoms within the system are hydrogens (1.0) which means there are also 12 electrons in the system respectively. Hydrogen was a perfect starting atom to use for this system because its energy levels and chemical makeup are relatively easy to use for our purposes of DFT. Bond lengths are also very easy to calculate (equilibrium distance of an alternating hydrogen chain). This system make up is created and made ready for use for DFT calculations within the file input_KSDFT.m.
  
# Creating a working periodic potential system (Hamiltonian/Central Finite Difference)
  One of the challenges of creating a strictly 1D periodic potential box or cell in solid state is accounting for Bloch's theorem. Which states that energy eigenstates for electron within a solid state can be seen as bloch waves. In periodic systems it is important to note that the solutions to our Schrodinger equations should be characterized by an integer number and a vector. In other words the periodicity of the system restricts the nonzero components of the ionic potential to reciprocal-lattice vectors. To account for this within my program I use simple central finite difference technique to construct a banded, symmetric Hamiltionian matrix to calculate my prospective eigenvalues and eigenvector. In my case I specifically use fourth-order Finite Difference constants to load my Hamiltonian matrix along with my step size around the selected cell size. We should expect to see a matrix that is symmetric therefore some of these constants should move to the right side of the cell if it begins outside of our boundary conditions. This step allowed us to get a plot with periodic density values for our electrons in the system.  
  
# Creating a working periodic potential system (Eigenvector + Eigenvalues)
  With a working, well-constructed hamiltonian we were now ready to proceed to obtaining our eigenvalues and eigenvector. This was extremely simple to do in MATLAB using the eig() function which returns both the solutions we were looking for. However it was important to normalize the eigenvector after sorting it correctly. 
  
# Evaluating our system and plotting (Still needs work to fully be accurate with 1D periodic system)
  We are now ready to perform the DFT calculation in order to obtain a plot displaying our atoms (displayed as red dots), density, external potential and our Kohn-Sham potential (KS potential). Density calculations are already fine tuned for a 1D periodic system however I am currently working on changing the external potential of the system in order to create a plot fully resembling a 1D periodic system.
  
# The 1D vs. 3D system Issue for External Potential (Fourier Transformation)
  As previously mentioned above our external potential is currently NOT periodic when we run the input program, why? Well this is a relatively easy issue to identify but a little tricky to solve. We are performing a fourier transformation to our system structure factor at a corresponding wave vector (defined mathematically as exp(i*G*R)) as well as the pseudopotential of the overall system. Pure coulumb potential is defined by Z/r where Z is the valence of the atom. In a 3D system we can multiply this integral (which is the difference of Coulumb potential and pseudopotential by 4*pi*r^2*dr, which correspondingly will cause the r in the Z/r expression to cancels and we eliminate the worry of the system diverging at r=0. However, in a 1D system we do not have a r^2 mathematically to work with which would cause Z to still be divided by r. What does this mean? This is a huge problem because our system would diverge at r=0 meaning simply removing the r^2 value to resemble the system as a 1D space is not sufficient. I am currently working on this step of the project.***
  
# Still to come:

  * K-point sampling for Brillouin zones (important contribution of Bloch's Theorem, consider only k within cell in reciprocal space)
  * Prove/Disprove Peierls' Theorem (states a one-dimensional equally spaced chain with one electron per ion is unstable)
  * Expand to other atoms
  
# Special Thanks 
  * Dr. Chen Huang Florida State University
  
  
  
  
