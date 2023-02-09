# Finite-element-method
## Physics Background
Solution of Time independent Schrodinger Eq using Finite element method. The system we solve here is Electrons-On-Helium system. As the name suggests,
the system is composed of an electron bound to a surface of superfluid Helium-4 (He). The binding surface potential consists of a short-range Pauli repulsion, keeping the electron out of the superfluid, and a long-range image-charge attraction, responsible for the electronic surface states. 
Since the atoms in the superfluid are adiabatic relative to the electron they can be modeled as "frozen". Thus simplification is used where the fluid is modeled by lattice atoms with the same density as the fluid atoms. My current research concentrate on the 

## Files in the repository 
### 2DTISE.py
For the 2 dimensional case the equation takes the form
$$(1) \\ H = -\frac{1}{2}(\partial_x^2 + \partial_z^2) + V$$
where the Hamiltonian is given in atomic units, ℏ=e=me=ϵ0=1.
The simple periodic potential 
$$(2) \\ V = \Theta(-z)\sin^{10}\left(\pi x/a\right)\sin^{10}\left(\pi z/a\right) - \Theta(z) \epsilon/z,$$
where ϵ is the relative dielectric constant of the superfluid, a is the lattice constant and Θ is the step function.  
This simple "toy" potential reproduces the same 1 eV barrier (or Gap in that case) and includes image potential. The solutions for z<0 can be generalized to 3D since we have x \leftright z symmetry between the x and z direction for z<0. 

### 3DTISE.py
Three-dimensional generalization of the 2DTISE file. The potential in this case is represented by a lattice of effective single electron potentials obtained from the interaction between the electron and the He-4 atoms, 
$$(3) \\ V(r) = \sum_{i \in cc} v^{He}(r-r_i),$$
where cc stands for cubic cell lattice and the single atom potential is (in atomic units), 
$$(4) \\ v_{He}(r) \approx a\left(b - \frac{1}{(r + 0.5)}\right) e^{-c \\ r} - \frac{\alpha}{2 r^4 + b}$$.
Here α is the Helium polarizability, a,b, and c can be derived from a more complicated form of effetive potential obtained by the pseudopotential method (Kestner (1967), Takada & Kohn (1987)) or Density-Functional-Theory (Whaley (1998)). The 1st term in Eq (4) comes from the attractive interaction of the electron with the nucleus and the Pauli repulsion (i.e. exclusion of Fermion statistics), and the second term comes mainly from the adiabatic dipole response of the Helium atom to the electric field of the electron.
