# Finite-element-method
## Physics Background
Solution of Time independent Schrodinger Eq using Finite element method. The system we solve here is Electrons-On-Helium system. As the name suggests,
the system is composed of an electron bound to a surface of superuid Helium-4 (He). The binding surface potential consists of a short-range Pauli repulsion, keeping the electron out of the superfluid, and a long-range image-charge attraction, responsible for the electronic surface states. 
Since the atoms in the superfluid are adiabatic with respet to the electron they  can be modeled as frozen. Thus simplification is used where the fluid is modeled by lattice atoms with the same density as the fluid atoms. My current research concentrate on the 

## Files in repository 
### 2DTISE.py
For the 2 domentional case the equation takes the form
$$(1) \; H = -\frac{1}{2}(\partial_x^2+ \partial_z^2) + V$$
where the Hamiltonian is given in atomic units, $\hbar = e = m_e = \epsilon_0 = 1.$
The simple periodic potential 
$$(2) V = \Theta(-z)\sin^{10}\left(\pi x/a\right)\sin^{10}\left(\pi z/a\right) - \Theta(z) \epsilon/z,$$
where $\epsilon$ is the relatice dielectric constant of the superfluid, a is the lattice constant and $\Theta$ is the step function.  
This toy potential reproduces the same 1 eV barrier (or Gap in that case) takes the form and includesbinding (for energies E < 0) image potantial. The solutions for z<0 can be genelralized to 3D since we have symmtery between the x and z direction for z<0.  

### 3DTISE.py
Three dimentional generalization of the 2DTISE file. The potential in this case is represented by lattice of effective single atom potentials the electcron experience from a single Helium-4 atoms, 
$$(3) \; V(r) = \sum_{i \in cc} v^{He}(r-r_i)$$,
where cc stands for cubic cell lattice and the single atom potential is (in atomic units), 
$$(4) \; v^{He}(r) \approx a(b - 1/r) e^{-c \; r} - \frac{\alpha}{2 r^4 + b}$$.
The simple form of the potential and the constants a,b,c can be derived from more complicated form obtained by the pseudopotential method (Kestner (1967), Takada & Kohn (1987)) or Density functional theory.

