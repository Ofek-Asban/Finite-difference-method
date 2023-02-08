import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import eigs
from scipy.sparse import kronsum
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from scipy.optimize import curve_fit

#Dimentions of the numerical solution are in Atomic units, i.e. e = hbar = K = m electron = 1. Thus all energies are in Hartrees and distances are in Bohrs.

'''a = lattice constant size in units of Bohr. 
d = discretization size (distance between dots)
MaxIdx = max index of eigenvalue, 
P = 0 or 1 (0 Boundary consition, i.e. infinte Wall) 1 (Periodic Boundary) in x direction 
'''
def H_2D(xmin, xmax, zmin, zmax, a, d, MaxIdx, P):  
    dx = dz = d*a; 
    x = np.arange(xmin * a, xmax * a + dx, dx) #Descretizing space in the x direction
    z = np.arange(zmin * a, zmax * a + dz, dz) #Descretizing space in the z direction   
    X, Z = np.meshgrid(x, z) # Set up the 2D grid of space
    V =  np.heaviside(-Z + a/2,0)*(np.sin(2*np.pi*X/(2*a))**10)*(np.sin(2*np.pi*(Z-a/2)/(2*a))**10)  - np.heaviside(Z-a/2,0)*(0.007/Z) #Here we make lattice potential in x,z directions, sin^10(x)sin^10(z) for z<0 and an attractive -1/z potential for z>0 (This is a toy model for electron on Helium that gives 1eV potential barrier at z = 0)
    Dx = diags([1, -2, 1, P, P], [-1, 0, 1, x.size-1, 1-x.size], shape=(x.size,x.size))/(dx**2); #2nd derivative matrix in real space for x direction
    Dz = diags([1, -2, 1, P, P], [-1, 0, 1, z.size-1, 1-z.size], shape=(z.size,z.size))/(dz**2); #2nd derivative matrix in real space for z direction
    T = - 0.5 * kronsum(Dx,Dz) #kron sum of sparce matrices to reduce a tensor problem to regular matrix eigenvalue problem
    H = T + diags(V.reshape(-1)) # The Hamiltonian operator in real space
    energies, states = eigs(H, k = MaxIdx, which='SM') # Calculationg the MaxIdx lowest eigenstates 
    idx = np.argsort(energies) #array of index values sorted according to increasing value
    energies = energies[idx] #sorting eigenvalues (energies) array
    states = states[:,idx] #sorting eigenvactors according of increasing value of their eigenvalues
    states = states.reshape((X.shape[0],X.shape[1],states.shape[1])) # Reshaping the eigenstate vector to take the form a f(z,x) function.
    states = states/np.sqrt(dx*dz) #eig func returns normalizaed eigenvectors according to the descrete rule sum_i (v_i).T (v_i) = 1 but in real space eigenstgates has units of 1/(length^d/2) according to the L-2 normalization in real space sum_{ij} (v[i,j]^dag) * v[i,j] dx dz = sum_{ij} (v[j,i]^*) * v[i,j] dx dz = 1 thus v is normalized with 1/((dx dz)^1/2), note that hbar = m = e = 1 (A.U.) thus dz and dz are dimensionless (in dimensions of Bohr radius)   
    return X, Z, energies, states, V, MaxIdx

#Global constants of the model
xmin = -3; xmax = 3; zmin = -2.2; zmax = 1/0.007; a = 7; d = 0.05; # d is the resolution of z axis (0.05 = 20 points per "a" lattice const), xmin, xman, zmin, zman are given in units of the lattice consant a = 7 Bohr
X, Z, energies, states, V, MaxIdx = H_2D(xmin, xmax, zmin, zmax, a, d, 300, 1)


#############################################################################################################################


#Plotting energies
plt.plot(np.arange(0, energies.size, 1), 27.2*np.real(energies), marker='o')

#Finding the indices of the eigenstate matrix that belong to re rigion around z<=0 
Z_origin_Idx = int(abs(zmin/d)); Z_Neg_Idx = Z_origin_Idx + int(zmin/d); Z_Pos_Idx = Z_origin_Idx + int((0.5)/d) + 1; #Z_origin_Idx is the index of z coordinate vector at the origin (z=0), Z_Neg_Idx (Z_Pos_Idx) is the index of z vector that multiples of lattice cinst "a" in z<0 (out z>0).


#plotting around z=0 real space 2D potential in 3D
plt.axes(projection='3d').plot_surface(X[Z_Neg_Idx : Z_Pos_Idx], Z[Z_Neg_Idx : Z_Pos_Idx], V[Z_Neg_Idx : Z_Pos_Idx,:], cmap=cm.coolwarm)
#Plotting Full real space 2D eigenstates in 3D (3rd dim is value of wave func)
plt.axes(projection='3d').plot_surface(X, Z, (-states[:,:,0]), cmap=cm.coolwarm)
#Plotting real space 2D eigenstates in 3D near NEAR THE Z=0 BORDER
plt.axes(projection='3d').plot_surface(X[Z_Neg_Idx : Z_Pos_Idx,:], Z[Z_Neg_Idx : Z_Pos_Idx, :], (-states[Z_Neg_Idx : Z_Pos_Idx, : , 286]), cmap=cm.coolwarm)

#Plotting the ground state (eigenvector of the loest wigenvalue) as a function of z for different x values represented by i
state_avg = np.zeros(Z_Pos_Idx-Z_Neg_Idx)
Area =  np.array([]) # initializing the area array 
int(abs(zmin/d))
for i in range(Z.shape[1]):
    EigenVecTail_i = -states[Z_Neg_Idx : Z_Pos_Idx, i, 0]
    plt.plot(Z[Z_Neg_Idx : Z_Pos_Idx, 0], EigenVecTail_i); # Z[z range, Choosing 0 WLOG] -> Taking only the tail of the state near z=0 for specific x value
    state_avg = state_avg + EigenVecTail_i #Summing all the different of the function values at different x to obtain an averaged value function
    area_i = np.sum(np.abs(EigenVecTail_i)) #Finding the area under the tail of the i'th wave function
    Area = np.append(Area, area_i) #Adding area_i to the i'th of the Area array for later use 
state_avg = state_avg/Z.shape[1] #Dividing the summed tails by the num of funcitions to get an avraged tail hight


EigenVec_MaxTail, EigenVec_MinTail = -states[Z_Neg_Idx : Z_Pos_Idx, np.argmax(Area), 0], -states[Z_Neg_Idx : Z_Pos_Idx, np.argmin(Area), 0] #Saving the function's X locations which has the maximum/minimum area the z direction

Z_tail = np.array([state_avg, EigenVec_MaxTail, EigenVec_MinTail]) #Saving the tail of wave func in z direction for 3 diff cases: 1. avg over x axis, 2. x value that has max tail in z direction, 3. x value that has min tail in z direction 
for idx in range(0,Z_tail.shape[0]):
    plt.plot(Z[Z_Neg_Idx : Z_Pos_Idx,0], Z_tail[idx], marker = 'o'); #plotting the min, max, and average
#Note that OutHeIdx should start where the He lattice potential starts (in this case we set the lattice potential to start at half lattice constant z = a/2 = 3.5 Bohr)

# defining to a simple exponential fit to the wave function. 
# b = sqrt(potential barrier -> V0) and c is the shift of the potential wall where the wave function starts to decay 
def func(x, a, b, c):
    return a * np.exp(-b * np.abs(x-c))

#Fiting the numerical wave function to a simple exponential function for 2 difference slices and the average on x axis 1. The average 2. Max 3. Min 
for idx in range(0,Z_tail.shape[0]):
    popt, pcov = curve_fit(func, Z[Z_Neg_Idx : Z_Pos_Idx,0], Z_tail[idx])
    #Plotting the fit function on top of the numerical wave funcion (as can be seen the variation of the exponential power is approximatly V0 z and we can neglect the dependence V(z) of the barrier)
    plt.plot(Z[Z_Neg_Idx : Z_Pos_Idx,0], Z_tail[idx], marker = 'o');
    plt.plot(Z[Z_Neg_Idx : Z_Pos_Idx,0], func(Z[:,i][Z_Neg_Idx : Z_Pos_Idx], *popt))
    plt.title("Max, Man and Average tail of ground state wave function for z<0 in region of Lattice potential fitter to simple exponential form");plt.xlabel("z (Bohr Radius)"); plt.ylabel("Wave Function (1/Bohr)");
    # b = sqrt((2m/hbar)*V0) = in A.U = sqrt(2*V0) -> V0 = b^2/2 = in eV = 27.2 * b^2/2. V0 is the average energy scale in which the wave func decay exponentially and accuratly it is the: V0 = gap hight - E of the particle, however E ~ - 8K and gap hight ~ 1eV >> E thus V0 ~ Gap hight ~ 1eV    
    V0 = 27.2*(popt[1])**2/2 #pop[1] is the b parameter in the fit function, in atomic units (m=hbar=1) it is sqrt(potential barrier), 27.2 is conversion from Hartrees to eV.
    print(V0)

#########################Looking on an X slice of the wave func to see if close states are symmetric and anti-symmetric########################

#Plotting the taild of the GS and 1st excited states  
dx = dz = d*a
x = np.arange(xmin * a, xmax * a + dx, dx)
fig = plt.figure()
ax = fig.gca(projection='3d')
surf1 = ax.plot_surface(X[Z_Neg_Idx : Z_Pos_Idx], Z[Z_Neg_Idx : Z_Pos_Idx], states[Z_Neg_Idx : Z_Pos_Idx,:,0], cmap=cm.coolwarm)
surf2 = ax.plot_surface(X[Z_Neg_Idx : Z_Pos_Idx], Z[Z_Neg_Idx : Z_Pos_Idx], states[Z_Neg_Idx : Z_Pos_Idx,:,1], cmap=cm.coolwarm)
fig.show()


#Slice along the X axis of the wave func (on specific values in the Z axis where the tail is) 
for i in range(229,230,1):
    avg_i = sum(states[Z_origin_Idx,:,i])/(states[Z_origin_Idx,:,i].size) #averaging in order to colaps the graphs on each other
    Shifted_state_i =  states[Z_origin_Idx,:,i] - avg_i 
    plt.plot(X[Z_origin_Idx], Shifted_state_i)
    
#Slice along the Z axis of the wave func (on specific values in the X axis) 
for i in range(4):
    plt.plot(Z[:,0], states[:,0,i]) #states[all Z values, X cut , n level], Note that the modulation of x is small and happends only at the tail of wave func at z <= a/2 (where a/2 = 3.5 aB) so at the scale of localization length on z axis we do not see this difference 
    
#Chacking if the states are orthomormal
np.sum(states[:,:,0]**2)*dx*dz # Normalization (summing over z and x axis)
np.sum(np.real(states[:,:,10]*states[:,:,11]))*dx*dz # Orthogonality

#Plotting the state and its x derivative as a function x for specific value of z in the lattice region near z<-0 
Xgradient_batch = np.gradient(states[Z_Neg_Idx : Z_Pos_Idx,:,0], dx, axis=1) #axis=0,(1) means that the derivative is w.r.t the Z (X)
surf = ax.plot_surface(X[Z_Neg_Idx : Z_Pos_Idx], Z[Z_Neg_Idx : Z_Pos_Idx], states[Z_Neg_Idx : Z_Pos_Idx,:,0], cmap=cm.coolwarm)
surf_grad = ax.plot_surface(X[Z_Neg_Idx : Z_Pos_Idx], Z[Z_Neg_Idx : Z_Pos_Idx], Xgradient_batch, cmap=cm.coolwarm)
fig.show()
 
# Calculating the momentum matrix element p = -i Dx where Dx is the derivative matrix.
X_origin_Idx = int(abs(xmin/d)); X_Neg_Idx = int(X_origin_Idx - 1/(2*d)) ; X_Pos_Idx = int(X_origin_Idx + 1/(2*d)) # Definind the range of X axis where we have unit cell ( a/2 x < -a/2  where a = 7)
Xgradient = np.real(np.gradient(states[:,X_Neg_Idx : X_Pos_Idx,0], dx, axis=1)) #Using the finite difference gradient function of numpy. Since we are solving for the periodic part of the Bloch wave func, any observable has to be integrated along ONLY one unit cell in the periodic dimention (X axis in our case is periodic while Z is not)
np.sum(np.real(-states[:,X_Neg_Idx : X_Pos_Idx,286]) * Xgradient)*dx*dz

