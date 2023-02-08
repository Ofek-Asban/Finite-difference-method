import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import eigs
from scipy.sparse import kronsum
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm


#Global parameteres of the potential. Parameteres based on psudopotential calculation of Helium-4 atom (Takada and Kohn's (1987))
A = 1.827; B = 2.595; C = 1.42; alpha = 1.38319; b0 = 2.52; 
#Fit to Kestner and Rice (1965) potential of He atom with parameters given by A = 2.486; B = 1.99; C = 1.45; alpha = 1.38319; b0 = 4.461; 


#Building a 3D lattice of atoms. Each atom is represented by the pseudo-potential V = a*(b - 1/r)* np.exp(-c*r) - alpha/(b0 + 2r**4)  where r = x**2 + y**2 + z**2
def V_3D_Lattice(xmin, xmax, ymin, ymax, zmin, zmax, dx ,dy, dz, a, X, Y, Z):
    I = np.arange(xmin * a, xmax * a, a); J= np.arange(ymin * a, ymax * a, a); K = np.arange(zmin * a, zmax * a, a);
    i = j = k = V_Latt = 0;
    for i in I:
        for j in J:
            for k in K:
                r = np.sqrt((X - i)**2 + (Y - j)**2 + (Z - k)**2) #Defining position vector on lattice (I,J,K) 
                V_Latt = V_Latt + A*(B - 1/(r + 0.36))* np.exp(-C * r) #The single atom potential as a finctin of r
    V = np.heaviside(-Z + a/2,0)*V_Latt + np.heaviside(Z - a/2,0)*(0.007/Z) #total potential of electron on Helium (He is latice at z<0 and and image potential -0.007/z for z>0)
    return V


#################Finite difference for 3D problem (verified for 3D potential well with infinit walls or PBC) ###############

def H_3D(xmin, xmax, ymin, ymax, zmin, zmax, a, MaxIdx, P):  
    dx = dy = dz = 0.2*a; #Discritization step of space (5 points per unit cell a = 7 Bohr)
    #Descretizing space in the x direction
    x = np.arange(xmin * a, xmax * a, dx); y = np.arange(ymin * a, ymax * a, dy); z = np.arange(zmin * a, zmax * a, dz);
    X, Z, Y = np.meshgrid(x, y, z) #Making 3D space from 3 arrays
    V = V_3D_Lattice(xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, dz, a, X, Y, Z) #Calls a funtion that creates a lattice of atoms with realistic pseudo-potential
    Dx = diags([1, -2, 1, P, P], [-1, 0, 1, x.size-1, 1-x.size], shape=(x.size,x.size))/(dx**2); # descrete form of 2nd derivative in real space
    Dy = diags([1, -2, 1, P, P], [-1, 0, 1, y.size-1, 1-y.size], shape=(y.size,y.size))/(dy**2);
    Dz = diags([1, -2, 1], [-1, 0, 1], shape=(z.size,z.size))/(dz**2);
    # In the scheme of making the wave funciton 2D matrix into a vector, the 2D derivative is Dx * 1 + 1 * Dy = Dx + Dy (simple generalization to 3D is Dx + Dy + Dz)
    T = - 0.5 * kronsum(Dx,kronsum(Dy,Dz)) #kron sum of 3 sparce matrices
    H = T + diags(V.reshape(-1))
    energies, states = eigs(H, k = MaxIdx, which='SM') # Calculationg the MaxIdx lowest eigenstates (which='SM' meaning smallest)
    idx = np.argsort(energies) #array of index values sorted according to increasing value
    energies = energies[idx] #sorting eigenvalue array
    states = states[:,idx] #sorting eigenvactors according of increasing value of their eigenvalues
    return X, Y, Z, energies, states, V, MaxIdx

X, Y, Z, energies3D, states3D, V3D, MaxIdx3D = H_3D(-1, 1, -1, 1, -2, 100 , 7, 300 , 1) # xmin = zmin = -1, and xmax = zmax = 1 since the problem is periodic in x,y directions it is enough to take one units cell   


315775*np.real(energies3D) #Energy values in Kelvin (315775 is conversion coefficient from Hartries to Kelvins)

#plot energies (eigenvalues) of the 3D Schrodinger equation
plt.plot(np.arange(0, MaxIdx3D, 1), 315775*np.real(energies3D), marker='o') #Energy (in Kelvin) as a func level number 


#Plot 3D wave function. The color represented by red/blue spectrum as high/low values of the potential
plt.axes(projection='3d').scatter(X, Y, Z, c=states3D[:,1].reshape(X.shape), cmap= cm.coolwarm, alpha = 0.1)

############################################################################################################################

