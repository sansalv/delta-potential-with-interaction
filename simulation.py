from re import X
import numpy as np
import matplotlib.pyplot as plt
import os

# Parameters for the simulation

N      = 81                      # Number of gridpoints on the 1D grid
N_part = 2                       # Number of particles
L      = 1.0                     # Length of the box
h      = L/N                     # Step size
Id     = np.diagflat(np.ones(N)) # Identity matrix

def kinetic_energy_matrix(index):
    """
    Computes the kinetic energy matrix for particle number "index". The kinetic energy matrix
    is always the same, this also tensors that matrix with identity.

    "Index" starts from 0.

    """
    # Diagflat creates diagonal matrices from an array. The second argument is how "far away" from the main diagonal the array is put.

    ke = np.diagflat(-2.0 * np.ones(N)) + np.diagflat(np.ones(N-1),-1) + np.diagflat(np.ones(N-1),1) # Single-particle finite difference matrix.
    ke = -0.5*(1/h**2)*ke # In the Schrödinger equation, we have -1/2 * d^2/dx^2, so we multiply with 0.5

    H_kin = np.zeros((N,N))
    if index == 0: # If we're dealing with the first particle, the initial matrix should be the kinetic energy matrix. Otherwise the identity.
        H_kin = ke
    else:
        H_kin = Id 
    for i in range(1, N_part): # Loop over all the particles, tensoring 1-particle identity or 1-particle kinetic energy in to the kinetic energy matrix.
        if i == index:
            H_kin = np.kron(H_kin, ke) # "Kron" is the finite-dimensional tensor product.
        else:
            H_kin = np.kron(H_kin, Id)

    return H_kin

def position_matrix(index):
    """
    Computes the operator x. We are using the position basis, so this operator is in fact
    diagonal. We use the same tensoring trick as above, and "index" has once again the same meaning.
    """
    xs = np.linspace(0,L,N) # Creates an array with values (0.0,0.01, ... 0.99). Divides L in to N pieces, in other words.
    x = np.diagflat(xs)

    H_x = np.zeros((N,N))

    if index == 0: # Same tricks as above
        H_x = x
    else:
        H_x = Id

    for i in range(1,N_part):
        if i == index:
            H_x = np.kron(H_x, x)
        else:
            H_x = np.kron(H_x, Id)
    
    return H_x


def harmonic_interaction_matrix(x1,x2,k):
    """
    Harmonic interaction matrix. This function only works for two-particle systems for now. The potential is 0.5*k*(x1-x2)^2.
    """
    assert N_part == 2
    k = k/h**2

    return 0.5*k*(x1-x2)@(x1-x2) # Recall that @ is matrix multiplication.

def dirac_delta_matrix(index, loc, v):
    """
    Dirac delta matrix. In the discrete approximation, the Dirac delta is 1/h in one spot on the diagonal, 0 elsewhere. 

    Index is again the particle index. Loc is the location in the array for the 1/h term.
    """

    diracs      = np.zeros(N)
    diracs[loc] = v*1/h

    dirac = np.diagflat(diracs)

    H_dirac = np.zeros((N,N))
    if index == 0:
        H_dirac = dirac
    else:
        H_dirac = Id

    for i in range(1,N_part):
        if i == index:
            H_dirac = np.kron(H_dirac, dirac)
        else:
            H_dirac = np.kron(H_dirac, Id)

    return H_dirac


def save_plots_with_infotxt(path, e, vec, P, v1, v2, k):
	"""
	Save plots and create info.txt file. 
	
	Variable "path" is the chosen directory path for the plots and info.txt file.
	"""
	txtfile = open(f"{path}/info.txt", "a") # Create and open info.txt file
	txtfile.truncate(0) # Empty the file from previous
	
	# Print the system info first
	print(f"N = {N}", file=txtfile)
	print(f"h = {np.round(h,3)}", file=txtfile)
	print(f"v1 = {v1}", file=txtfile)
	print(f"v2 = {v2}", file=txtfile)
	print(f"k = {k}/h**2", file=txtfile)
	print("----------------------------------", file=txtfile)
	
	# Difference parameter "diff" is used to group degenerate energy states in the txtfile
	# If Dirac's delta or interaction is nonzero, we switch to group approximate degenerate states (for example diff = 3)
	diff = 0.00002
	if abs(v1)+abs(v2)+abs(k) != 0:
		diff = 3
	
	# Loop for saving plots into the directory
	# Also, listing one particle mode energies into "modes" list
	modes = []
	for i in range(P):
		vec_grid = vec[:,i].reshape(N,N)
		# np.abs(vec_grid)**2
		plt.imshow(vec_grid, cmap="RdBu", origin="lower")
		plt.savefig(f"{path}/plot{i+1}.png")
		
		if (abs(e[i]-e[i-1])>diff and abs(e[i]-e[i+1])>diff) or i==0:
			modes.append(e[i]/2)
		
	print_vectors_with_modes(e, P, modes, diff, txtfile)
	txtfile.close()
	print(open(f"{path}/info.txt", "r").read()) # Print the finished txtfile also into the terminal

def print_vectors_with_modes(e, P, modes, diff, txtfile):
	e1 = np.round(e,1)
	print(f"Mode energies:\n{np.round(modes,1)}", file=txtfile)
	print("----------------------------------", file=txtfile)
	for i in range(0,P):
		pair_found = False # These loops find the energy sums of the modes
		for j1 in range(len(modes)):
			for j2 in range(len(modes)):
				if abs(modes[j1]+modes[j2]-e[i]) < diff: # energy = mode1 + mode2
					pair_found = True
					if abs(e[i]-e[i-1])>diff and i!=0:
						print("----------------------------------", file=txtfile)
					print(f"Vector #{i+1}: Energy = {e1[i]} = ({j1+1},{j2+1})", file=txtfile)
					break
			if pair_found:
				break
				
def main():
    # Main simulation. Create a system, use above functions, plot something.
	P = int(input("Number of eigenvalue plots: "))
	print("Forming Hamiltonian.. ")
	H = kinetic_energy_matrix(0) + kinetic_energy_matrix(1) # Forming the full matrix as the sum of the above terms.

	v1 = float(input("Particle 1 Dirac's delta potential strength:\nv1 =  "))
	if v1 != 0:
		H += dirac_delta_matrix(0,(N-1)//2, v1)
	v2 = float(input("Particle 2 Dirac's delta potential strength:\nv2 =  "))
	if v2 != 0:
		H += dirac_delta_matrix(1,(N-1)//2, v2)
	
	k = float(input("Harmonic interaction strength:\nk =  "))
	if k != 0:
		H += harmonic_interaction_matrix(position_matrix(0), position_matrix(1), k)

	print("Solving the eigenproblem.. ")
	e, vec = np.linalg.eigh(H) # vec stores the eigenvectors. vec[:,0] is the lowest energy, v[:,1] second lowest, etc.
	print("----------------------------------")
	
	save_plots_with_infotxt("plots/test_with_delta_and_interaction", e, vec, P, v1, v2, k) # Save plots and create info.txt file
	
if __name__ == "__main__":
    main()

	
	

