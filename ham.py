import numpy as np

# Parameters
L = 6  # Number of sites
t1 = 0.5  # Hopping amplitude for a-particles
t2 = 0.5  # Hopping amplitude for b-particles
t12 = 0.5  # Coupling between a- and b-particles

phi=np.pi/16.0
def construct_hamiltonian(L, t1, t2, t12, phi):
    # Hamiltonian size: 2L (L for a and L for b)
    H = np.zeros((2 * L, 2 * L), dtype=complex)
    
    # Hopping terms for a-particles
    for j in range(L - 1):
        H[j, j + 1] = t1*1j *np.exp(1j * phi/L)
        H[j + 1, j] = t1*(-1j)*np.exp(-1j * phi/L)
    # PBC term for a-particles
    H[L - 1, 0] = t1*1j * np.exp(1j * phi/L)
    H[0, L - 1] = t1*(-1j) * np.exp(-1j * phi/L)
    
    # Hopping terms for b-particles
    for j in range(L, 2 * L - 1):
        H[j, j + 1] = -t2*1j * np.exp(1j * phi/L)
        H[j + 1, j] = -t2*(-1j)  * np.exp(-1j * phi/L)
    # PBC term for b-particles
    H[2 * L - 1, L] = -t2*(1j) * np.exp(1j * phi/L)
    H[L, 2 * L - 1] = -t2*(-1j) * np.exp(-1j * phi/L)
    
    # Coupling terms between a and b particles
    for j in range(L - 1):
        H[j, j + L+1] = t12
        H[j + L+1, j] = t12
   
        
    # PBC coupling terms
    H[L - 1, L] = t12
    H[L , L-1] = t12
    
    for j in range(1,L):
        H[j, j + L-1] = t12
        H[j + L-1, j] = t12
          
        
    # PBC coupling terms
    H[0,2*L-1] = t12
    H[2*L-1,0] = t12

    
    return H


H = construct_hamiltonian(L, t1, t2, t12, phi)
energies = np.linalg.eigvalsh(H)  # Compute eigenvalues (sorted)
    

print(energies[0])



