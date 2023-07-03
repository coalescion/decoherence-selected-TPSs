###################################################################################################################
###
###  Time Evolution 
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###  The main function of this module is to time evolve a reduced density matrix 
###  
####################################################################################################################


from numpy.linalg import norm
from qutip import *


# =============================== CHARACTERISTIC TIME ===============================

def characteristic_time(Hamiltonian):
    '''Calculates the characteristic time for a given Hamiltonian.
    
        Parameters: 
        Hamiltonian: The hamiltonian matrix.
        
        Returns: The characteristic time.'''

    #Time inversely related to energy
    return 1/(norm(Hamiltonian, 2))
    

# =============================== TIME EVOLUTION ===============================

def time_evolution(initial_state, H, t):
    '''Time evolution on density matrix.
    
        Parameters:
        initial state: State to evolve.
        H: Hamiltonian to use in evolution.
        t: Time to evolve to.
        
        Returns: Evolved state.'''
    
    te_state = ((1j)*H*t).expm() * initial_state * (-(1j)*H*t).expm()
        
    return te_state