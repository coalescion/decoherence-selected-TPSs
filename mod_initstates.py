###################################################################################################################
###
###  Creating Initial States
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###  The main function of this module is to create the initial states of the env and the sys  
###  -- takes in two variables:
###  ----> n: number of coupled qubits
###  -- outputs collective initial state 
###
####################################################################################################################


import numpy as np 
from qutip import *

def init_states(n):
    '''Parameters:
        n: Number of environmental Qubits.

        Returns: Initial state to be evolved forward in time.'''

    t_env_i = 1/(2**n)*tensor([qeye(2)]*n)
    env_i = (1/(2**n))*np.identity(2**n)

    tensor_initial_states = []
    initial_states = []
    for i in range(2):
        sys_i = np.zeros((2, 2), int)
        sys_i[i][i] = 1
        sys_i = Qobj(sys_i)
        tensor_state = tensor(sys_i, t_env_i)
        tensor_initial_states.append(tensor_state)
        
        sys_i = np.zeros((2, 2), int)
        sys_i[i][i] = 1
        state = Qobj(np.kron(sys_i, env_i))
        initial_states.append(state)

    return tensor_initial_states, initial_states
