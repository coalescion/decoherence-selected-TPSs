###################################################################################################################
###
###  Scoring
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###  The main function of this module is to score the different tensor factrizations of the Hilbert Space  
###
####################################################################################################################


import math
import numpy as np 
import mod_GGMM
from qutip import *

def construct_unitary(nVar, thetas, GGMMs):
    '''Constructs unitary matrix from thetas.
    
        Parameters:
        thetas: Theta coefficients.
        
        Returns: Unitary matrix.'''
        
    param_GGMM = []
    for i in range(len(thetas)):
        param_GGMM.append(Qobj(thetas[i]*GGMMs[i].full(), dims=[[2]*nVar,[2]*nVar]))
    
    sumG = Qobj(sum(param_GGMM))
    scrambler = (-1j*sumG).expm()
    return scrambler 

def construct_localUnitary(n, thetas):
    '''
        Parameters:
        thetas: List of 3 element lists. Each 3 element list corresponds to one qubit's SU(2) factor.
        
        Returns: Unitary matrix.
        '''
    
    allLambdas = []
    selfLambdas = mod_GGMM.construct_ggmm_sub(2)
    for i in range(n):
        allLambdas += selfLambdas

    sumList = []
    tensorList = []
    a = 0
    for i in range(n):
        for j in range(a, a+3):
            factor = thetas[j]*allLambdas[j]
            sumList.append(factor)
            a += 1
        factor = Qobj(-1j*sum(sumList)).expm()
        tensorList.append(factor)
    return tensor(tensorList)