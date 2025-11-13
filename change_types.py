#!/usr/bin/env python3
"""
Python script to change atom types in LAMMPS based on coordination and random probability.
This script is called from LAMMPS using fix python/invoke or similar mechanism.
"""

import sys
import random

def change_types(lmp, nlocal, tag, type_array, coord_A, coord_B, prob, seed):
    """
    Change atom types based on coordination and random probability.
    
    Args:
        lmp: LAMMPS instance
        nlocal: Number of local atoms
        tag: Atom tags
        type_array: Current atom types (will be modified)
        coord_A: Coordination numbers for type 2 atoms
        coord_B: Coordination numbers for type 3 atoms
        prob: Probability of state change
        seed: Random seed
    """
    random.seed(seed)
    
    for i in range(nlocal):
        current_type = int(type_array[i])
        
        # Type 2 atoms (patches_A): change to 3 if coordinated and random < prob
        if current_type == 2:
            if coord_A[i] > 0.1 and random.random() < prob:
                type_array[i] = 3
        
        # Type 3 atoms (patches_B): change to 2 if coordinated and random < prob
        elif current_type == 3:
            if coord_B[i] > 0.1 and random.random() < prob:
                type_array[i] = 2
        
        # Type 1 atoms: no change
        # else: type stays 1

