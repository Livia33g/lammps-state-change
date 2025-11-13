#!/usr/bin/env python3
"""
Script to analyze LAMMPS trajectory and detect state changes (type changes).
"""

import sys
import re
from collections import defaultdict

def parse_trajectory(filename):
    """Parse LAMMPS trajectory file and extract atom types at each timestep."""
    timesteps = {}
    current_timestep = None
    atoms = []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            
            # Check for timestep
            if line == "ITEM: TIMESTEP":
                # Save previous timestep if exists
                if current_timestep is not None and atoms:
                    timesteps[current_timestep] = atoms.copy()
                    atoms = []
                # Read next line for timestep value
                current_timestep = int(f.readline().strip())
                continue
            
            # Check for atom data
            if line == "ITEM: ATOMS id mol type x y z":
                # Read atoms until next ITEM or end of file
                for atom_line in f:
                    atom_line = atom_line.strip()
                    if not atom_line:
                        continue
                    if atom_line.startswith("ITEM:"):
                        # Put this line back for next iteration
                        f.seek(f.tell() - len(atom_line.encode()) - 1)
                        break
                    parts = atom_line.split()
                    if len(parts) >= 3:
                        atom_id = int(parts[0])
                        mol_id = int(parts[1])
                        atom_type = int(parts[2])
                        atoms.append({
                            'id': atom_id,
                            'mol': mol_id,
                            'type': atom_type
                        })
        
        # Save last timestep
        if current_timestep is not None and atoms:
            timesteps[current_timestep] = atoms.copy()
    
    return timesteps

def analyze_state_changes(timesteps):
    """Analyze state changes by comparing atom types across timesteps."""
    if len(timesteps) < 2:
        print("Need at least 2 timesteps to detect changes")
        return
    
    # Get all timesteps sorted
    sorted_timesteps = sorted(timesteps.keys())
    
    print(f"Found {len(sorted_timesteps)} timesteps: {sorted_timesteps[:10]}...")
    print()
    
    # Track type changes for each atom
    atom_changes = defaultdict(list)
    
    # Compare consecutive timesteps
    for i in range(len(sorted_timesteps) - 1):
        t1 = sorted_timesteps[i]
        t2 = sorted_timesteps[i + 1]
        
        atoms1 = {a['id']: a for a in timesteps[t1]}
        atoms2 = {a['id']: a for a in timesteps[t2]}
        
        # Find atoms that changed type
        changes = []
        for atom_id in atoms1:
            if atom_id in atoms2:
                type1 = atoms1[atom_id]['type']
                type2 = atoms2[atom_id]['type']
                mol_id = atoms1[atom_id]['mol']
                
                if type1 != type2:
                    changes.append({
                        'atom_id': atom_id,
                        'mol_id': mol_id,
                        'old_type': type1,
                        'new_type': type2,
                        'from_step': t1,
                        'to_step': t2
                    })
                    atom_changes[atom_id].append((t1, t2, type1, type2))
        
        if changes:
            print(f"State changes detected between timestep {t1} and {t2}:")
            for change in changes:
                print(f"  Atom {change['atom_id']} (mol {change['mol_id']}): "
                      f"type {change['old_type']} -> {change['new_type']}")
            print()
    
    # Summary statistics
    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)
    
    total_changes = sum(len(changes) for changes in atom_changes.values())
    atoms_with_changes = len(atom_changes)
    
    print(f"Total number of state changes: {total_changes}")
    print(f"Number of atoms that changed type: {atoms_with_changes}")
    print(f"Total number of atoms: {len(timesteps[sorted_timesteps[0]])}")
    
    # Count changes by type transition
    type_transitions = defaultdict(int)
    for atom_id, changes in atom_changes.items():
        for t1, t2, old_type, new_type in changes:
            type_transitions[(old_type, new_type)] += 1
    
    print("\nType transitions:")
    for (old_type, new_type), count in sorted(type_transitions.items()):
        print(f"  Type {old_type} -> Type {new_type}: {count} times")
    
    # Show atoms that changed most frequently
    if atom_changes:
        print("\nAtoms with most changes:")
        sorted_atoms = sorted(atom_changes.items(), key=lambda x: len(x[1]), reverse=True)
        for atom_id, changes in sorted_atoms[:10]:
            print(f"  Atom {atom_id}: {len(changes)} changes")
            for t1, t2, old_type, new_type in changes[:5]:  # Show first 5
                print(f"    Step {t1}->{t2}: {old_type}->{new_type}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 check_state_changes.py <trajectory_file>")
        sys.exit(1)
    
    trajectory_file = sys.argv[1]
    print(f"Analyzing trajectory file: {trajectory_file}")
    print("=" * 60)
    print()
    
    timesteps = parse_trajectory(trajectory_file)
    analyze_state_changes(timesteps)

