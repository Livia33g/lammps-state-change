#!/usr/bin/env python3
"""
Analyze state changes and verify they only occur when dimers form.
This script:
1. Parses LAMMPS trajectory to track atom types (state changes)
2. Calculates COM distances between molecules
3. Identifies dimers (COM distance < 3.0)
4. Correlates state changes with dimer formation
"""

import sys
import numpy as np
from collections import defaultdict

def parse_trajectory(filename):
    """Parse LAMMPS trajectory file and extract atom data at each timestep."""
    timesteps = {}
    current_timestep = None
    atoms = []
    box_size = None
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        
        # Check for timestep
        if line == "ITEM: TIMESTEP":
            # Save previous timestep if exists
            if current_timestep is not None and atoms:
                timesteps[current_timestep] = {
                    'atoms': atoms.copy(),
                    'box_size': box_size
                }
                atoms = []
            # Read next line for timestep value
            i += 1
            if i < len(lines):
                current_timestep = int(lines[i].strip())
            i += 1
            continue
        
        # Check for box bounds
        if line.startswith("ITEM: BOX BOUNDS"):
            # Read box dimensions
            bounds = []
            for j in range(3):
                i += 1
                if i < len(lines):
                    bounds_line = lines[i].strip().split()
                    if len(bounds_line) >= 2:
                        bounds.append((float(bounds_line[0]), float(bounds_line[1])))
            if len(bounds) == 3:
                box_size = bounds[0][1] - bounds[0][0]  # Assume cubic box
            i += 1
            continue
        
        # Check for atom data
        if line == "ITEM: ATOMS id mol type x y z":
            i += 1
            # Read atoms until next ITEM or end of file
            while i < len(lines):
                atom_line = lines[i].strip()
                if not atom_line:
                    i += 1
                    continue
                if atom_line.startswith("ITEM:"):
                    # Process this line in next iteration
                    break
                parts = atom_line.split()
                if len(parts) >= 6:
                    atom_id = int(parts[0])
                    mol_id = int(parts[1])
                    atom_type = int(parts[2])
                    x, y, z = float(parts[3]), float(parts[4]), float(parts[5])
                    atoms.append({
                        'id': atom_id,
                        'mol': mol_id,
                        'type': atom_type,
                        'x': x,
                        'y': y,
                        'z': z
                    })
                i += 1
            continue
        
        i += 1
    
    # Save last timestep
    if current_timestep is not None and atoms:
        timesteps[current_timestep] = {
            'atoms': atoms.copy(),
            'box_size': box_size
        }
    
    return timesteps

def calculate_com(atoms, mol_id, box_size):
    """Calculate center of mass for a molecule."""
    mol_atoms = [a for a in atoms if a['mol'] == mol_id]
    if not mol_atoms:
        return None
    
    com_x = sum(a['x'] for a in mol_atoms) / len(mol_atoms)
    com_y = sum(a['y'] for a in mol_atoms) / len(mol_atoms)
    com_z = sum(a['z'] for a in mol_atoms) / len(mol_atoms)
    
    return np.array([com_x, com_y, com_z])

def pbc_distance(r1, r2, box_size):
    """Calculate minimum image distance with periodic boundary conditions."""
    dr = r2 - r1
    dr -= box_size * np.round(dr / box_size)
    return np.linalg.norm(dr)

def find_dimers(atoms, box_size, dimer_cutoff=3.0):
    """Find dimers based on COM distance."""
    # Get unique molecule IDs
    mol_ids = sorted(set(a['mol'] for a in atoms))
    
    # Calculate COM for each molecule
    mol_coms = {}
    for mol_id in mol_ids:
        com = calculate_com(atoms, mol_id, box_size)
        if com is not None:
            mol_coms[mol_id] = com
    
    # Find dimers (pairs with COM distance < cutoff)
    dimers = []
    for i, mol1_id in enumerate(mol_ids):
        if mol1_id not in mol_coms:
            continue
        for mol2_id in mol_ids[i+1:]:
            if mol2_id not in mol_coms:
                continue
            dist = pbc_distance(mol_coms[mol1_id], mol_coms[mol2_id], box_size)
            if dist < dimer_cutoff:
                dimers.append((mol1_id, mol2_id, dist))
    
    return dimers, mol_coms

def get_patch_types(atoms, mol_id):
    """Get patch types for a molecule."""
    patches = [a for a in atoms if a['mol'] == mol_id and a['type'] in [2, 3]]
    if not patches:
        return []
    # All patches in a molecule should have the same type (they change together)
    return [p['type'] for p in patches]

def analyze_state_changes_with_dimers(timesteps, dimer_cutoff=3.0):
    """Analyze state changes and correlate with dimer formation."""
    sorted_timesteps = sorted(timesteps.keys())
    
    print(f"Analyzing {len(sorted_timesteps)} timesteps")
    print(f"Dimer cutoff: {dimer_cutoff}")
    print("=" * 70)
    print()
    
    # Track state changes
    state_changes = []  # List of (timestep, mol_id, old_type, new_type, is_dimer)
    
    # Track dimers at each timestep
    dimers_by_timestep = {}
    
    # Compare consecutive timesteps
    for i in range(len(sorted_timesteps) - 1):
        t1 = sorted_timesteps[i]
        t2 = sorted_timesteps[i + 1]
        
        data1 = timesteps[t1]
        data2 = timesteps[t2]
        atoms1 = data1['atoms']
        atoms2 = data2['atoms']
        box_size = data1['box_size'] or data2['box_size']
        
        # Find dimers at timestep t2 (after potential state change)
        dimers, mol_coms = find_dimers(atoms2, box_size, dimer_cutoff)
        dimers_by_timestep[t2] = dimers
        
        # Create sets of dimers for quick lookup
        dimer_molecules = set()
        for mol1, mol2, dist in dimers:
            dimer_molecules.add(mol1)
            dimer_molecules.add(mol2)
        
        # Find molecules that changed type
        mol_types1 = {}
        mol_types2 = {}
        
        for mol_id in set(a['mol'] for a in atoms1):
            patch_types1 = get_patch_types(atoms1, mol_id)
            patch_types2 = get_patch_types(atoms2, mol_id)
            
            if patch_types1 and patch_types2:
                # Use the first patch type (all patches in molecule have same type)
                mol_types1[mol_id] = patch_types1[0]
                mol_types2[mol_id] = patch_types2[0]
        
        # Detect state changes
        for mol_id in mol_types1:
            if mol_id not in mol_types2:
                continue
            
            type1 = mol_types1[mol_id]
            type2 = mol_types2[mol_id]
            
            # Only track changes between type 2 and 3 (patch types)
            if type1 != type2 and type1 in [2, 3] and type2 in [2, 3]:
                is_in_dimer = mol_id in dimer_molecules
                state_changes.append({
                    'timestep': t2,
                    'mol_id': mol_id,
                    'old_type': type1,
                    'new_type': type2,
                    'is_dimer': is_in_dimer,
                    'dimer_partner': None
                })
                
                # Find dimer partner if in dimer
                if is_in_dimer:
                    for mol1, mol2, dist in dimers:
                        if mol1 == mol_id:
                            state_changes[-1]['dimer_partner'] = mol2
                            state_changes[-1]['dimer_distance'] = dist
                            break
                        elif mol2 == mol_id:
                            state_changes[-1]['dimer_partner'] = mol1
                            state_changes[-1]['dimer_distance'] = dist
                            break
    
    # Analysis
    print("=" * 70)
    print("STATE CHANGE ANALYSIS")
    print("=" * 70)
    print()
    
    total_changes = len(state_changes)
    changes_in_dimers = sum(1 for sc in state_changes if sc['is_dimer'])
    changes_not_in_dimers = total_changes - changes_in_dimers
    
    print(f"Total state changes detected: {total_changes}")
    print(f"  • State changes while in dimer: {changes_in_dimers} ({100*changes_in_dimers/max(total_changes,1):.1f}%)")
    print(f"  • State changes while NOT in dimer: {changes_not_in_dimers} ({100*changes_not_in_dimers/max(total_changes,1):.1f}%)")
    print()
    
    if total_changes > 0:
        print("✅ VERIFICATION:")
        if changes_not_in_dimers == 0:
            print("   ✅ ALL state changes occurred while molecules were in dimers!")
            print("   ✅ State change mechanism is working correctly!")
        else:
            print(f"   ⚠️  {changes_not_in_dimers} state changes occurred while NOT in dimers")
            print("   ⚠️  This may indicate the coordination cutoff is too loose")
            print()
            print("   State changes NOT in dimers:")
            for sc in state_changes:
                if not sc['is_dimer']:
                    print(f"      Timestep {sc['timestep']}: Mol {sc['mol_id']} "
                          f"({sc['old_type']} -> {sc['new_type']})")
    
    print()
    print("=" * 70)
    print("DETAILED STATE CHANGES")
    print("=" * 70)
    print()
    
    if state_changes:
        print(f"{'Timestep':<12} {'Mol ID':<8} {'Type':<12} {'In Dimer':<10} {'Partner':<10} {'Distance':<10}")
        print("-" * 70)
        for sc in state_changes[:20]:  # Show first 20
            partner_str = str(sc['dimer_partner']) if sc['dimer_partner'] else "N/A"
            dist_str = f"{sc['dimer_distance']:.3f}" if 'dimer_distance' in sc else "N/A"
            dimer_str = "YES" if sc['is_dimer'] else "NO"
            print(f"{sc['timestep']:<12} {sc['mol_id']:<8} "
                  f"{sc['old_type']}->{sc['new_type']:<8} {dimer_str:<10} "
                  f"{partner_str:<10} {dist_str:<10}")
        if len(state_changes) > 20:
            print(f"... and {len(state_changes) - 20} more")
    else:
        print("No state changes detected in this trajectory.")
    
    print()
    print("=" * 70)
    print("DIMER STATISTICS")
    print("=" * 70)
    print()
    
    # Count dimers over time
    dimer_counts = [len(dimers_by_timestep.get(t, [])) for t in sorted_timesteps]
    if dimer_counts:
        avg_dimers = np.mean(dimer_counts)
        max_dimers = max(dimer_counts)
        print(f"Average number of dimers: {avg_dimers:.2f}")
        print(f"Maximum number of dimers: {max_dimers}")
        print(f"Timesteps with dimers: {sum(1 for d in dimer_counts if d > 0)} / {len(dimer_counts)}")
    
    return state_changes, dimers_by_timestep

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 analyze_state_changes_with_dimers.py <trajectory_file> [dimer_cutoff]")
        print("  trajectory_file: LAMMPS dump file")
        print("  dimer_cutoff: COM distance threshold for dimers (default: 3.0)")
        sys.exit(1)
    
    trajectory_file = sys.argv[1]
    dimer_cutoff = float(sys.argv[2]) if len(sys.argv) > 2 else 3.0
    
    print(f"Analyzing trajectory: {trajectory_file}")
    print(f"Dimer cutoff: {dimer_cutoff}")
    print()
    
    timesteps = parse_trajectory(trajectory_file)
    state_changes, dimers = analyze_state_changes_with_dimers(timesteps, dimer_cutoff)

