"""
Example submission for problem-001-ksat-advanced.

This demonstrates ONE possible encoding/design/decode strategy. You are free
to use a completely different approach.

STRATEGY OVERVIEW:
- Variables → particles with two states (A=unset, C=set to True)
- Clauses → "catalyst" particles (B) that bind to unsatisfied clauses
- State-change policy: When a variable particle touches a clause catalyst,
  flip the variable to satisfy that clause
- Decode: Count final states of variable particles

This is a simplified example. A real solution would need more sophisticated
logic to handle conflicts and optimize work.
"""

import os
import json
from pathlib import Path

API_VERSION = "v1"


class StateChangeSolution:
    """
    Example implementation for abstract 3-SAT problem.
    
    WARNING: This is a minimal working example. Real solutions should:
    - Handle clause conflicts more intelligently
    - Optimize work (minimize unnecessary flips)
    - Use more sophisticated decoding (e.g., check clause satisfaction)
    """

    def __init__(self, problem_def, work_dir):
        self.problem_def = problem_def
        self.work_dir = Path(work_dir)
        
        # Extract SAT instance
        self.n_vars = problem_def['sat_instance']['n_variables']
        self.clauses = problem_def['sat_instance']['clauses']
        self.n_clauses = len(self.clauses)
        
        # Create working directories
        self.sim_dir = self.work_dir / 'simulation'
        self.gen_dir = self.work_dir / 'generated'
        self.sim_dir.mkdir(parents=True, exist_ok=True)
        self.gen_dir.mkdir(parents=True, exist_ok=True)

    def encode(self):
        """
        Encode 3-SAT as patchy particles:
        - 6 variable particles (types 2-7): can be A (unset) or C (True)
        - 5 clause catalyst particles (type 8): bind to unsatisfied clauses
        - 1 core per particle (type 1): excluded volume
        """
        n_vars = self.n_vars
        n_clauses = self.n_clauses
        
        # Total particles: n_vars (variables) + n_clauses (clauses) = 11
        # Each particle: 1 core + 1 patch = 2 atoms
        n_particles = n_vars + n_clauses
        n_atoms = n_particles * 2
        
        # Write LAMMPS data file
        data_file = self.sim_dir / 'data.ksat'
        with open(data_file, 'w') as f:
            f.write(f"# 3-SAT encoding: {n_vars} variables, {n_clauses} clauses\n\n")
            f.write(f"{n_atoms} atoms\n")
            f.write(f"{n_particles} bonds\n\n")
            f.write("2 atom types  # 1=core, 2-7=var patches, 8=clause patches\n")
            f.write("1 bond types\n\n")
            f.write("-10.0 10.0 xlo xhi\n")
            f.write("-10.0 10.0 ylo yhi\n")
            f.write("-10.0 10.0 zlo zhi\n\n")
            f.write("Atoms # full\n\n")
            
            atom_id = 1
            # Variable particles (types 2-7)
            for i in range(n_vars):
                # Core
                f.write(f"{atom_id} 1 {i+1} 0.0 {i*2.0} 0.0 0.0\n")
                atom_id += 1
                # Patch (initially A=type 2)
                f.write(f"{atom_id} 2 {i+1} 0.0 {i*2.0} 0.5 0.0\n")
                atom_id += 1
            
            # Clause particles (type 8)
            for i in range(n_clauses):
                # Core
                f.write(f"{atom_id} 1 {n_vars+i+1} 0.0 {n_vars*2.0 + i*2.0} 0.0 0.0\n")
                atom_id += 1
                # Patch
                f.write(f"{atom_id} 8 {n_vars+i+1} 0.0 {n_vars*2.0 + i*2.0} 0.5 0.0\n")
                atom_id += 1
            
            f.write("\nBonds\n\n")
            bond_id = 1
            for mol_id in range(1, n_particles + 1):
                f.write(f"{bond_id} 1 {mol_id*2-1} {mol_id*2}\n")
                bond_id += 1
        
        # Write LAMMPS input script
        input_file = self.sim_dir / 'in.ksat'
        with open(input_file, 'w') as f:
            f.write("""# 3-SAT state-change simulation

units lj
atom_style full
bond_style harmonic

read_data data.ksat

# Pair interactions
pair_style lj/cut 2.5
pair_coeff 1 1 1.0 1.0 2.5  # core-core repulsion
pair_coeff 1 2 0.0 1.0 2.5  # core-var neutral
pair_coeff 1 8 0.0 1.0 2.5  # core-clause neutral
pair_coeff 2 8 1.0 0.5 0.7  # var-clause attraction (when clause unsatisfied)

# State-change fix
fix sc all state/change/ksat

# Dynamics
velocity all create 0.5 12345
fix nvt all nvt temp 0.5 0.5 0.1

# Output
thermo_style custom step temp pe ke etotal
thermo 1000
dump 1 all custom 1000 dump.ksat.lammpstrj id type x y z

run 100000
""")
        
        return {
            "n_particles": n_particles,
            "n_variables": n_vars,
            "n_clauses": n_clauses,
            "encoding_notes": "Variables as particles with A/C states, clauses as catalysts"
        }

    def design_policy(self, system_meta):
        """
        Generate a simple state-change fix that flips variables when they
        touch clause catalysts.
        """
        cpp_file = self.gen_dir / 'fix_state_change_ksat.cpp'
        h_file = self.gen_dir / 'fix_state_change_ksat.h'
        
        # Write header
        with open(h_file, 'w') as f:
            f.write("""#ifndef LMP_FIX_STATE_CHANGE_KSAT_H
#define LMP_FIX_STATE_CHANGE_KSAT_H

#include "fix.h"

namespace LAMMPS_NS {

class FixStateChangeKSAT : public Fix {
 public:
  FixStateChangeKSAT(class LAMMPS *, int, char **);
  ~FixStateChangeKSAT() {}
  int setmask();
  void init();
  void post_force(int);

 private:
  double cutoff;
  int var_type;  // type 2 (variable patches)
  int clause_type;  // type 8 (clause catalysts)
};

}

#endif
""")
        
        # Write implementation
        with open(cpp_file, 'w') as f:
            f.write("""#include "fix_state_change_ksat.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "error.h"
#include <cmath>

using namespace LAMMPS_NS;

FixStateChangeKSAT::FixStateChangeKSAT(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR, "fix state/change/ksat requires cutoff");
  cutoff = utils::numeric(FLERR, arg[3], false, lmp);
  var_type = 2;
  clause_type = 8;
}

int FixStateChangeKSAT::setmask() {
  return POST_FORCE;
}

void FixStateChangeKSAT::init() {
  // Check that types exist
}

void FixStateChangeKSAT::post_force(int vflag) {
  double **x = atom->x;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  
  // Simple logic: if a variable patch (type 2) is near a clause patch (type 8),
  // flip the variable to type 4 (C = True)
  for (int i = 0; i < nlocal; i++) {
    if (type[i] != var_type) continue;
    
    for (int j = 0; j < nlocal; j++) {
      if (type[j] != clause_type) continue;
      
      double dx = x[i][0] - x[j][0];
      double dy = x[i][1] - x[j][1];
      double dz = x[i][2] - x[j][2];
      double rsq = dx*dx + dy*dy + dz*dz;
      
      if (rsq < cutoff*cutoff) {
        // Flip variable to satisfied state
        type[i] = 4;  // C = True
        break;
      }
    }
  }
}
""")
        
        return {
            "fix_files": [
                "fix_state_change_ksat.cpp",
                "fix_state_change_ksat.h"
            ]
        }

    def decode(self):
        """
        Decode assignment by reading final particle states from trajectory.
        This is a simplified example - real solutions should verify clause satisfaction.
        """
        # In a real implementation, you would:
        # 1. Read dump.ksat.lammpstrj
        # 2. Identify variable particles (molecules 1-6)
        # 3. Check their final patch types (2=A/unset, 4=C/True)
        # 4. Map to truth assignment
        
        # For this example, return a placeholder assignment
        # Real code would parse the trajectory file
        assignment = [1] * self.n_vars  # Placeholder: all True
        
        # Verify this assignment satisfies all clauses
        satisfied = True
        for clause in self.clauses:
            clause_sat = False
            for lit in clause:
                var_idx = abs(lit) - 1  # Convert to 0-indexed
                if lit > 0:
                    clause_sat = clause_sat or (assignment[var_idx] == 1)
                else:
                    clause_sat = clause_sat or (assignment[var_idx] == 0)
            if not clause_sat:
                satisfied = False
                break
        
        return {
            "assignment": assignment,
            "satisfied": satisfied,
            "diagnostics": {
                "note": "This is a minimal example. Real solutions must parse trajectory files."
            }
        }

