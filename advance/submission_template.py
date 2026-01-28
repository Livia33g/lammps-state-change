"""
Template for advanced state-change competition submissions.

Participants should:
  - Copy this entire file.
  - Fill in the methods of StateChangeSolution.
  - Paste the final script into the competition Google Form.

Do NOT change:
  - API_VERSION
  - The name of the class (StateChangeSolution)
  - The method names or their required arguments
"""

API_VERSION = "v1"


class StateChangeSolution:
    """
    User-implemented solution for an open-design state-change problem.

    The competition backend will:
      1. Load the problem definition from problem.json (as a dict).
      2. Create a dedicated working directory (work_dir) on the cluster.
      3. Instantiate this class: StateChangeSolution(problem_def, work_dir).
      4. Call encode(), then design_policy(), then (after simulation) decode().
    """

    def __init__(self, problem_def, work_dir):
        """
        Parameters
        ----------
        problem_def : dict
            Parsed JSON describing the abstract problem instance
            (e.g. a k-SAT formula). The exact schema will be documented
            in the problem README.

        work_dir : str or Path-like
            Directory that you may write files into. You MUST NOT read
            or write outside this directory.
        """
        self.problem_def = problem_def
        self.work_dir = work_dir

    # ------------------------------------------------------------------
    # PHASE 1: ENCODE
    # ------------------------------------------------------------------
    def encode(self):
        """
        Encode the abstract problem into a concrete LAMMPS system.

        You MUST:
          - Create a directory:
                {work_dir}/simulation/
          - Inside that directory, write at least:
                in.system      # LAMMPS input script
                data.system    # LAMMPS data file
          - You may write additional files as needed.

        Returns
        -------
        dict
            Metadata about the constructed system. This will be passed as
            the 'system_meta' argument to design_policy(). Keep it
            lightweight and JSON-serializable, e.g.:

                {
                  "n_particles": 1234,
                  "variable_mapping": {...},
                  "notes": "optional text"
                }
        """
        raise NotImplementedError("encode() must be implemented by the competitor")

    # ------------------------------------------------------------------
    # PHASE 2: DESIGN POLICY
    # ------------------------------------------------------------------
    def design_policy(self, system_meta):
        """
        Generate the C++ state-change fix code for this encoding.

        You MUST:
          - Create a directory:
                {work_dir}/generated/
          - Inside that directory, write one or more C++ / header files,
            for example:
                fix_state_change_competitor.cpp
                fix_state_change_competitor.h
          - Follow the naming rules given in the competition guide.

        Returns
        -------
        dict
            Must contain a list of C++/header filenames relative to the
            'generated/' directory, e.g.:

                {
                  "fix_files": [
                    "fix_state_change_competitor.cpp",
                    "fix_state_change_competitor.h"
                  ]
                }

            These files will be copied into the LAMMPS source tree and
            used during compilation.
        """
        raise NotImplementedError("design_policy() must be implemented by the competitor")

    # ------------------------------------------------------------------
    # PHASE 3: DECODE
    # ------------------------------------------------------------------
    def decode(self):
        """
        Decode a solution from the LAMMPS trajectory.

        The competition backend will:
          - Compile LAMMPS with your generated fix files.
          - Run a simulation using the files in {work_dir}/simulation/.

        At this point you may:
          - Read LAMMPS output from {work_dir}/simulation/
            (e.g. dump files, log.lammps, custom output).
          - Perform any analysis needed to infer a solution to the
            abstract problem (e.g. a satisfying assignment to a k-SAT
            instance).

        Returns
        -------
        dict
            Must be fully JSON-serializable. A typical structure is:

                {
                  "assignment": [0, 1, 1, 0, ...],
                  "satisfied": true,
                  "diagnostics": {
                    "num_state_changes": 42,
                    "energy_series": [...],
                    "notes": "optional"
                  }
                }

            - 'assignment' should have one entry per variable in the
              problem definition (schema documented with the problem).
            - 'satisfied' is your claim about whether the formula is
              satisfied by this assignment.
            - 'diagnostics' can contain any extra scalar / small-array
              data you want to expose on the leaderboard.
        """
        raise NotImplementedError("decode() must be implemented by the competitor")


