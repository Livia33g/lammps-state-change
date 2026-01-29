"""
Local checker for advanced state-change submissions.

Usage (from the repo root or the advance/ directory):

    python advance/check_submission.py \\
        --submission path/to/my_submission.py \\
        --problem problems/problem-001-ksat-advanced/problem.json

This script performs **structural** checks only:
  - Imports the submission module.
  - Validates API_VERSION and the StateChangeSolution class.
  - Runs encode() and design_policy() in a temporary work directory.
  - Confirms that required folders / files exist.
  - Optionally runs decode() on an empty or mock simulation directory.

It does NOT compile LAMMPS or run a full-length simulation; those steps
are handled by the competition infrastructure on the cluster.
"""

import argparse
import importlib.util
import json
import re
import shutil
import tempfile
from pathlib import Path


REQUIRED_API_VERSION = "v1"
CLASS_NAME = "StateChangeSolution"

def _read_text(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8", errors="replace")
    except Exception:
        return path.read_text(errors="replace")

def _check_lammps_input_has_masses(sim_dir: Path) -> None:
    """
    Heuristic check: require at least one `mass <type> <value>` in at least one in.* file.
    This prevents a very common sandbox failure:
      "Not all per-type masses are set ..."
    """
    in_files = list(sim_dir.glob("in.*"))
    for p in in_files:
        txt = _read_text(p)
        if re.search(r'^\s*mass\s+\d+\s+[-+0-9.eE]+\s*$', txt, flags=re.MULTILINE):
            return
    raise ValueError(
        "LAMMPS input is missing `mass <type> <value>` commands. "
        "For atom_style atomic (and most styles), you MUST set per-type masses "
        "before `velocity`/integration."
    )

def _check_fix_registration_pattern(gen_dir: Path) -> None:
    """
    Enforce the registration pattern that matches the moderator LAMMPS tree:
      - FixStyle(...) must appear in the HEADER under `#ifdef FIX_CLASS`
      - Do NOT include fix_style.h (not present in moderator build)
    """
    headers = sorted(gen_dir.glob("fix_state_change_*.h"))
    sources = sorted(gen_dir.glob("fix_state_change_*.cpp"))

    if not headers or not sources:
        # structural checker already handles missing files, but keep message crisp
        return

    for h in headers:
        txt = _read_text(h)
        if "fix_style.h" in txt:
            raise ValueError(
                f"{h.name} includes fix_style.h, which is not available on the moderator build. "
                "Use the `#ifdef FIX_CLASS` / `FixStyle(...)` header pattern instead."
            )
        if "#ifdef FIX_CLASS" not in txt or "FixStyle(" not in txt:
            raise ValueError(
                f"{h.name} does not contain the required registration block. "
                "Your header must include `#ifdef FIX_CLASS` and a `FixStyle(state/change/<name>,FixClass);` line."
            )
        # Ensure FixStyle is in the FIX_CLASS section (rough heuristic)
        parts = txt.split("#else", 1)
        if parts:
            fix_class_block = parts[0]
            if "FixStyle(" not in fix_class_block:
                raise ValueError(
                    f"{h.name} has FixStyle but not in the `#ifdef FIX_CLASS` block. "
                    "Move FixStyle(...) above the #else."
                )

    for cpp in sources:
        txt = _read_text(cpp)
        if "fix_style.h" in txt:
            raise ValueError(
                f"{cpp.name} includes fix_style.h, which is not available on the moderator build. "
                "Do not include it; register via the header's `#ifdef FIX_CLASS` block."
            )

def load_submission_module(path: Path):
    if not path.exists():
        raise FileNotFoundError(f"Submission file not found: {path}")

    spec = importlib.util.spec_from_file_location("submission_module", str(path))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load submission module from {path}")

    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[arg-type]
    return module


def ensure_api(module):
    api_version = getattr(module, "API_VERSION", None)
    if api_version != REQUIRED_API_VERSION:
        raise ValueError(
            f"API_VERSION mismatch: expected '{REQUIRED_API_VERSION}', "
            f"got '{api_version}'"
        )

    if not hasattr(module, CLASS_NAME):
        raise ValueError(f"Missing required class '{CLASS_NAME}' in submission.")

    cls = getattr(module, CLASS_NAME)
    if not callable(cls):
        raise TypeError(f"{CLASS_NAME} is not callable (not a class).")

    return cls


def run_encode_and_design(SolutionClass, problem_def: dict):
    work_dir = Path(tempfile.mkdtemp(prefix="state_change_advance_")).resolve()
    print(f"[checker] Using temporary work_dir: {work_dir}")

    try:
        sol = SolutionClass(problem_def, str(work_dir))

        # ----------------------------- encode -----------------------------
        print("[checker] Running encode() ...")
        system_meta = sol.encode()
        if not isinstance(system_meta, dict):
            raise TypeError("encode() must return a dict.")

        sim_dir = work_dir / "simulation"
        if not sim_dir.is_dir():
            raise FileNotFoundError(
                f"encode() must create directory: {sim_dir}"
            )

        in_files = list(sim_dir.glob("in.*"))
        data_files = list(sim_dir.glob("data.*"))
        if not in_files:
            raise FileNotFoundError(
                f"No LAMMPS input files found in {sim_dir} (expected 'in.*')."
            )
        if not data_files:
            raise FileNotFoundError(
                f"No LAMMPS data files found in {sim_dir} (expected 'data.*')."
            )
        _check_lammps_input_has_masses(sim_dir)

        # ------------------------- design_policy -------------------------
        print("[checker] Running design_policy() ...")
        policy_info = sol.design_policy(system_meta)
        if not isinstance(policy_info, dict):
            raise TypeError("design_policy() must return a dict.")

        fix_files = policy_info.get("fix_files")
        if not isinstance(fix_files, list) or not fix_files:
            raise ValueError(
                "design_policy() must return a dict with a non-empty "
                "'fix_files' list."
            )

        gen_dir = work_dir / "generated"
        if not gen_dir.is_dir():
            raise FileNotFoundError(
                f"design_policy() must create directory: {gen_dir}"
            )

        for fname in fix_files:
            if not isinstance(fname, str):
                raise TypeError("Each entry in 'fix_files' must be a string.")
            # Basic naming rule: must start with fix_state_change_
            if not fname.startswith("fix_state_change_"):
                raise ValueError(
                    f"Invalid fix filename '{fname}'. "
                    "Files must start with 'fix_state_change_'."
                )
            fpath = gen_dir / fname
            if not fpath.exists():
                raise FileNotFoundError(
                    f"Listed fix file does not exist: {fpath}"
                )

        # Additional static checks that match common moderator sandbox failures
        _check_fix_registration_pattern(gen_dir)

        print("[checker] encode() and design_policy() passed structural checks.")

        # ----------------------------- decode -----------------------------
        # Optional: we call decode() once to verify the return type.
        print("[checker] Running decode() (basic return-type check) ...")
        decoded = sol.decode()
        if not isinstance(decoded, dict):
            raise TypeError("decode() must return a dict.")

        # Check required keys if present
        if "assignment" in decoded:
            assignment = decoded["assignment"]
            if not isinstance(assignment, list) or not all(
                isinstance(x, (int, bool)) for x in assignment
            ):
                raise TypeError(
                    "decoded['assignment'] must be a list of ints/bools."
                )

        if "satisfied" in decoded and not isinstance(decoded["satisfied"], bool):
            raise TypeError("decoded['satisfied'] must be a bool.")

        # Check JSON serializability
        try:
            json.dumps(decoded)
        except TypeError as e:
            raise TypeError(
                f"decode() must return a JSON-serializable dict. Error: {e}"
            )

        print("[checker] decode() return value looks valid and serializable.")

    finally:
        # Always clean up temp directory to avoid clutter
        try:
            shutil.rmtree(work_dir)
        except Exception:
            pass


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Structural checker for advanced state-change submissions."
    )
    parser.add_argument(
        "--submission",
        required=True,
        help="Path to the competitor's submission.py file.",
    )
    parser.add_argument(
        "--problem",
        required=True,
        help="Path to an example problem.json to pass into the solution.",
    )

    args = parser.parse_args(argv)

    submission_path = Path(args.submission).resolve()
    problem_path = Path(args.problem).resolve()

    if not problem_path.exists():
        raise FileNotFoundError(f"Problem file not found: {problem_path}")

    problem_def = json.loads(problem_path.read_text())
    if not isinstance(problem_def, dict):
        raise TypeError("Problem JSON must decode to a dict.")

    module = load_submission_module(submission_path)
    SolutionClass = ensure_api(module)

    run_encode_and_design(SolutionClass, problem_def)

    print("\n[checker] All structural checks passed ✅")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


