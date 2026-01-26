#!/usr/bin/env python3
"""
Workflow management system for competition submissions.

Inspired by signac-flow - action-based workflow with status tracking.

Usage:
    workflow.py status <problem_id>           # Show current state
    workflow.py validate <problem_id>         # Validate pending submissions
    workflow.py generate <problem_id>         # Generate C++ fixes
    workflow.py compile <problem_id>          # Compile LAMMPS
    workflow.py simulate <problem_id>         # Run simulations
    workflow.py analyze <problem_id>          # Analyze results
    workflow.py cleanup <problem_id>          # Remove intermediate files
    workflow.py update-leaderboard <problem_id>  # Update public leaderboard
    workflow.py run-all <problem_id>          # Execute all pending actions
    workflow.py process-new <problem_id>      # Process only new submissions

Actions are idempotent - won't repeat completed work.
"""

import argparse
import json
import os
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Set
from dataclasses import dataclass, asdict
from enum import Enum

try:
    from tqdm import tqdm
except ImportError:
    print("Warning: tqdm not installed. Install with: pip install tqdm")
    # Fallback implementation
    class tqdm:
        def __init__(self, iterable=None, desc=None, total=None, **kwargs):
            self.iterable = iterable
            self.desc = desc
            self.total = total or (len(iterable) if iterable else 0)
            self.n = 0

        def __iter__(self):
            for item in self.iterable:
                yield item
                self.update(1)

        def update(self, n=1):
            self.n += n
            if self.desc:
                print(f"{self.desc}: {self.n}/{self.total}", end='\r')

        def __enter__(self):
            return self

        def __exit__(self, *args):
            if self.desc:
                print()


class ActionStatus(Enum):
    """Status of an action for a submission."""
    PENDING = "pending"
    RUNNING = "running"
    COMPLETE = "complete"
    FAILED = "failed"
    SKIPPED = "skipped"


@dataclass
class SubmissionStatus:
    """Track status of all actions for a submission."""
    submission_id: str
    username: str
    timestamp: str
    problem_id: str

    # Action statuses
    validated: str = ActionStatus.PENDING.value
    generated: str = ActionStatus.PENDING.value
    compiled: str = ActionStatus.PENDING.value
    simulated: str = ActionStatus.PENDING.value
    analyzed: str = ActionStatus.PENDING.value
    cleaned: str = ActionStatus.PENDING.value

    # Results
    score: Optional[float] = None
    error_message: Optional[str] = None

    # Timestamps
    created_at: Optional[str] = None
    completed_at: Optional[str] = None


class WorkflowManager:
    """Manage submission workflow for a problem."""

    def __init__(self, problem_id: str, private_repo: Optional[Path] = None):
        self.problem_id = problem_id

        if private_repo is None:
            private_repo = Path(os.environ.get(
                'SUBMISSIONS_PRIVATE_REPO',
                Path.home() / 'lammps-state-change-private'
            ))

        self.private_repo = Path(private_repo)
        self.problem_dir = self.private_repo / 'submissions-archive' / problem_id
        self.workflow_dir = self.problem_dir / '.workflow'
        self.status_file = self.workflow_dir / 'status.json'
        self.scores_file = self.workflow_dir / 'all_scores.json'

        # Create directories
        self.workflow_dir.mkdir(parents=True, exist_ok=True)

        # Load or initialize status
        self.statuses: Dict[str, SubmissionStatus] = {}
        self.load_status()

    def load_status(self):
        """Load workflow status from disk."""
        if self.status_file.exists():
            with open(self.status_file) as f:
                data = json.load(f)
                self.statuses = {
                    sid: SubmissionStatus(**status_data)
                    for sid, status_data in data.items()
                }

    def save_status(self):
        """Save workflow status to disk."""
        data = {
            sid: asdict(status)
            for sid, status in self.statuses.items()
        }
        with open(self.status_file, 'w') as f:
            json.dump(data, f, indent=2)

    def scan_submissions(self) -> List[Path]:
        """Find all submission directories."""
        if not self.problem_dir.exists():
            return []

        submissions = []
        for path in self.problem_dir.iterdir():
            if path.is_dir() and not path.name.startswith('.'):
                # Check for required files
                if (path / 'policy.json').exists() and (path / 'submission.json').exists():
                    submissions.append(path)

        return sorted(submissions)

    def get_submission_id(self, submission_path: Path) -> str:
        """Get unique submission ID from path."""
        return submission_path.name

    def get_username(self, submission_path: Path) -> str:
        """Extract username from submission directory name."""
        # Format: username_timestamp
        return submission_path.name.split('_')[0]

    def register_submission(self, submission_path: Path):
        """Register a new submission in the workflow."""
        submission_id = self.get_submission_id(submission_path)

        if submission_id not in self.statuses:
            username = self.get_username(submission_path)

            # Extract timestamp from directory name
            parts = submission_path.name.split('_')
            timestamp = '_'.join(parts[1:]) if len(parts) > 1 else datetime.now().isoformat()

            self.statuses[submission_id] = SubmissionStatus(
                submission_id=submission_id,
                username=username,
                timestamp=timestamp,
                problem_id=self.problem_id,
                created_at=datetime.now().isoformat()
            )
            self.save_status()

    def update_action_status(self, submission_id: str, action: str, status: ActionStatus, error: Optional[str] = None):
        """Update status of an action."""
        if submission_id in self.statuses:
            setattr(self.statuses[submission_id], action, status.value)
            if error:
                self.statuses[submission_id].error_message = error
            if status == ActionStatus.COMPLETE and action == "analyzed":
                self.statuses[submission_id].completed_at = datetime.now().isoformat()
            self.save_status()

    def get_pending_actions(self, action: str) -> List[tuple]:
        """Get submissions that need this action."""
        pending = []

        for submission_id, status in self.statuses.items():
            # Check if this action is pending and dependencies are met
            if getattr(status, action) == ActionStatus.PENDING.value:
                if self._check_dependencies(status, action):
                    submission_path = self.problem_dir / submission_id
                    if submission_path.exists():
                        pending.append((submission_id, submission_path))

        return pending

    def _check_dependencies(self, status: SubmissionStatus, action: str) -> bool:
        """Check if dependencies for an action are satisfied."""
        dependencies = {
            'validated': [],
            'generated': ['validated'],
            'compiled': ['generated'],
            'simulated': ['compiled'],
            'analyzed': ['simulated'],
            'cleaned': ['analyzed'],
        }

        if action not in dependencies:
            return True

        for dep in dependencies[action]:
            if getattr(status, dep) != ActionStatus.COMPLETE.value:
                return False

        return True

    def print_status(self):
        """Print current workflow status."""
        print(f"\n{'='*80}")
        print(f"  Workflow Status: {self.problem_id}")
        print(f"{'='*80}\n")

        if not self.statuses:
            print("No submissions found.\n")
            return

        # Count by status
        total = len(self.statuses)
        by_action = {
            'validated': {'complete': 0, 'pending': 0, 'failed': 0},
            'generated': {'complete': 0, 'pending': 0, 'failed': 0},
            'compiled': {'complete': 0, 'pending': 0, 'failed': 0},
            'simulated': {'complete': 0, 'pending': 0, 'failed': 0},
            'analyzed': {'complete': 0, 'pending': 0, 'failed': 0},
        }

        for status in self.statuses.values():
            for action in by_action:
                action_status = getattr(status, action)
                if action_status == ActionStatus.COMPLETE.value:
                    by_action[action]['complete'] += 1
                elif action_status == ActionStatus.FAILED.value:
                    by_action[action]['failed'] += 1
                else:
                    by_action[action]['pending'] += 1

        print(f"Total submissions: {total}\n")
        print(f"{'Action':<15} {'Complete':<10} {'Pending':<10} {'Failed':<10}")
        print(f"{'-'*50}")
        for action, counts in by_action.items():
            print(f"{action:<15} {counts['complete']:<10} {counts['pending']:<10} {counts['failed']:<10}")

        print(f"\n{'Next available actions:'}")
        print(f"{'-'*50}")

        # Show what can be done next
        actions_available = []
        for action in ['validated', 'generated', 'compiled', 'simulated', 'analyzed', 'cleaned']:
            pending = self.get_pending_actions(action)
            if pending:
                actions_available.append(f"  • {action}: {len(pending)} submissions ready")

        if actions_available:
            for line in actions_available:
                print(line)
        else:
            print("  All actions complete! ✓")

        print()

    def validate_submissions(self, concurrent: bool = False):
        """Validate all pending submissions."""
        pending = self.get_pending_actions('validated')

        if not pending:
            print("No submissions to validate.")
            return

        print(f"\nValidating {len(pending)} submissions...")

        for submission_id, submission_path in tqdm(pending, desc="Validating"):
            try:
                # Run validation
                result = subprocess.run(
                    ['python3', 'tools/validate_submission.py', str(submission_path)],
                    capture_output=True,
                    text=True,
                    timeout=30
                )

                if result.returncode == 0:
                    self.update_action_status(submission_id, 'validated', ActionStatus.COMPLETE)
                else:
                    self.update_action_status(
                        submission_id,
                        'validated',
                        ActionStatus.FAILED,
                        error=result.stderr
                    )
            except Exception as e:
                self.update_action_status(
                    submission_id,
                    'validated',
                    ActionStatus.FAILED,
                    error=str(e)
                )

        print(f"✓ Validation complete\n")

    def generate_cpp(self, concurrent: bool = False):
        """Generate C++ fixes for validated submissions."""
        pending = self.get_pending_actions('generated')

        if not pending:
            print("No submissions to generate C++.")
            return

        print(f"\nGenerating C++ for {len(pending)} submissions...")

        for submission_id, submission_path in tqdm(pending, desc="Generating"):
            try:
                policy_file = submission_path / 'policy.json'
                output_dir = submission_path / 'generated'
                output_dir.mkdir(exist_ok=True)

                result = subprocess.run(
                    ['python3', 'core/generators/generate_fix_from_policy.py',
                     str(policy_file), str(output_dir)],
                    capture_output=True,
                    text=True,
                    timeout=60
                )

                if result.returncode == 0:
                    self.update_action_status(submission_id, 'generated', ActionStatus.COMPLETE)
                else:
                    self.update_action_status(
                        submission_id,
                        'generated',
                        ActionStatus.FAILED,
                        error=result.stderr
                    )
            except Exception as e:
                self.update_action_status(
                    submission_id,
                    'generated',
                    ActionStatus.FAILED,
                    error=str(e)
                )

        print(f"✓ C++ generation complete\n")

    def compile_lammps(self, concurrent: bool = False):
        """Compile LAMMPS with generated fixes (sequential - one at a time)."""
        pending = self.get_pending_actions('compiled')

        if not pending:
            print("No submissions to compile.")
            return

        print(f"\nCompiling LAMMPS for {len(pending)} submissions...")
        print("Note: Compilation is sequential to avoid conflicts\n")

        for submission_id, submission_path in tqdm(pending, desc="Compiling"):
            try:
                generated_dir = submission_path / 'generated'

                # Find generated C++ files
                cpp_files = list(generated_dir.glob('fix_state_change_*.cpp'))
                h_files = list(generated_dir.glob('fix_state_change_*.h'))

                if not cpp_files or not h_files:
                    raise FileNotFoundError("Generated C++ files not found")

                # Get LAMMPS source directory
                lammps_src = Path(os.environ.get('LAMMPS_SRC', ''))
                if not lammps_src.exists():
                    raise FileNotFoundError("LAMMPS_SRC not set or doesn't exist")

                # Copy files to LAMMPS src
                for cpp_file in cpp_files:
                    subprocess.run(['cp', str(cpp_file), str(lammps_src)], check=True)
                for h_file in h_files:
                    subprocess.run(['cp', str(h_file), str(lammps_src)], check=True)

                # Compile
                result = subprocess.run(
                    ['make', 'mpi', '-j8'],
                    cwd=str(lammps_src),
                    capture_output=True,
                    text=True,
                    timeout=600  # 10 minute timeout
                )

                if result.returncode == 0:
                    self.update_action_status(submission_id, 'compiled', ActionStatus.COMPLETE)
                else:
                    self.update_action_status(
                        submission_id,
                        'compiled',
                        ActionStatus.FAILED,
                        error=result.stderr[:500]  # Truncate long errors
                    )
            except Exception as e:
                self.update_action_status(
                    submission_id,
                    'compiled',
                    ActionStatus.FAILED,
                    error=str(e)
                )

        print(f"✓ Compilation complete\n")

    def run_simulations(self, concurrent: int = 4):
        """Run LAMMPS simulations (concurrent with job limit)."""
        pending = self.get_pending_actions('simulated')

        if not pending:
            print("No submissions to simulate.")
            return

        print(f"\nRunning simulations for {len(pending)} submissions...")
        print(f"Concurrent jobs: {concurrent}\n")

        # Check if SLURM is available
        has_slurm = subprocess.run(['which', 'sbatch'], capture_output=True).returncode == 0

        if has_slurm:
            self._run_simulations_slurm(pending)
        else:
            self._run_simulations_local(pending, concurrent)

    def _run_simulations_slurm(self, pending: List[tuple]):
        """Submit simulations to SLURM."""
        job_ids = {}

        for submission_id, submission_path in tqdm(pending, desc="Submitting"):
            try:
                # Submit SLURM job
                result = subprocess.run(
                    ['sbatch', 'tools/evaluate_submission.slurm',
                     self.problem_id, submission_id],
                    capture_output=True,
                    text=True
                )

                if result.returncode == 0:
                    # Extract job ID
                    job_id = result.stdout.strip().split()[-1]
                    job_ids[submission_id] = job_id
                    self.update_action_status(submission_id, 'simulated', ActionStatus.RUNNING)
                else:
                    self.update_action_status(
                        submission_id,
                        'simulated',
                        ActionStatus.FAILED,
                        error=result.stderr
                    )
            except Exception as e:
                self.update_action_status(
                    submission_id,
                    'simulated',
                    ActionStatus.FAILED,
                    error=str(e)
                )

        print(f"\n✓ Submitted {len(job_ids)} jobs to SLURM")
        print(f"Monitor with: squeue -u $USER\n")

        # Save job IDs
        jobs_file = self.workflow_dir / 'slurm_jobs.json'
        with open(jobs_file, 'w') as f:
            json.dump(job_ids, f, indent=2)

    def _run_simulations_local(self, pending: List[tuple], concurrent: int):
        """Run simulations locally (sequential for now)."""
        for submission_id, submission_path in tqdm(pending, desc="Simulating"):
            try:
                generated_dir = submission_path / 'generated'

                # Find LAMMPS input file
                input_files = list(generated_dir.glob('in.*.lammps'))
                if not input_files:
                    raise FileNotFoundError("LAMMPS input file not found")

                input_file = input_files[0]
                lammps_bin = Path(os.environ.get('LAMMPS_SRC', '')) / '..' / 'src' / 'lmp_mpi'

                # Run LAMMPS
                with open(generated_dir / 'lammps_stdout.log', 'w') as stdout:
                    with open(generated_dir / 'lammps_stderr.log', 'w') as stderr:
                        result = subprocess.run(
                            ['mpirun', '-np', '4', str(lammps_bin), '-in', str(input_file)],
                            cwd=str(generated_dir),
                            stdout=stdout,
                            stderr=stderr,
                            timeout=86400  # 24 hour timeout
                        )

                if result.returncode == 0:
                    self.update_action_status(submission_id, 'simulated', ActionStatus.COMPLETE)
                else:
                    self.update_action_status(
                        submission_id,
                        'simulated',
                        ActionStatus.FAILED,
                        error="LAMMPS exited with error"
                    )
            except Exception as e:
                self.update_action_status(
                    submission_id,
                    'simulated',
                    ActionStatus.FAILED,
                    error=str(e)
                )

        print(f"✓ Simulations complete\n")

    def analyze_results(self, concurrent: bool = False):
        """Analyze simulation results."""
        pending = self.get_pending_actions('analyzed')

        if not pending:
            print("No simulations to analyze.")
            return

        print(f"\nAnalyzing {len(pending)} simulations...")

        for submission_id, submission_path in tqdm(pending, desc="Analyzing"):
            try:
                generated_dir = submission_path / 'generated'
                analysis_dir = generated_dir / 'analysis'
                analysis_dir.mkdir(exist_ok=True)

                # Find trajectory and log files
                dump_files = list(generated_dir.glob('dump.*.lammpstrj'))
                log_file = generated_dir / 'lammps_stdout.log'
                stderr_file = generated_dir / 'lammps_stderr.log'

                if not dump_files or not log_file.exists():
                    raise FileNotFoundError("Trajectory or log files not found")

                # Run analysis
                result = subprocess.run(
                    ['python3', 'core/benchmark/score_policy_from_timeseries.py',
                     '--dump', str(dump_files[0]),
                     '--thermo', str(log_file),
                     '--events', str(stderr_file) if stderr_file.exists() else '',
                     '--output', str(analysis_dir)],
                    capture_output=True,
                    text=True,
                    timeout=300
                )

                if result.returncode == 0:
                    # Load score
                    scores_file = analysis_dir / 'scores.json'
                    if scores_file.exists():
                        with open(scores_file) as f:
                            scores = json.load(f)
                            score = scores.get('work_per_yield', float('inf'))
                            self.statuses[submission_id].score = score

                    self.update_action_status(submission_id, 'analyzed', ActionStatus.COMPLETE)
                else:
                    self.update_action_status(
                        submission_id,
                        'analyzed',
                        ActionStatus.FAILED,
                        error=result.stderr[:500]
                    )
            except Exception as e:
                self.update_action_status(
                    submission_id,
                    'analyzed',
                    ActionStatus.FAILED,
                    error=str(e)
                )

        print(f"✓ Analysis complete\n")

    def cleanup_intermediate_files(self):
        """Remove intermediate files after successful analysis."""
        eligible = []

        for submission_id, status in self.statuses.items():
            if status.analyzed == ActionStatus.COMPLETE.value and status.cleaned == ActionStatus.PENDING.value:
                submission_path = self.problem_dir / submission_id
                if submission_path.exists():
                    eligible.append((submission_id, submission_path))

        if not eligible:
            print("No submissions to clean.")
            return

        print(f"\nCleaning up {len(eligible)} submissions...")

        for submission_id, submission_path in tqdm(eligible, desc="Cleaning"):
            try:
                generated_dir = submission_path / 'generated'

                # Remove large intermediate files
                files_to_remove = [
                    'dump.*.lammpstrj',      # Trajectory files (large)
                    'fix_state_change_*.cpp', # C++ source (can regenerate)
                    'fix_state_change_*.h',   # C++ header (can regenerate)
                    'data.*.lammps',         # Data file (can regenerate)
                    'in.*.lammps',           # Input script (can regenerate)
                    'lammps_stdout.log',     # LAMMPS stdout (keep stderr for events)
                ]

                for pattern in files_to_remove:
                    for file in generated_dir.glob(pattern):
                        file.unlink()

                self.update_action_status(submission_id, 'cleaned', ActionStatus.COMPLETE)

            except Exception as e:
                self.update_action_status(
                    submission_id,
                    'cleaned',
                    ActionStatus.FAILED,
                    error=str(e)
                )

        print(f"✓ Cleanup complete\n")

    def update_leaderboard(self):
        """Update public leaderboard with best scores per user."""
        print("\nUpdating leaderboard...")

        # Collect all scores by user
        user_scores: Dict[str, List[tuple]] = {}  # username -> [(score, submission_id, timestamp)]

        for submission_id, status in self.statuses.items():
            if status.analyzed == ActionStatus.COMPLETE.value and status.score is not None:
                username = status.username
                if username not in user_scores:
                    user_scores[username] = []
                user_scores[username].append((status.score, submission_id, status.timestamp))

        # Keep only best score per user (lower is better for work_per_yield)
        best_submissions = {}
        for username, scores in user_scores.items():
            scores.sort()  # Sort by score (ascending)
            best_score, best_id, timestamp = scores[0]
            best_submissions[username] = (best_score, best_id, timestamp)

        print(f"Found {len(best_submissions)} users with submissions")

        # Update leaderboard using existing tool
        for username, (score, submission_id, timestamp) in tqdm(best_submissions.items(), desc="Updating"):
            submission_path = self.problem_dir / submission_id
            try:
                subprocess.run(
                    ['python3', 'tools/update_leaderboard.py',
                     '--submission', str(submission_path)],
                    check=True,
                    capture_output=True
                )
            except Exception as e:
                print(f"Warning: Failed to update leaderboard for {username}: {e}")

        print(f"✓ Leaderboard updated with best scores\n")

    def run_all_pending(self, concurrent: int = 4):
        """Execute all pending actions in order."""
        print("\nRunning all pending actions...")

        actions = [
            ('validate', lambda: self.validate_submissions()),
            ('generate', lambda: self.generate_cpp()),
            ('compile', lambda: self.compile_lammps()),
            ('simulate', lambda: self.run_simulations(concurrent)),
            ('analyze', lambda: self.analyze_results()),
            ('cleanup', lambda: self.cleanup_intermediate_files()),
            ('update-leaderboard', lambda: self.update_leaderboard()),
        ]

        for action_name, action_func in actions:
            pending = self.get_pending_actions(action_name)
            if pending:
                print(f"\n{'-'*60}")
                print(f"  Action: {action_name} ({len(pending)} submissions)")
                print(f"{'-'*60}")
                action_func()
            else:
                print(f"✓ {action_name}: nothing to do")

        print("\n✓ All actions complete!")
        self.print_status()

    def process_new_only(self, concurrent: int = 4):
        """Process only submissions that haven't started processing yet."""
        new_submissions = []

        for submission_id, status in self.statuses.items():
            # Check if this is a new submission (validated is still pending)
            if status.validated == ActionStatus.PENDING.value:
                submission_path = self.problem_dir / submission_id
                if submission_path.exists():
                    new_submissions.append(submission_id)

        if not new_submissions:
            print("No new submissions to process.")
            return

        print(f"\nFound {len(new_submissions)} new submissions")
        print("Processing new submissions only...\n")

        self.run_all_pending(concurrent)


def main():
    parser = argparse.ArgumentParser(
        description='Workflow management for competition submissions',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('action', choices=[
        'status', 'validate', 'generate', 'compile', 'simulate',
        'analyze', 'cleanup', 'update-leaderboard', 'run-all', 'process-new'
    ], help='Action to perform')

    parser.add_argument('problem_id', help='Problem ID (e.g., problem-001-dimer-ksat)')
    parser.add_argument('--concurrent', '-j', type=int, help='Number of parallel jobs')

    args = parser.parse_args()

    # Create workflow manager
    wf = WorkflowManager(args.problem_id)

    # Scan for submissions and register new ones
    submissions = wf.scan_submissions()
    for sub_path in submissions:
        wf.register_submission(sub_path)

    # Execute action
    if args.action == 'status':
        wf.print_status()

    elif args.action == 'validate':
        wf.validate_submissions(concurrent=args.concurrent)
        wf.print_status()

    elif args.action == 'generate':
        wf.generate_cpp(concurrent=args.concurrent)
        wf.print_status()

    elif args.action == 'compile':
        wf.compile_lammps()
        wf.print_status()

    elif args.action == 'simulate':
        concurrent = args.concurrent or 4
        wf.run_simulations(concurrent)
        wf.print_status()

    elif args.action == 'analyze':
        wf.analyze_results()
        wf.print_status()

    elif args.action == 'cleanup':
        wf.cleanup_intermediate_files()
        wf.print_status()

    elif args.action == 'update-leaderboard':
        wf.update_leaderboard()

    elif args.action == 'run-all':
        concurrent = args.concurrent or 4
        wf.run_all_pending(concurrent)

    elif args.action == 'process-new':
        concurrent = args.concurrent or 4
        wf.process_new_only(concurrent)


if __name__ == '__main__':
    main()
