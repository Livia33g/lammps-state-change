#!/usr/bin/env python3
"""
Update leaderboard CSV from evaluation results.

This script reads scores.json from an evaluated submission and updates
the leaderboard CSV file for the corresponding problem.

Usage:
    python3 tools/update_leaderboard.py \
        --problem problem-001-dimer-ksat \
        --username alice \
        --scores submissions-private/problem-001-dimer-ksat/alice_2026-01-23/generated/analysis/scores.json

    # Or auto-detect paths:
    python3 tools/update_leaderboard.py \
        --submission submissions-private/problem-001-dimer-ksat/alice_2026-01-23/

Features:
    - Automatically updates or appends leaderboard entries
    - Sorts by primary metric (work_per_yield by default)
    - Preserves existing entries
    - Creates leaderboard if it doesn't exist
    - Validates scores before updating
"""

import argparse
import csv
import json
import sys
from pathlib import Path
from datetime import datetime
from typing import Dict, Any, List, Optional


def load_scores(scores_path: Path) -> Dict[str, Any]:
    """Load scores from JSON file."""
    try:
        with open(scores_path) as f:
            scores = json.load(f)
        return scores
    except FileNotFoundError:
        print(f"Error: Scores file not found: {scores_path}", file=sys.stderr)
        sys.exit(1)
    except json.JSONDecodeError as e:
        print(f"Error: Invalid JSON in scores file: {e}", file=sys.stderr)
        sys.exit(1)


def load_problem_config(problem_id: str) -> Dict[str, Any]:
    """Load problem.json to get scoring configuration."""
    problem_path = Path("problems") / problem_id / "problem.json"

    if not problem_path.exists():
        # Use defaults if problem.json not found
        return {
            "scoring": {
                "primary_metric": "work_per_yield",
                "sort_order": "ascending"
            }
        }

    with open(problem_path) as f:
        return json.load(f)


def validate_scores(scores: Dict[str, Any]) -> bool:
    """Validate that scores contain required fields."""
    required_fields = ['final_yield', 'flip_count', 'cumulative_work']

    for field in required_fields:
        if field not in scores:
            print(f"Error: Missing required field '{field}' in scores", file=sys.stderr)
            return False

        if not isinstance(scores[field], (int, float)):
            print(f"Error: Field '{field}' must be numeric", file=sys.stderr)
            return False

    # Check for reasonable values
    if scores['final_yield'] < 0 or scores['final_yield'] > 1:
        print(f"Warning: final_yield={scores['final_yield']} outside [0,1]", file=sys.stderr)

    return True


def calculate_derived_metrics(scores: Dict[str, Any]) -> Dict[str, Any]:
    """Calculate work_per_yield if not already present."""
    if 'work_per_yield' not in scores and scores['final_yield'] > 0:
        scores['work_per_yield'] = scores['cumulative_work'] / scores['final_yield']
    elif scores['final_yield'] == 0:
        scores['work_per_yield'] = float('inf')

    return scores


def read_leaderboard(leaderboard_path: Path) -> List[Dict[str, Any]]:
    """Read existing leaderboard CSV."""
    if not leaderboard_path.exists():
        return []

    with open(leaderboard_path, 'r') as f:
        reader = csv.DictReader(f)
        return list(reader)


def write_leaderboard(leaderboard_path: Path, entries: List[Dict[str, Any]],
                     fieldnames: List[str]) -> None:
    """Write leaderboard CSV."""
    leaderboard_path.parent.mkdir(parents=True, exist_ok=True)

    with open(leaderboard_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(entries)


def format_value(value: Any) -> str:
    """Format value for CSV (round floats to reasonable precision)."""
    if isinstance(value, float):
        if value == float('inf'):
            return 'inf'
        return f"{value:.4f}"
    return str(value)


def update_leaderboard(problem_id: str, username: str, scores: Dict[str, Any],
                      date: Optional[str] = None, anonymize: bool = False) -> None:
    """Update leaderboard with new scores."""

    # Load problem configuration
    problem_config = load_problem_config(problem_id)
    primary_metric = problem_config['scoring'].get('primary_metric', 'work_per_yield')
    sort_order = problem_config['scoring'].get('sort_order', 'ascending')

    # Determine leaderboard path
    leaderboard_path = Path("problems") / problem_id / "leaderboard.csv"

    # Read existing entries
    entries = read_leaderboard(leaderboard_path)

    # Calculate derived metrics
    scores = calculate_derived_metrics(scores)

    # Create new entry
    if date is None:
        date = datetime.now().strftime("%Y-%m-%d")

    new_entry = {
        'username': username,
        'final_yield': format_value(scores.get('final_yield', 0)),
        'work_per_yield': format_value(scores.get('work_per_yield', float('inf'))),
        'flip_count': format_value(scores.get('flip_count', 0)),
        'cumulative_work': format_value(scores.get('cumulative_work', 0)),
        'date': date
    }

    # Add optional metrics if present
    if 'time_to_threshold' in scores:
        new_entry['time_to_threshold'] = format_value(scores['time_to_threshold'])

    # Check if username already exists - update or append
    existing_index = None
    for i, entry in enumerate(entries):
        if entry.get('username') == username:
            existing_index = i
            break

    if existing_index is not None:
        print(f"Updating existing entry for {username}")
        entries[existing_index] = new_entry
    else:
        print(f"Adding new entry for {username}")
        entries.append(new_entry)

    # Sort by primary metric
    reverse = (sort_order == 'descending')

    def sort_key(entry):
        value = entry.get(primary_metric, '')
        if value == 'inf':
            return float('inf')
        try:
            return float(value)
        except (ValueError, TypeError):
            return float('inf')

    entries.sort(key=sort_key, reverse=reverse)

    # Anonymize if requested
    if anonymize:
        for i, entry in enumerate(entries, 1):
            entry['rank'] = str(i)
            entry.pop('username', None)
        fieldnames = ['rank', 'final_yield', 'work_per_yield', 'flip_count',
                     'cumulative_work', 'date']
    else:
        fieldnames = ['username', 'final_yield', 'work_per_yield', 'flip_count',
                     'cumulative_work', 'date']
        if any('time_to_threshold' in e for e in entries):
            fieldnames.insert(-1, 'time_to_threshold')

    # Write updated leaderboard
    write_leaderboard(leaderboard_path, entries, fieldnames)

    print(f"✓ Updated leaderboard: {leaderboard_path}")
    print(f"  Total entries: {len(entries)}")
    print(f"  Primary metric: {primary_metric} ({sort_order})")

    # Show top 5
    print("\nTop 5:")
    for i, entry in enumerate(entries[:5], 1):
        user = entry.get('username', f"Rank {entry.get('rank', i)}")
        metric_value = entry.get(primary_metric, 'N/A')
        print(f"  {i}. {user}: {primary_metric}={metric_value}")


def main():
    parser = argparse.ArgumentParser(
        description='Update competition leaderboard from evaluation results',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Update from scores.json
  python3 tools/update_leaderboard.py \\
      --problem problem-001-dimer-ksat \\
      --username alice \\
      --scores submissions-private/problem-001/alice/generated/analysis/scores.json

  # Auto-detect from submission directory
  python3 tools/update_leaderboard.py \\
      --submission submissions-private/problem-001/alice_2026-01-23/

  # Anonymize leaderboard (remove usernames, show ranks)
  python3 tools/update_leaderboard.py \\
      --problem problem-001-dimer-ksat \\
      --username alice \\
      --scores path/to/scores.json \\
      --anonymize
        """
    )

    parser.add_argument('--problem', help='Problem ID (e.g., problem-001-dimer-ksat)')
    parser.add_argument('--username', help='Participant username')
    parser.add_argument('--scores', help='Path to scores.json')
    parser.add_argument('--submission', help='Path to submission directory (auto-detect paths)')
    parser.add_argument('--date', help='Submission date (YYYY-MM-DD, default: today)')
    parser.add_argument('--anonymize', action='store_true',
                       help='Anonymize leaderboard (remove usernames)')

    args = parser.parse_args()

    # Auto-detect mode
    if args.submission:
        submission_path = Path(args.submission)

        if not submission_path.exists():
            print(f"Error: Submission directory not found: {submission_path}", file=sys.stderr)
            sys.exit(1)

        # Extract problem and username from path
        # Expected: submissions-private/problem-id/username_timestamp/
        parts = submission_path.parts

        if 'submissions-private' in parts:
            idx = parts.index('submissions-private')
            if len(parts) > idx + 2:
                problem_id = parts[idx + 1]
                username_timestamp = parts[idx + 2]
                username = username_timestamp.split('_')[0]
            else:
                print("Error: Cannot parse submission path", file=sys.stderr)
                sys.exit(1)
        else:
            print("Error: Expected path format: submissions-private/problem-id/username/", file=sys.stderr)
            sys.exit(1)

        # Find scores.json
        scores_path = submission_path / "generated" / "analysis" / "scores.json"

        if not scores_path.exists():
            print(f"Error: scores.json not found at {scores_path}", file=sys.stderr)
            print("Make sure evaluation has completed successfully", file=sys.stderr)
            sys.exit(1)

    else:
        # Manual mode
        if not args.problem or not args.username or not args.scores:
            parser.print_help()
            print("\nError: --problem, --username, and --scores are required (or use --submission)",
                  file=sys.stderr)
            sys.exit(1)

        problem_id = args.problem
        username = args.username
        scores_path = Path(args.scores)

    # Load and validate scores
    scores = load_scores(scores_path)

    if not validate_scores(scores):
        print("Error: Invalid scores", file=sys.stderr)
        sys.exit(1)

    # Update leaderboard
    update_leaderboard(problem_id, username, scores, args.date, args.anonymize)


if __name__ == '__main__':
    main()
