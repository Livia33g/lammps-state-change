#!/usr/bin/env python3
"""
Validate a submission before evaluation.

Checks:
- Required files present (policy.json, submission.json)
- JSON syntax valid
- Policy conforms to schema
- Parameters in allowed ranges
- No suspicious content

Usage:
    python validate_submission.py submissions/problem-001/username/
"""

import json
import sys
from pathlib import Path
import argparse


def validate_submission(submission_dir: Path) -> bool:
    """Validate a submission directory. Returns True if valid."""

    errors = []
    warnings = []

    print(f"Validating submission: {submission_dir}")

    # Check required files
    required_files = ['policy.json', 'submission.json']
    for fname in required_files:
        fpath = submission_dir / fname
        if not fpath.exists():
            errors.append(f"Missing required file: {fname}")

    # Optional files
    params_file = submission_dir / 'params.json'
    has_params = params_file.exists()

    if errors:
        print("\n❌ ERRORS:")
        for err in errors:
            print(f"  - {err}")
        return False

    # Load and validate JSON syntax
    try:
        with open(submission_dir / 'policy.json') as f:
            policy = json.load(f)
        print("✓ policy.json: Valid JSON")
    except json.JSONDecodeError as e:
        errors.append(f"policy.json: Invalid JSON - {e}")
        return False

    try:
        with open(submission_dir / 'submission.json') as f:
            submission = json.load(f)
        print("✓ submission.json: Valid JSON")
    except json.JSONDecodeError as e:
        errors.append(f"submission.json: Invalid JSON - {e}")
        return False

    if has_params:
        try:
            with open(params_file) as f:
                params = json.load(f)
            print("✓ params.json: Valid JSON")
        except json.JSONDecodeError as e:
            errors.append(f"params.json: Invalid JSON - {e}")
            return False

    # Validate policy structure
    if 'policy_name' not in policy:
        errors.append("policy.json: Missing 'policy_name' field")
    if 'state_transitions' not in policy:
        errors.append("policy.json: Missing 'state_transitions' field")
    elif not isinstance(policy['state_transitions'], list):
        errors.append("policy.json: 'state_transitions' must be an array")
    elif len(policy['state_transitions']) == 0:
        errors.append("policy.json: 'state_transitions' cannot be empty")

    # Validate submission metadata
    if 'problem_id' not in submission:
        errors.append("submission.json: Missing 'problem_id' field")
    if 'team_name' not in submission:
        errors.append("submission.json: Missing 'team_name' field")

    # Check parameter ranges (basic validation)
    if has_params:
        for key, value in params.items():
            if not isinstance(value, (int, float)):
                warnings.append(f"params.json: {key} is not a number")
            if isinstance(value, (int, float)) and value < 0:
                warnings.append(f"params.json: {key} is negative (may be invalid)")

    # Check for suspicious content
    policy_str = json.dumps(policy)
    suspicious_patterns = ['rm -rf', 'system(', 'exec(', '__import__', 'eval(']
    for pattern in suspicious_patterns:
        if pattern in policy_str:
            errors.append(f"Suspicious content detected: {pattern}")

    # Print results
    if warnings:
        print("\n⚠ WARNINGS:")
        for warn in warnings:
            print(f"  - {warn}")

    if errors:
        print("\n❌ VALIDATION FAILED:")
        for err in errors:
            print(f"  - {err}")
        return False

    print("\n✅ VALIDATION PASSED")
    return True


def main():
    parser = argparse.ArgumentParser(description='Validate competition submission')
    parser.add_argument('submission_dir', help='Path to submission directory')
    args = parser.parse_args()

    submission_dir = Path(args.submission_dir)
    if not submission_dir.exists():
        print(f"Error: Directory not found: {submission_dir}")
        sys.exit(1)

    valid = validate_submission(submission_dir)
    sys.exit(0 if valid else 1)


if __name__ == '__main__':
    main()
