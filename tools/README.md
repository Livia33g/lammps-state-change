# Moderator Tools

This directory contains automation scripts for processing submissions and maintaining the competition.

---

## 📋 Tool Overview

| Tool | Purpose | Usage |
|------|---------|-------|
| `validate_submission.py` | Validate submission format and security | Run before accepting submissions |
| `process_new_submissions.sh` | Batch process submissions from inbox | Run periodically to handle new submissions |
| `evaluate_submission.slurm` | Evaluate single submission on cluster | Called by process script or manually |
| `update_leaderboard.py` | Update leaderboard from scores | Run after evaluation completes |

---

## 🔧 Individual Tools

### 1. validate_submission.py

**Purpose:** Validate submission structure, JSON syntax, and security before processing.

**Usage:**
```bash
python3 tools/validate_submission.py submissions-inbox/alice_problem-001/
```

**Checks:**
- Required files present (`policy.json`, `submission.json`)
- Valid JSON syntax
- Schema compliance
- Parameter ranges
- Suspicious content detection
- File size limits

**Exit codes:**
- `0`: Validation passed
- `1`: Validation failed

**Example output:**
```
✓ Found policy.json
✓ Found submission.json
✓ Valid JSON syntax
✓ Policy schema valid
✓ Parameters within allowed ranges
✓ No security issues detected

Submission is valid!
```

---

### 2. process_new_submissions.sh

**Purpose:** Automatically process all submissions in the inbox directory.

**Usage:**
```bash
# Process all submissions in submissions-inbox/
bash tools/process_new_submissions.sh
```

**Directory structure:**
```
submissions-inbox/
├── alice_problem-001/
│   ├── policy.json
│   ├── params.json
│   └── submission.json
└── bob_problem-001/
    └── ...
```

**Workflow:**
1. Scans `submissions-inbox/` for new submissions
2. Validates each submission
3. Moves valid submissions to `submissions-private/`
4. Submits to cluster via SLURM
5. Archives processed submissions to `.processed/`

**Output:**
```
