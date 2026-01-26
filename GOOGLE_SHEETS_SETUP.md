# Google Sheets Submission Setup

**Use Google Sheets for submissions** - Easiest for both you and participants!

---

## 🎯 Why Google Sheets?

✅ **Participants:** Just fill out a form - no email, no file attachments
✅ **You:** Export CSV → run script → only new submissions extracted
✅ **Automatic tracking:** Script remembers which rows it already processed
✅ **No manual work:** Just export and run the script periodically

---

## 🚀 Quick Setup (5 Minutes)

### Step 1: Create Google Form

1. Go to https://forms.google.com
2. Create new form: **"Competition Submission - Problem 001"**

### Step 2: Add Form Fields

**Field 1: Username** (Short answer, Required)
```
Your username for the leaderboard
Example: alice
```

**Field 2: Problem ID** (Dropdown, Required)
```
Which problem are you submitting for?
Options:
  - problem-001-dimer-ksat
  - problem-002-octahedron
  - problem-003-hamiltonian
```

**Field 3: policy.json** (Paragraph, Required)
```
Paste the ENTIRE contents of your policy.json file here

IMPORTANT: Copy the entire file including all { } brackets.
Must be valid JSON!

Example:
{
  "policy_name": "my_policy",
  "check_frequency": 100,
  "state_transitions": [...]
}
```

**Field 4: params.json** (Paragraph, Optional)
```
If this problem allows custom parameters, paste your params.json here.
Otherwise leave blank.

Example:
{
  "morse_depth_AB": 1.2,
  "morse_alpha": 5.0
}
```

**Field 5: Email** (Short answer, Required, Email validation)
```
Your email for result notifications
```

### Step 3: Configure Form Settings

1. **Settings** → **Responses**:
   - ✅ Collect email addresses
   - ✅ Limit to 1 response (uncheck if you want multiple submissions per user)
   - ✅ Allow response editing (so users can update their submission)

2. **Settings** → **General**:
   - Anyone with the link can respond (or require sign-in)

### Step 4: Get Shareable Link

1. Click **Send** button
2. Click link icon 🔗
3. Copy the link
4. Share with participants!

### Step 5: Link to Google Sheets

1. Go to **Responses** tab
2. Click the Google Sheets icon
3. Create a new spreadsheet or select existing
4. Responses will auto-populate the sheet!

---

## 📥 Processing Submissions

### Daily Workflow (2 Commands)

```bash
# 1. Export Google Sheet as CSV
#    In Google Sheets: File → Download → Comma Separated Values (.csv)
#    Save to ~/Downloads/responses.csv

# 2. Extract ONLY new submissions (auto-tracks which rows already processed)
python3 tools/extract_from_sheets.py ~/Downloads/responses.csv problem-001-dimer-ksat

# Output:
#   ==========================================
#     Extract from Google Sheets
#   ==========================================
#     CSV file: ~/Downloads/responses.csv
#     Problem:  problem-001-dimer-ksat
#     Already processed: 5 rows
#
#   Total rows in CSV: 8
#
#   Processing 3 new rows...
#
#     ✓ Row 6: alice → ...submissions-inbox/alice_problem-001
#     ✓ Row 7: diana → ...submissions-inbox/diana_problem-001
#     ✓ Row 8: bob → ...submissions-inbox/bob_problem-001
#
#   ==========================================
#     Summary
#   ==========================================
#     New rows:     3
#     Extracted:    3
#     Failed:       0
#
#   ✓ Tracking saved: 8 total rows processed

# 3. Process the new submissions
tools/workflow process-new problem-001-dimer-ksat -j 4

# Done! The 3 new submissions are processed, the 5 old ones are skipped.
```

### Weekly Workflow (Even Easier)

```bash
# Monday: Export sheet
# File → Download → CSV → ~/Downloads/week1.csv

# Run extraction
python3 tools/extract_from_sheets.py ~/Downloads/week1.csv problem-001-dimer-ksat

# Process
tools/workflow process-new problem-001-dimer-ksat -j 8

# That's it for the week!
```

---

## 🔄 How Auto-Tracking Works

The script maintains a tracking file that remembers which rows it's already extracted:

```
~/lammps-state-change-private/
└── .sheets-tracking/
    └── responses.json           # Tracks processed rows
```

**`responses.json` content:**
```json
{
  "processed_rows": [2, 3, 4, 5, 6, 7, 8],
  "last_updated": "2026-01-26T14:30:22"
}
```

**How it works:**

1. **First run:** Processes all rows, saves row numbers to `responses.json`
2. **Second run:** Only processes rows NOT in `responses.json`
3. **Third run:** Only processes new rows added since last run
4. **And so on...**

**You never have to manually track which submissions you've already extracted!**

---

## 📊 Example: Real Usage

### Day 1: Launch Competition

```bash
# Create Google Form → share link with participants
# 5 people submit

# Export CSV
python3 tools/extract_from_sheets.py ~/Downloads/responses.csv problem-001-dimer-ksat

# Output:
#   Processing 5 new rows...
#   ✓ Row 2: alice
#   ✓ Row 3: bob
#   ✓ Row 4: charlie
#   ✓ Row 5: diana
#   ✓ Row 6: eve
#   Extracted: 5

# Process
tools/workflow process-new problem-001-dimer-ksat -j 4
```

### Day 3: More Submissions

```bash
# 3 more people submitted (rows 7, 8, 9 in sheet)
# Export same CSV again

python3 tools/extract_from_sheets.py ~/Downloads/responses.csv problem-001-dimer-ksat

# Output:
#   Already processed: 5 rows
#   Processing 3 new rows...
#   ✓ Row 7: frank
#   ✓ Row 8: grace
#   ✓ Row 9: henry
#   Extracted: 3  ← Only the NEW ones!

# Process
tools/workflow process-new problem-001-dimer-ksat -j 4
```

### Day 5: Alice Resubmits (Updated Her Form Response)

```bash
# Alice edited her response (row 2 updated)
# Export CSV

python3 tools/extract_from_sheets.py ~/Downloads/responses.csv problem-001-dimer-ksat

# Output:
#   Already processed: 8 rows
#   Processing 0 new rows...  ← No NEW rows
#   Extracted: 0

# Wait, Alice updated her submission but it's not being extracted!
# Solution: Use --force to reprocess specific rows OR
# Better: Alice should submit a NEW form response (new row)
```

**Important:** If users edit their responses, it updates the existing row. The script won't detect this as "new".

**Solutions:**
1. **Recommended:** Tell users to submit a NEW response (not edit old one)
2. **Alternative:** Use `--force` to reprocess all rows (then delete duplicates manually)

---

## 🛠️ Advanced Usage

### See What Would Be Extracted (Dry Run)

```bash
# Don't actually extract, just show what WOULD be extracted
python3 tools/extract_from_sheets.py ~/Downloads/responses.csv problem-001 --dry-run

# Output:
#   Processing 3 new rows...
#   Would extract row 7: frank
#   Would extract row 8: grace
#   Would extract row 9: henry
#   Dry run complete - no files created
```

### Force Reprocess All Rows

```bash
# Ignore tracking, reprocess everything
python3 tools/extract_from_sheets.py ~/Downloads/responses.csv problem-001 --force

# Output:
#   Processing 8 new rows...
#   [Reprocesses all rows]
#   Extracted: 8
```

**Warning:** This will create duplicate submissions in inbox. Use carefully!

### Reset Tracking (Start Fresh)

```bash
# Delete tracking file to reset
rm ~/lammps-state-change-private/.sheets-tracking/responses.json

# Next run will process all rows again
python3 tools/extract_from_sheets.py ~/Downloads/responses.csv problem-001
```

---

## 📝 Google Form Template (Copy This)

### Form Title
```
Molecular Computing Competition - Submit Your Solution
```

### Form Description
```
Submit your policy design for evaluation!

Instructions:
1. Fill in your username (this will appear on the leaderboard)
2. Select which problem you're solving
3. Copy your ENTIRE policy.json file and paste below
4. (Optional) If the problem allows parameters, paste params.json
5. Enter your email for result notifications

Results will be sent within 48 hours!

Questions? Open an issue: https://github.com/YourOrg/lammps-state-change/issues
```

### Field Templates

**Username:**
```
Type: Short answer
Label: Username
Help text: Choose a username for the leaderboard (e.g., alice)
Required: Yes
```

**Problem ID:**
```
Type: Dropdown
Label: Problem ID
Options:
  - problem-001-dimer-ksat
  - problem-002-octahedron
  - problem-003-hamiltonian
Required: Yes
```

**Policy JSON:**
```
Type: Paragraph
Label: policy.json
Help text: Paste the ENTIRE contents of your policy.json file here.
           Must be valid JSON (check with jsonlint.com if unsure).

           Example format:
           {
             "policy_name": "my_awesome_policy",
             "check_frequency": 100,
             "state_transitions": [...]
           }
Required: Yes
```

**Parameters JSON:**
```
Type: Paragraph
Label: params.json (optional)
Help text: Only if this problem allows custom parameters.
           Leave blank otherwise.

           Example:
           {
             "morse_depth_AB": 1.2,
             "morse_alpha": 5.0
           }
Required: No
```

**Email:**
```
Type: Short answer
Label: Email
Help text: For result notifications
Validation: Email address
Required: Yes
```

### Confirmation Message
```
Thank you for your submission!

We've received your solution and it will be evaluated within 48 hours.

You'll receive an email with:
- Your scores
- Your leaderboard position
- Links to the public leaderboard

Want to improve your solution? You can submit again anytime!

Good luck! 🚀
```

---

## ✅ Advantages Over Email

| Feature | Email | Google Sheets |
|---------|-------|---------------|
| **Participant ease** | Need to attach files | Just paste JSON |
| **Your extraction time** | 2 min per submission | 10 sec total |
| **Auto-tracking** | Manual | Automatic |
| **Multiple problems** | Separate emails | One form with dropdown |
| **Validation** | None | Can add form validation |
| **History** | Email archive | Spreadsheet |
| **Bulk export** | One by one | Export all at once |
| **Resubmissions** | New email | Edit response or new row |

---

## 🔧 Column Name Customization

If you want to use different column names in your Google Form, edit the script:

**File:** `tools/extract_from_sheets.py`

**Line ~58:**
```python
# Current (default):
username = row.get('Username', row.get('username', '')).strip()
problem_id_from_form = row.get('Problem ID', row.get('problem_id', problem_id)).strip()
policy_json_str = row.get('policy.json', row.get('Policy JSON', '')).strip()

# Customize to match your form:
username = row.get('Your Name', '').strip()  # If you named it "Your Name"
problem_id_from_form = row.get('Which Problem?', problem_id).strip()
policy_json_str = row.get('Paste Policy Here', '').strip()
```

---

## 🚀 Quick Reference

```bash
# Export Google Sheet: File → Download → CSV

# Extract only new submissions
python3 tools/extract_from_sheets.py responses.csv problem-001-dimer-ksat

# Process
tools/workflow process-new problem-001-dimer-ksat -j 4

# That's it!
```

---

## ❓ FAQ

**Q: What if someone submits twice?**
A: If they submit a NEW response (not edit), you'll get 2 rows. Both will be processed. Leaderboard will show their best score automatically.

**Q: What if someone edits their response?**
A: The existing row is updated, but the script won't detect it as "new". Better to have them submit a new response.

**Q: Can I use this for multiple problems?**
A: Yes! Either use one form with a "Problem ID" dropdown, or create separate forms (and CSVs) per problem.

**Q: What if extraction fails for a row?**
A: The script shows error details and continues to next row. That row is NOT marked as processed, so next run will retry it.

**Q: Can I delete the tracking file?**
A: Yes, but next run will reprocess ALL rows (creating duplicates in inbox).

**Q: How do I know which submissions are new without running the script?**
A: Compare total rows in CSV with number in `~/.sheets-tracking/responses.json`

---

## 🎉 Summary

**Setup once:** Create Google Form (5 minutes)
**Share:** Give participants the form link
**Process:** Export CSV + run script (2 commands, 1 minute total)
**Automatic:** Script tracks which rows already processed

**No more manual email extraction!** 🚀
