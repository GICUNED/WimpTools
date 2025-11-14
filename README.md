# WimpTools

## Automated Repository Maintenance

This repository includes an automated maintenance workflow that helps keep the repository clean by removing inactive and merged branches.

### How It Works

The maintenance workflow runs weekly (every Sunday at 3 AM UTC) and can also be triggered manually. It performs the following actions:

1. **Identifies merged branches**: Detects branches that have been merged into the default branch via a pull request
2. **Identifies stale branches**: Finds branches with no commits in the last 90 days (configurable)
3. **Protects active branches**: Skips branches that:
   - Are protected (main, master, develop by default)
   - Have open pull requests
   - Are the repository's default branch

### Default Configuration

- **DRY_RUN**: `true` (safe mode - only reports what would be deleted, no actual deletions)
- **MAINTENANCE_DAYS**: `90` (days of inactivity before considering a branch stale)
- **PROTECTED_BRANCHES**: `main,master,develop`

### Manual Execution

You can manually trigger the maintenance workflow from the Actions tab:

1. Go to the "Actions" tab in the GitHub repository
2. Select "Git Repository Maintenance" workflow
3. Click "Run workflow"
4. Optionally configure:
   - **Dry run mode**: Keep as `true` to preview changes without deleting
   - **Maintenance days**: Adjust the inactivity threshold
   - **Protected branches**: Modify the list of branches to protect

### Testing the Workflow

Before enabling actual deletions, it's recommended to:

1. Run the workflow in dry-run mode (default) to see what would be deleted
2. Review the workflow logs to verify the branches identified for deletion
3. If satisfied, run manually with `dry_run: false` to perform actual cleanup
4. Monitor the results before relying on the scheduled automatic runs

### Precautions

⚠️ **Important Notes**:

- The workflow does NOT handle branches with rewritten history or force pushes differently
- Always review dry-run logs before disabling dry-run mode
- Consider protecting any long-lived feature branches by adding them to PROTECTED_BRANCHES
- The workflow uses GITHUB_TOKEN for authentication, which has write access to the repository

### Workflow File Location

The maintenance workflow is defined in `.github/workflows/git-maintenance.yml` and uses the Python script `scripts/maintain_repo.py`.
