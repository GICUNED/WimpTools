#!/usr/bin/env python3
"""
Git Repository Maintenance Script

This script automates the cleanup of inactive and merged branches in a GitHub repository.
It identifies branches that are either:
1. Merged into the default branch via a pull request
2. Stale (no commits in the last MAINTENANCE_DAYS days)

The script will NOT delete branches that:
- Are protected (listed in PROTECTED_BRANCHES)
- Are the repository's default branch
- Have open pull requests

Environment Variables:
- GITHUB_TOKEN: GitHub authentication token (required)
- GITHUB_REPOSITORY: Repository in 'owner/repo' format (required)
- DRY_RUN: If 'true', only report what would be deleted (default: 'true')
- MAINTENANCE_DAYS: Days of inactivity before considering a branch stale (default: 90)
- PROTECTED_BRANCHES: Comma-separated list of branches to protect (default: 'main,master,develop')
"""

import os
import sys
from datetime import datetime, timedelta, timezone
from github import Github, GithubException
from dateutil import parser as date_parser


def log(message, level="INFO"):
    """Print a timestamped log message."""
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    print(f"[{timestamp}] [{level}] {message}")


def get_env_var(name, default=None, required=False):
    """Get an environment variable with optional default and required check."""
    value = os.environ.get(name, default)
    if required and not value:
        log(f"Missing required environment variable: {name}", "ERROR")
        sys.exit(1)
    return value


def parse_bool(value):
    """Parse a boolean value from string."""
    if isinstance(value, bool):
        return value
    return str(value).lower() in ('true', '1', 'yes', 'on')


def main():
    # Read configuration from environment
    github_token = get_env_var("GITHUB_TOKEN", required=True)
    repository_name = get_env_var("GITHUB_REPOSITORY", required=True)
    dry_run = parse_bool(get_env_var("DRY_RUN", "true"))
    maintenance_days = int(get_env_var("MAINTENANCE_DAYS", "90"))
    protected_branches_str = get_env_var("PROTECTED_BRANCHES", "main,master,develop")
    
    # Parse protected branches
    protected_branches = set(b.strip() for b in protected_branches_str.split(",") if b.strip())
    
    log("=" * 60)
    log("Git Repository Maintenance Script")
    log("=" * 60)
    log(f"Repository: {repository_name}")
    log(f"Dry Run: {dry_run}")
    log(f"Maintenance Days: {maintenance_days}")
    log(f"Protected Branches: {', '.join(sorted(protected_branches))}")
    log("=" * 60)
    
    # Initialize GitHub client
    try:
        gh = Github(github_token)
        repo = gh.get_repo(repository_name)
        log(f"Connected to repository: {repo.full_name}")
    except GithubException as e:
        log(f"Failed to connect to repository: {e}", "ERROR")
        sys.exit(1)
    
    # Get default branch
    default_branch = repo.default_branch
    log(f"Default branch: {default_branch}")
    
    # Add default branch to protected branches
    protected_branches.add(default_branch)
    
    # Calculate cutoff date for stale branches
    cutoff_date = datetime.now(timezone.utc) - timedelta(days=maintenance_days)
    log(f"Stale cutoff date: {cutoff_date.strftime('%Y-%m-%d %H:%M:%S UTC')}")
    log("=" * 60)
    
    # Get repository owner for PR filtering
    owner = repository_name.split('/')[0]
    
    # Counters for summary
    total_branches = 0
    skipped_protected = 0
    skipped_active_prs = 0
    merged_branches = 0
    stale_branches = 0
    deleted_branches = 0
    errors = 0
    
    # Iterate through all branches
    try:
        branches = repo.get_branches()
        log("Analyzing branches...")
        log("")
        
        for branch in branches:
            total_branches += 1
            branch_name = branch.name
            
            # Skip protected branches
            if branch_name in protected_branches:
                log(f"⊗ {branch_name}: Protected branch, skipping")
                skipped_protected += 1
                continue
            
            try:
                # Check for open pull requests with this branch as head
                open_prs = list(repo.get_pulls(state='open', head=f"{owner}:{branch_name}"))
                if open_prs:
                    log(f"⊗ {branch_name}: Has {len(open_prs)} open PR(s), skipping")
                    skipped_active_prs += 1
                    continue
                
                # Check if branch was merged via a pull request
                closed_prs = list(repo.get_pulls(
                    state='closed',
                    head=f"{owner}:{branch_name}",
                    base=default_branch
                ))
                
                is_merged = False
                for pr in closed_prs:
                    if pr.merged:
                        is_merged = True
                        log(f"✓ {branch_name}: Merged via PR #{pr.number} on {pr.merged_at.strftime('%Y-%m-%d')}")
                        merged_branches += 1
                        break
                
                if not is_merged:
                    # Check if branch is stale (no recent commits)
                    commit = branch.commit
                    # Get commit date from commit author date
                    if commit.commit.author and commit.commit.author.date:
                        commit_date = commit.commit.author.date
                        # Ensure commit_date is timezone-aware
                        if commit_date.tzinfo is None:
                            commit_date = commit_date.replace(tzinfo=timezone.utc)
                        
                        if commit_date < cutoff_date:
                            days_old = (datetime.now(timezone.utc) - commit_date).days
                            log(f"✓ {branch_name}: Stale ({days_old} days old, last commit: {commit_date.strftime('%Y-%m-%d')})")
                            stale_branches += 1
                            is_merged = True  # Mark for deletion
                        else:
                            log(f"⊗ {branch_name}: Active (last commit {commit_date.strftime('%Y-%m-%d')})")
                            continue
                    else:
                        log(f"⊗ {branch_name}: Cannot determine commit date, skipping")
                        continue
                
                # Delete the branch if not in dry run mode
                if is_merged:
                    if dry_run:
                        log(f"  [DRY RUN] Would delete branch: {branch_name}")
                    else:
                        try:
                            ref = repo.get_git_ref(f"heads/{branch_name}")
                            ref.delete()
                            log(f"  ✗ DELETED: {branch_name}", "WARNING")
                            deleted_branches += 1
                        except GithubException as e:
                            log(f"  Failed to delete {branch_name}: {e}", "ERROR")
                            errors += 1
            
            except GithubException as e:
                log(f"  Error processing {branch_name}: {e}", "ERROR")
                errors += 1
                continue
        
        log("")
        log("=" * 60)
        log("SUMMARY")
        log("=" * 60)
        log(f"Total branches analyzed: {total_branches}")
        log(f"Protected branches (skipped): {skipped_protected}")
        log(f"Branches with open PRs (skipped): {skipped_active_prs}")
        log(f"Merged branches identified: {merged_branches}")
        log(f"Stale branches identified: {stale_branches}")
        
        if dry_run:
            log(f"Branches that WOULD be deleted: {merged_branches + stale_branches}")
            log("")
            log("*** DRY RUN MODE - No branches were actually deleted ***")
        else:
            log(f"Branches deleted: {deleted_branches}")
            if errors > 0:
                log(f"Errors encountered: {errors}", "WARNING")
        
        log("=" * 60)
        log("Maintenance completed successfully")
        
    except GithubException as e:
        log(f"Failed to retrieve branches: {e}", "ERROR")
        sys.exit(1)


if __name__ == "__main__":
    main()
