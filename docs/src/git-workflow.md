# Git workflow in VS Code

This manual describes a practical workflow for contributing to AeroBeams_UVLM.jl with Git and VS Code. It assumes that you are working on Windows and that Julia 1.11.5 is installed.

The repository is hosted at [github.com/jptpsantos/AeroBeams_UVLM.jl](https://github.com/jptpsantos/AeroBeams_UVLM.jl).

## 1. The basic idea

Git records the history of a project. It lets you make changes, review them, return to an earlier state, and share work with other people.

There are four places to keep in mind:

- **Working tree:** the files currently open in your VS Code folder.
- **Staging area:** the changes selected for the next commit.
- **Local repository:** the commits stored on your computer.
- **Remote repository:** the shared copy on GitHub, usually called `origin`.

A **commit** is a named checkpoint. A **branch** is an independent line of commits. A **pull request** asks other people to review and merge a branch into another branch.

The usual direction of work is:

```text
GitHub (remote) -> your local branch -> edit files -> stage -> commit -> push -> GitHub
```

## 2. One-time setup

### Install the tools

Install:

- Git for Windows: <https://git-scm.com/download/win>
- Julia: version 1.9 or newer below 2.0 is compatible with this project; Julia 1.11.5 is expected to work.
- The VS Code **Julia** extension.

Restart VS Code after installing Git so that the integrated terminal can find the `git` command.

### Check the installations

Open the VS Code terminal with `Ctrl+` and run:

```powershell
git --version
julia --version
```

If `git` is not recognized, close and reopen VS Code. If it still is not recognized, install Git for Windows or add its installation directory to your Windows PATH.

### Configure your Git identity

Git stores your name and email in each commit. Use the email associated with your GitHub account when possible:

```powershell
git config --global user.name "Your Name"
git config --global user.email "your.email@example.com"
```

Check the configuration:

```powershell
git config --global --list
```

This does not authenticate you to GitHub. When GitHub asks for credentials, use GitHub's browser sign-in, a credential manager, or an SSH key. Never put a password or access token inside a repository file or a command that will be saved in shell history.

## 3. Open the project correctly

In VS Code, use **File > Open Folder** and open the folder that contains `Project.toml`, `src`, and `test`:

```text
AeroBeams_UVLM.jl
```

Opening a parent folder or only the `src` folder can make Git and Julia commands behave unexpectedly.

The VS Code terminal should start in the repository root. Verify it with:

```powershell
Get-Location
Test-Path .git
Test-Path .\Project.toml
```

The first command prints the current folder. The two `Test-Path` commands should print `True`.

## 4. Understand the repository state

The most useful Git command is:

```powershell
git status
```

It tells you:

- your current branch;
- whether it is ahead of or behind GitHub;
- which files are modified;
- which files are staged;
- which files are untracked.

A compact version is:

```powershell
git status --short --branch
```

The two columns in the short output describe the file state. The first column is the staging area and the second is the working tree. For example:

```text
 M src/Beam.jl
M  test/staticStructuralTests.jl
?? test/myNewTest.jl
```

These mean:

- ` M`: modified but not staged;
- `M `: modified and staged;
- `??`: a new file not yet tracked by Git.

## 5. Start every change safely

The shared branch in this repository is `main`. Keep it stable and create a branch for each task.

First update your local `main`:

```powershell
git switch main
git pull --ff-only origin main
```

Then create a task branch:

```powershell
git switch -c fix/short-description
```

Examples of branch names:

```text
fix/cantilever-boundary-condition
feature/new-gust-model
docs/git-workflow
```

Confirm the branch:

```powershell
git branch --show-current
```

A branch is cheap to create. Use one even for a small change because it keeps unfinished work away from `main`.

## 6. Make a change in this project

The main code is in `src/`. Tests are in `test/`. Examples and research scripts are in `dev/` and `test/examples/`.

Important project files include:

- `Project.toml`: package name, dependencies, and compatibility rules;
- `src/`: package implementation;
- `test/runtests.jl`: test entry point;
- `test/`: structural, aerodynamic, and aeroelastic tests;
- `docs/`: documentation source.

Before running Julia code, activate the repository environment:

```powershell
julia --project=.
```

Inside the Julia prompt, install or update the dependencies recorded for this environment:

```julia
using Pkg
Pkg.instantiate()
```

You can then load the package:

```julia
using AeroBeams
```

For an interactive development session, `Revise` can reload source changes after you edit them:

```julia
using Revise
using AeroBeams
```

Do not edit files inside Julia's global package depot to change this repository. Make changes in this checkout, normally under `src/` or `test/`.

## 7. Test before committing

Run the complete package test suite from the repository root:

```powershell
julia --project=. -e "using Pkg; Pkg.test()"
```

A successful run should finish without test errors. Tests may take some time because this package includes structural and aerodynamic analyses.

When working on one area, you may run a relevant test file directly:

```powershell
julia --project=. test/staticStructuralTests.jl
```

Use the complete `Pkg.test()` command before opening a pull request. If a test fails, save the failure output, inspect the first relevant error, and do not commit generated output just because it was produced during the run.

After changing dependencies in `Project.toml`, run:

```powershell
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

Review dependency changes carefully. `Manifest.toml`, if present, records the resolved versions and may also change.

## 8. Review your changes

Check the files that changed:

```powershell
git status
git diff
```

`git diff` shows unstaged changes. To inspect staged changes:

```powershell
git diff --cached
```

Useful history commands are:

```powershell
git log --oneline --decorate -10
git show --stat HEAD
git diff main...HEAD
```

Before staging, check for accidental changes such as debug prints, large output files, credentials, or unrelated formatting.

## 9. Stage and commit

Stage only the files belonging to this task:

```powershell
git add src/Beam.jl test/staticStructuralTests.jl
```

For a new documentation file:

```powershell
git add docs/src/git-workflow.md docs/src/index.md
```

Review exactly what will be committed:

```powershell
git status
git diff --cached
```

Create a commit:

```powershell
git commit -m "Add Git workflow documentation"
```

Good commit messages are short, specific, and written as an action. Examples:

```text
Fix beam boundary condition assembly
Add regression test for gust response
Document Julia development workflow
```

A commit should normally contain one coherent change. If you staged something by mistake, remove it from the staging area without deleting your file:

```powershell
git restore --staged path/to/file.jl
```

## 10. Push your branch and open a pull request

The first push connects your local branch to a branch on GitHub:

```powershell
git push -u origin fix/short-description
```

After that, future pushes only need:

```powershell
git push
```

Open the repository on GitHub and create a pull request from your branch into `main`. Include:

- what changed;
- why it changed;
- which tests you ran;
- any known limitations or follow-up work.

Do not work directly on `main` unless the project maintainer specifically asks you to.

## 11. Keep a branch up to date

If `main` receives new commits while your pull request is open, update your branch:

```powershell
git fetch origin
git switch main
git pull --ff-only origin main
git switch fix/short-description
git merge main
```

If there are no conflicts, run the tests again and push:

```powershell
julia --project=. -e "using Pkg; Pkg.test()"
git push
```

`git fetch` downloads information without changing your files. `git pull` downloads and integrates changes into the currently checked-out branch.

## 12. Resolve a merge conflict

A conflict means Git could not automatically combine two edits. Check the affected files:

```powershell
git status
```

Open each conflicted file in VS Code. Conflict markers look like this:

```text
<<<<<<< HEAD
Your branch's version
=======
The other branch's version
>>>>>>> main
```

Choose or combine the correct code, then delete all conflict markers. Run the relevant tests. Mark the resolved file:

```powershell
git add path/to/resolved-file.jl
```

When every conflict is staged, finish the merge:

```powershell
git commit
```

If you started the merge and want to cancel it before committing:

```powershell
git merge --abort
```

Do not guess when resolving scientific or numerical code. Compare both versions and rerun the tests that exercise the affected behavior.

## 13. Undo common mistakes

### Edit is not staged

Nothing special is required. Continue editing or discard it deliberately.

### Staged the wrong file

```powershell
git restore --staged path/to/file
```

### Want to discard local edits in one file

This permanently removes uncommitted changes in that file:

```powershell
git restore path/to/file
```

Check `git diff` first. Do not use this command on work you may need later.

### Want to change the last commit

If the commit has not been pushed:

```powershell
git add path/to/file
git commit --amend
```

Avoid amending a commit that other people have already based work on.

### Need to pause unfinished work

You can make a temporary local commit, or use a stash:

```powershell
git stash push -m "unfinished beam change"
git stash list
git stash pop
```

Inspect the result after `git stash pop`; it can also produce conflicts.

Never use `git reset --hard` or force-push unless you understand exactly which commits and files will be lost and a maintainer has agreed when shared history is involved.

## 14. VS Code interface

The Source Control panel is opened with `Ctrl+Shift+G`.

- The branch selector in the lower-left corner creates or switches branches.
- The Source Control panel lists modified and untracked files.
- Clicking a file opens a side-by-side diff.
- The `+` button stages a file.
- The checkmark creates a commit from staged files.
- The sync or push control sends your branch to GitHub.
- The three-dot menu contains fetch, pull, branch, stash, and other commands.

The interface and terminal run the same Git operations. Use whichever makes the operation easiest to understand, and check the Source Control diff before committing.

## 15. A complete daily example

```powershell
# Begin from the latest shared code
git switch main
git pull --ff-only origin main

# Create an isolated branch
git switch -c fix/improve-test-message

# Edit files in VS Code

# Activate the Julia project and test
julia --project=. -e "using Pkg; Pkg.test()"

# Review and commit
git diff
git status
git add path/to/changed-file.jl
git diff --cached
git commit -m "Improve test message"

# Share the branch
git push -u origin fix/improve-test-message
```

## 16. Quick reference

```powershell
git status                         # Show current state
git branch --show-current          # Show current branch
git switch main                    # Switch branches
git switch -c feature/name         # Create and switch to a branch
git pull --ff-only origin main     # Update main safely
git fetch origin                   # Download remote information
git diff                           # Review unstaged changes
git add path/to/file               # Stage a file
git restore --staged path/to/file  # Unstage a file
git diff --cached                  # Review staged changes
git commit -m "Message"           # Create a commit
git push                           # Upload commits
git log --oneline -10              # View recent history

julia --project=. -e "using Pkg; Pkg.test()"  # Run package tests
```

When uncertain, stop before committing and inspect `git status`, `git diff`, and `git diff --cached`. Those three commands usually make the next safe step clear.
