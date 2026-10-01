# Release Checklist

This checklist must be completely validated before final tagging (`v1.0`) and release.

## Integrity and Pre-Registration
- [ ] PROTOCOL.md frozen at Gate 1 with git tag `protocol-v1` prior to model fitting
- [ ] All deviations logged in DEVIATIONS.md with dates and rationales
- [ ] No fabricated, simulated, or hardcoded values in `results/` or manuscript
- [ ] No causal language used; banned terms absent: "state-of-the-art", "robust", "clinical utility", "novel", "validated"
- [ ] Explicit label included on all documents: "preliminary, preprint in preparation"

## Code and Verification
- [ ] Pinned `requirements.txt` installs cleanly in clean venv
- [ ] `pytest` passes with 0 failures or warnings
- [ ] CI workflow runs and passes on GitHub Actions
- [ ] `python scripts/run_all.py` executes end-to-end deterministically
- [ ] Clean-room reproduction verified on fresh clone and environment within numerical tolerance

## Data and Manifests
- [ ] Raw data gitignored; no raw data files in repository
- [ ] `results/data_manifest.json` committed with SHA256 hashes, source URLs, access dates
- [ ] No file > 5MB committed

## Claims and Reporting
- [ ] Every numeric claim in README, manuscript, and abstract traced in CLAIMS.md
- [ ] TRIPOD-AI checklist completed honestly in `checklists/TRIPOD-AI.md`
- [ ] All citations marked `[CITATION TO VERIFY]` verified against peer-reviewed records
- [ ] LaTeX manuscript compiles cleanly with 0 undefined citations or missing references
- [ ] `paper_overleaf.zip` generated and verified
