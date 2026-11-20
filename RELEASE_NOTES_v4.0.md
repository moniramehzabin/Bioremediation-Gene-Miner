# Bioremediation Gene Miner v4.0

Release date: 2026-11-20

- Add a Windows launcher for selecting annotated GenBank files and running Gene Miner.
- Add foreground InterPro recovery with saved job IDs, cached results, bounded read retries and individual failure handling.
- Support InterProScan v6 REST submission and validated result downloads.
- Include the 413-sequence reference snapshot and matching annotation rules.
- Add PET assay-control evidence and Review routing for exact assay-negative PET control matches.

Requires Python 3 with Tkinter and a Windows DIAMOND executable.
Open windows_launcher/Run Gene Miner.cmd for a new analysis.
Open windows_launcher/Resume InterPro.cmd to continue saved InterPro scans.

This release remains a research prototype.
