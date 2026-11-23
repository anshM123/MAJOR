# Running the extra Monte Carlo jobs on your own machine (Windows 11, via WSL2)

These jobs compute **direct crystal points**: the same physics code as `ads4/code/spinfunc/hc_points.py`, unchanged. They are independent of
each other, so they parallelise across cores. Each one needs **≈ 0.7 GB RAM** (measured peak 0.66 GB for a full-size point) and takes a
few minutes on one core. They read **no** cached data, so nothing beyond the git repo is needed.

## 1. One-time setup (PowerShell as Administrator, then the Ubuntu terminal)
```
wsl --install -d Ubuntu
```
Reboot if asked, open "Ubuntu" from the Start menu, then:
```
sudo apt update && sudo apt install -y python3-venv git
python3 -m venv ~/venv && source ~/venv/bin/activate
pip install numpy scipy mpmath
git clone https://github.com/anshM123/MAJOR.git && cd MAJOR && git checkout claude/eloquent-cori-nftfjy
```
(The code uses Linux-only file locking, so run it inside WSL, not in native Windows Python.)

## 2. Run
```
source ~/venv/bin/activate && cd ~/MAJOR/ads4/ext
python3 run_ext.py 12
```
`12` = number of parallel workers. Rule of thumb: **workers ≈ (free RAM in GB) / 0.8**, at most (cores − 2). With 10 GB to spare, use 12.
WSL2 by default gets half of the machine's RAM; to give it more, create `C:\Users\<you>\.wslconfig` with
`[wsl2]` / `memory=16GB` / `processors=22` and run `wsl --shutdown`.

The run can be stopped (Ctrl-C) and restarted any time; finished points are skipped. There are 217 jobs; expect roughly 1–3 hours with
12 workers.

## 3. Send the results back
```
python3 pack.py
git add results/ext_results.json && git commit -m "ext: direct points from external run" && git push
```
Only that one JSON file is committed (the raw `out/` directory is git-ignored).

## What the jobs are (`jobs.txt`)
* **A** — direct bulk rows (34 nodes, x = 0.05–0.775) for the near-gap states J = 7, 23, 27, 31 at Δ_g = 1000, the ones still at
  model level; they get the same treatment that moved J = 12, 16, 21, 29 to error-budgeted.
* **B** — J = 29 small-ν window (ν = 10–17.5, n3 contour), 16 extra independent seeds, to close its margin without the signed window.
* **C** — exploratory: J = 400 at its twist edge, in the near-gap band 345 ≤ J < 857 that no theorem covers.
