# Module0 Analysis Toolkit

This repository contains a small ROOT-based program that reads hit-level data from a tree, organizes the hits per detector layer and eta side, clusters nearby hits, and reports simple efficiency numbers. The code is intentionally split into a few focused components so you can tweak the behaviour without having to read an entire monolithic source file.

## Build and Run

`compile.sh` shows the exact build/run command:

```bash
g++ -o plot_eff eff_plot_from_hit.cxx cluster_utils.cxx `root-config --cflags --libs`
./plot_eff your_file.root
```

Adjust the root file glob as needed (the script currently runs over `*.root`). The executable expects the standard branches (`nHits`, `hit_channel`, `hit_time1`, `hit_time2`, `hit_rise`, `hit_bcid`) to exist inside a tree called `tree`.

## Code Layout

* `eff_plot_from_hit.cxx` — The main analysis program. It owns the ROOT I/O, event loop, trigger handling, layer bookkeeping, and final reporting (hit summaries and efficiencies). Edit here when you need to:
  * Change which branches are read or how events are filtered.
  * Modify the trigger definition or the summary output.
  * Adjust how efficiencies are computed.

* `layer.h` — Defines the `Layer` class that stores hits separately for eta1/eta2, enforces "first leading edge only", and computes the time-over-threshold including the special wrapping correction. Edit here if you need to:
  * Capture additional per-hit quantities.
  * Change how TOT is calculated or corrected.

* `cluster.h` — Defines the lightweight `Cluster` struct (channels, earliest time, TOT from that earliest hit, center channel, etc.). Extend this to store more derived quantities per cluster.

* `cluster_utils.h/.cxx` — Holds `LayerClusters` and the `buildClusters` helper that sorts hits, groups neighbours (±1 channel and ≤5 ns apart), and selects the cluster timing/TOT. Update this file when you want to tweak the clustering criteria or how clusters are constructed.

## Typical Edits

* **New detector layers or eta logic** → Update `eff_plot_from_hit.cxx` where layers are instantiated and populated.
* **Different clustering rules** → Modify `cluster_utils.cxx::buildClusters`.
* **Alternative TOT handling** → Adjust the correction logic inside `Layer::addHit` (in `layer.h`).
* **Additional summary metrics** → Extend the per-event printout or efficiency block at the end of `eff_plot_from_hit.cxx`.

Keeping the responsibilities separated this way should make it straightforward to evolve the analysis without breaking unrelated pieces. If you add new source files, remember to list them in `compile.sh`.