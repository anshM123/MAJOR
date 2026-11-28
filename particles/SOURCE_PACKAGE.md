# Point-particle research source package

Research checkpoint: 28 November 2026.

The `particles/` directory contains all new proofs, source audits, five numerical/symbolic checkers, the 372-facet catalogue, archived validation outputs, and publication-ready figures. Start with `particles/README.md`; its status table distinguishes physical HRT theorems from mathematical functionals and states the open collision realization problem.

The accompanying `particle_cone_research.pdf` is the complete 70-page report. `particle_cone_research.tex` uses relative figure paths and can be rebuilt from this directory with a standard TeX installation. The separate `particle_cone_extension.patch` adds only `particles/` to the upstream repository at the base revision recorded in `particles/SOURCE_MANIFEST.json`; it leaves previous `strips/` work untouched.

From `particles/`, install `requirements.txt` and run the commands in its README. The checkers contact no remote services. No downloaded research papers are bundled. The results rely on the repository's standing vacuum theorems; this package is not an independent audit of that entire foundation.
