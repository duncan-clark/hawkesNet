# Development: September 22, 2026

* Added explicit CS-1 (`cs_mode="independent"`) and CS-2 (`cs_mode="size_conditional"`) mark laws with matched simulation and cached likelihoods. The earlier all-subset tilt remains `cs_mode="joint"`.
* CS-1 uses frozen logistic single-edge change-statistic weights and exact Poisson collapsing, without subset enumeration. CS-2 draws a truncated-Poisson edge count and enumerates only its size class using whole-update statistics.
* Added a fitting guard for CS-2's unidentified edge-count coefficient, stable normalization and exact nonempty conditioning, invalid-update checks, and exact categorical support at probability zero/one.
* Added runnable examples and matched-law local smoke checks. The previous 5,000-event CS results concern the joint compatibility law.

# hawkesNet 0.1.0

* Initial release of `hawkesNet`, a package for fitting and simulating Hawkes processes on growing networks.
* Support for Barabasi-Albert (BA) and Change Statistics (CS/ERNM) mark PMFs.
* Support for homogeneous and inhomogeneous background rates (KDE-based).
* Goodness-of-fit (GOF) diagnostics for network statistics (degree, ESP, geodesic distance, waiting times).
* Parallelization support for intensity cache building and GOF simulations.
* Visualization tools for network growth and KDE background estimation.
