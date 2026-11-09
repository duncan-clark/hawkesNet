# hawkesNet 0.2.0

* Added a small-network birth/turnover reference API with separate initial- and later-topology coefficients, likelihood fitting and a feedback information check. A generic Poisson-proposal union interface supports arbitrary caller-defined edit groups.

* Added the CRAN complete-update likelihood, fitting and forward-simulation interfaces. Release and edit rates can share a two-component update mixture; the paired law sums over single-edit and migration-pair proposal configurations.
* Added shared maintainer, author and contributor effects to the primary CS and TERGM specifications. The manuscript's transitive-triple, open-diamond and shared-consumer motifs are reserved for goodness-of-fit checks.
* Added calendar timing with separate fast/slow release-type effects and separate initial- and later-topology coefficients, with clipped later changes.
* Reworked onboarding and small runnable examples. The documentation distinguishes the Poisson-set proposal construction from the historical direct size-conditional CS law and identifies the scope of addition-only growth interfaces.
* Corrected adaptive thinning for inhomogeneous backgrounds, observation-window KDE normalization, tied-time handling and experimental dyadic score/amount calculations. Added input validation for native CRAN data and explicit macOS/Linux source support.
* Retired superseded research scripts and local development files. The installed CRAN workflow includes the guide and primary refit script. Research scripts used by paper diagnostics retain their checkout paths.

## Earlier development work (September 22, 2026)

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
