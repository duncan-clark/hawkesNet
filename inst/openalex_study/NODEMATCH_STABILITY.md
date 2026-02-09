# Why the nodeMatch fit can be unstable

Summary of likely causes and things to try.

## 1. **Weak identification when gender is mostly "unknown"**

- If most or all nodes have `gender == "unknown"`, the **nodeMatch** statistic has almost no variation (nearly all edges are unknown–unknown).
- The **vertex_categorical** parameters (multinomial P(female), P(male); P(unknown) = 1 − sum) are identified mainly by the observed gender of *new* nodes at each event. If almost all observed values are "unknown", the likelihood is barely sensitive to the female/male probabilities, so the optimizer can drift or hit boundaries and look "unstable" or "scrambled".
- **Mitigation:** Ensure gender has enough non-unknown variation (e.g. run gender prediction so many nodes are female/male), or **skip the nodeMatch step** when e.g. >90% of nodes are "unknown".

## 2. **n−1 parametrization and boundaries**

- `vertex_categorical` uses an n−1 parametrization: two free parameters (e.g. female, male) in (0, 1) with sum < 1 (unknown = reference). They are bounded in `[eps, 1−eps]`.
- Nelder–Mead can push one parameter to the boundary; the other then becomes poorly determined. `repair_vertex_categorical_params` then clamps and can make displayed estimates look odd.
- **Mitigation:** Reparametrize to an unconstrained scale (e.g. log-odds for female and male vs unknown) so the optimizer never sits on the boundary; transform back for display and simulation. Alternatively, try **L-BFGS-B** with the existing bounds so the optimizer respects constraints without collapsing the simplex.

## 3. **Scale and parscale**

- Structural **CS_params** and **vertex_categorical** are on different scales. If the log-likelihood is much more sensitive to CS_params than to vertex_categorical, Nelder–Mead may make many steps that change only structure, and then large relative steps in vertex_categorical can look erratic.
- **Mitigation:** Tune **parscale** so that typical steps in vertex_categorical (e.g. 0.1) are comparable in *effect on the log-likelihood* to typical steps in CS_params. You can try slightly larger parscale for the vertex_categorical entries (e.g. 0.2) so the simplex moves more in that direction.

## 4. **Initialization**

- The same initial **vertex_categorical** (e.g. female=0.1, male=0.5) is used every time. If the data favor very different proportions (e.g. almost all unknown, or 0.4/0.4), the simplex may contract in a bad direction before finding a good region.
- **Mitigation:** Initialize from **observed proportions** of female/male/unknown among nodes (or among new nodes at event times), with a small perturbation. Or run **multiple random starts** (e.g. perturb female/male around 0.33/0.33) and keep the fit with highest log-likelihood.

## 5. **Flat likelihood in the nodeMatch direction**

- If the structural terms already fit the data well and the **nodeMatch** coefficient is small, the likelihood can be flat in the (nodeMatch coefficient, vertex_categorical) direction. Then the optimizer may not converge to a unique mode and estimates can vary across runs.
- **Mitigation:** Fix **vertex_categorical** at observed (or pre-estimated) proportions and fit only the nodeMatch coefficient; if that is stable, then free vertex_categorical with that start. Alternatively, add a **weak prior** or penalty to keep vertex_categorical near observed proportions.

## 6. **Nelder–Mead and iteration limit**

- Nelder–Mead does not use gradients; with flat or boundary-sensitive regions the simplex can shrink slowly or wander. With a tight **MAX_ITER** the fit may stop before settling.
- **Mitigation:** Increase **MAX_ITER** for the nodeMatch fit (e.g. 5000 → 8000), or try **L-BFGS-B** with the current bounds for the nodeMatch run only.

## Practical checklist

- [ ] Check gender distribution: if almost all "unknown", skip nodeMatch or run gender prediction first.
- [ ] Initialize vertex_categorical from observed node proportions (with small perturbation).
- [ ] Tune parscale for vertex_categorical (e.g. 0.15–0.2) so it is not under-scaled relative to CS_params.
- [ ] Consider L-BFGS-B for the nodeMatch fit (bounds already defined in `build_optim_bounds`).
- [ ] Increase MAX_ITER for nodeMatch.
- [ ] Optional: fix vertex_categorical at observed proportions and fit only the nodeMatch coefficient to check stability.
