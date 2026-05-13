# Formula references

Every physics formula implemented in this project, with its provenance.

**Primary source** for the IC construction and integrator is

> Naab, T. (2006). *N-body Simulations and Galaxy Formation* (practicum
> manual, version 0.0, Munich, 11.01.2006). Hereafter **"manual"**.
> Local copy: `../project_3_ideas/nbody_manual.pdf`.

Equation numbers prefixed `2.x` and `3.x` refer to the manual.

Literature citations match the manual's bibliography unless noted.
Where I had to *derive* a formula (because the manual uses it without
spelling out the derivation), the derivation is shown below.

Units throughout: `G = 1`, in the system of Hernquist (1992, 1993b) defined
on manual page 18 (h = 1, M_d = 1, M_h = 5.8, M_b = 1/3, z_0 = 0.2,
γ = 1, r_c = 10, a = 0.1).

---

## 1. `forces.py`

### `gravitational_accelerations`

```
a_i = sum_{j != i}  m_j (x_j - x_i) / (|x_i - x_j|^2 + eps^2)^(3/2)
```

- **Source:** manual eq 2.1 (force), divided by m_i to get acceleration.
- **Original softening method:** Plummer-softened gravity; see
  Aarseth (1963) for the original idea and Hockney & Eastwood (1981) for
  the modern treatment. Manual eq 2.6 shows the equivalent smoothed-density
  interpretation.

---

## 2. `diagnostics.py`

### `kinetic_energy`

```
KE = sum_i (1/2) m_i |v_i|^2
```

- **Source:** elementary classical mechanics; no citation needed.

### `potential_energy`

```
U = - sum_{i<j} m_i m_j / sqrt(|x_i - x_j|^2 + eps^2)        (G = 1)
```

- **Source:** manual eq 2.5 gives the single-particle potential
  `Phi(x_i) = -sum_{j != i} G m_j / sqrt(r^2 + eps^2)`.
- **Derivation:** total energy is `U = (1/2) sum_i m_i Phi_i`, which sums
  each pair (i, j) twice; rewriting as `sum_{i<j}` removes the factor of
  1/2.

### `angular_momentum`

```
L = sum_i m_i (r_i x v_i)
```

- **Source:** elementary classical mechanics.

---

## 3. `integrator.py`

### `velocity_verlet_step` (KDK form)

```
v_{n+1/2} = v_n + (dt/2) a_n
x_{n+1}   = x_n + dt v_{n+1/2}
a_{n+1}   = a(x_{n+1})
v_{n+1}   = v_{n+1/2} + (dt/2) a_{n+1}
```

- **Source:** manual eq 2.12 (the leapfrog/Verlet method).
  The manual writes it in the staggered-time form
  `x^(n+1) = x^(n) + v^(n+1/2) dt`, `v^(n+1/2) = v^(n-1/2) + F(x^(n))/m dt`.
  The kick-drift-kick (KDK) form used here is mathematically equivalent
  (Hockney & Eastwood 1981) and stores velocities at integer steps,
  which is more convenient for energy diagnostics.
- **Properties:** 2nd-order accurate, time-reversible, symplectic.
  See Hockney & Eastwood (1981) §5; cited at manual eq 2.12.

---

## 4. `initial_conditions.py`

### Disk density

```
rho_d(R, z) = (M_d / (4 pi h^2 z_0)) * exp(-R/h) * sech^2(z/z_0)
```

- **Source:** manual eq 2.38.
- **Original references in the manual:** Freeman (1970) for the exponential
  radial profile; Bahcall & Soneira (1980), Spitzer (1942) for the sech^2
  vertical profile.

### Disk radial CDF (for sampling)

```
P(<R) = 1 - (1 + R/h) exp(-R/h)
```

- **Derivation:** the surface density is the z-integral of eq 2.38,
  `Sigma(R) = (M_d / (2 pi h^2)) exp(-R/h)` (since
  `integral sech^2(z/z_0) dz = 2 z_0`). The radial probability density
  is `dP/dR = 2 pi R Sigma / M_d = (R/h^2) exp(-R/h)`. Integrating
  from 0 to R gives the formula above. Standard exponential-disk
  result; not numbered explicitly in the manual.

### Disk vertical CDF (for sampling)

```
F(z) = (1 + tanh(z/z_0)) / 2
```

- **Derivation:** `integral sech^2(u) du = tanh(u)`. Normalising the
  sech^2 density and integrating from -infty gives the form above.
- **Inverse:** `z = z_0 * atanh(2u - 1)`, used in the sampler.

### Bulge density (Hernquist 1990)

```
rho_b(r) = (M_b a) / (2 pi r (r + a)^3)
```

- **Source:** manual eq 2.67 (cited as Hernquist 1990).
- **Note:** the manual labels this `rho_h` due to a typo; in context
  (Section 2.3.2 stellar bulge), this is the bulge density.

### Bulge cumulative mass

```
M_b(r) = M_b * r^2 / (r + a)^2
```

- **Source:** manual eq 2.69 (Hernquist 1990).

### Bulge inverse CDF

```
r = a * sqrt(u) / (1 - sqrt(u)),     u = M_b(r) / M_b
```

- **Derivation:** invert `u = r^2 / (r + a)^2`. Standard Hernquist
  sampling trick; not in the manual but trivial.

### Halo density

```
rho_h(r) = (M_h alpha / (2 pi^(3/2) r_c))
           * exp(-r^2 / r_c^2) / (r^2 + gamma^2)
```

- **Source:** manual eq 2.57 (phenomenological cored isothermal sphere
  with exponential truncation).
- **Note on heritage:** this is a Hernquist (1993a) modification of the
  isothermal sphere with a Gaussian-cutoff; the manual cites it without
  a primary reference but Hernquist 1993a Appendix A is the canonical
  source.

### Halo normalisation

```
alpha = 1 / (1 - sqrt(pi) * q * exp(q^2) * (1 - erf(q))),    q = gamma / r_c
```

- **Source:** manual eq 2.59.
- **Derivation:** require `integral 4 pi r^2 rho_h(r) dr = M_h` over
  [0, infty). The integral evaluates to
  `(2 alpha M_h / sqrt(pi)) * I(q)` where
  `I(q) = integral_0^infty u^2 e^(-u^2) / (u^2 + q^2) du =
   (sqrt(pi)/2) [1 - sqrt(pi) q e^(q^2) (1 - erf(q))]`.
  Setting the total equal to M_h gives the alpha above.

### Halo CDF (numerical)

- Built by trapezoidal integration of `4 pi r^2 rho_h(r)` on a fine
  radial grid. Inverted by `np.interp`. Standard numerical procedure.

### Isotropic sphere sampling

```
cos(theta) ~ U(-1, 1),     phi ~ U(0, 2*pi)
```

- **Source:** classical result for uniform points on a sphere. Any
  Monte Carlo textbook (e.g. Press et al., *Numerical Recipes*).

---

## 5. `potential.py`

### Spherical-shell averaging of the disk density

```
<rho_d>(r) = (1/(4 pi)) integral_{S^2} rho_d(R(theta, phi), z(theta, phi)) dOmega
           = (1/2)     integral_0^pi   rho_d(r sin theta, r cos theta) sin theta dtheta
```

- **Source:** definition of the spherical-shell average of an
  axisymmetric function. Used by manual page 19 / fig 2.1 caption,
  which explicitly notes the resulting v_c is underestimated by ~ 15 %.
  Heritage: Hernquist (1993a) approximation used by the manual.

### Total enclosed mass

```
M(<r) = integral_0^r 4 pi r'^2 rho_total(r') dr'
```

- **Source:** standard; Newton's first theorem (Binney & Tremaine 1987,
  §2.1).

### Circular speed squared

```
v_c^2(R) = G M(<R) / R                     (spherical approximation)
```

- **Source:** spherical Newton's first theorem; Binney & Tremaine 1987
  eq 2-22. Manual page 19 / fig 2.1 caption notes the ~15 % underestimate
  this approximation introduces for a flattened disk.

### Epicyclic frequency squared (spherical limit)

```
kappa^2(R) = G M(<R) / R^3 + 4 pi G rho_total(R)
```

- **Source manual:** eq 2.49 states the general form
  `kappa^2 = (3/R) dPhi/dR + d^2Phi/dR^2` (Binney & Tremaine 1987 §3.2).
- **Derivation (spherical):** with `dPhi/dR = G M(<R)/R^2`,
  `d^2Phi/dR^2 = G dM/dR / R^2 - 2 G M / R^3 = 4 pi G rho(R) - 2 G M / R^3`.
  Substituting gives
  `kappa^2 = 3 G M/R^3 + 4 pi G rho - 2 G M/R^3 = G M/R^3 + 4 pi G rho`.

### Radial potential derivative

```
dPhi/dr = G M(<r) / r^2
```

- **Source:** spherical Newton's first theorem. Used in the bulge / halo
  Jeans equation below.

### Disk surface density (exact, analytic)

```
Sigma_d(R) = (M_d / (2 pi h^2)) * exp(-R/h)
```

- **Derivation:** `Sigma_d(R) = integral rho_d(R, z) dz`, where the
  z-integral of `sech^2(z/z_0)` is `2 z_0`. The 2 z_0 cancels with the
  1/z_0 in the prefactor of eq 2.38. Not numbered in the manual but a
  standard exponential-disk result (Freeman 1970).

### Isotropic spherical Jeans dispersion (bulge, halo)

```
sigma_r^2(r) = (1 / rho(r)) integral_r^infty rho(r') * (dPhi/dr')(r') * dr'
```

- **Source:** manual eq 2.63 (isotropic limit of eq 2.61).
- **Derivation:** integrate eq 2.61 with beta = 0 (isotropic) and
  apply the boundary condition `rho v_r^2 -> 0` at infinity (the system
  is finite-mass). The `dPhi/dr` here uses the *total* mass distribution,
  per manual text below eq 2.64 ("Phi includes the self-gravity of the
  halo and of all other components").

---

## 6. `velocities.py`

### Disk vertical dispersion

```
sigma_z^2(R) = pi G Sigma_d(R) z_0
```

- **Source manual:** eq 2.40.
- **Derivation:** for an isothermal self-gravitating slab,
  `rho(z) = rho_0 sech^2(z/z_0)`. The vertical Jeans equation
  `sigma_z^2 d(ln rho)/dz = -dPhi/dz` combined with the slab Poisson
  equation `d^2 Phi/dz^2 = 4 pi G rho` yields
  `sigma_z^2 = 2 pi G rho_0 z_0^2 = pi G Sigma_d z_0`, using
  `Sigma_d = 2 rho_0 z_0`.
- **Original literature:** Spitzer (1942), reviewed in Binney &
  Tremaine 1987 §4.3 and Bahcall & Soneira (1980).
- **Caveat:** this assumes disk *self-gravity only*. The bulge and halo
  also contribute to dPhi/dz, so the true vertical-equilibrium
  sigma_z is slightly larger. Manual page 17 (Hernquist 1993a method)
  uses eq 2.40 as an approximation; the disk settles after a "short
  evolution" (page 14).

### Toomre Q stability prescription

```
sigma_R(R_crit) = Q * 3.36 * G * Sigma_d(R_crit) / kappa(R_crit),    Q = 1.5
R_crit = 2.4 * h
```

- **Source manual:** eq 2.47 (critical dispersion), eq 2.48
  (definition of Q), value Q = 1.5 stated on page 17.
- **Original literature:** Toomre (1964); the constant 3.36 is the
  exact value from solving the dispersion relation for axisymmetric
  perturbations in a stellar disk (see also Binney & Tremaine 1987
  §6.2.3).
- **R_crit:** manual page 16 ("we assume a critical radius of
  R_crit = 2.4 h. This corresponds to the solar radius R_sun if the
  model is scaled to the Milky Way.")

### Disk radial dispersion (R-scaling)

```
sigma_R^2(R) = sigma_R^2(R_crit) * Sigma_d(R) / Sigma_d(R_crit)
            propto exp(-R/h)
```

- **Source manual:** eq 2.39 (`v_R^2 propto exp(-R/h)`).
- **Original literature:** observationally established by Freeman (1970),
  Lewis & Freeman (1989), Bottema (1993); see manual page 14.

### Disk azimuthal dispersion (epicyclic)

```
sigma_phi^2(R) = sigma_R^2(R) * kappa^2(R) / (4 Omega^2(R))
```

- **Source manual:** eq 2.53.
- **Original literature:** epicyclic approximation, Binney & Tremaine
  1987 §3.2.3.

### Disk mean azimuthal velocity (asymmetric drift)

```
<v_phi>^2(R) = v_c^2(R) + sigma_R^2(R) * (1 - kappa^2/(4 Omega^2) - 2R/h)
```

- **Source manual:** eq 2.56 in the form
  `<v_phi>^2 - v_c^2 = sigma_R^2 [1 - kappa^2/(4 Omega^2) - R/h +
  d(ln sigma_R^2)/d(ln R) + (R/sigma_R^2) d(v_R v_z)/dz]`.
- **Derivation:** under the model assumptions of this code,
  `sigma_R^2 propto exp(-R/h)`, so
  `d(ln sigma_R^2)/d(ln R) = R * d(ln sigma_R^2)/dR = -R/h`.
  At z = 0, mixed moments `<v_R v_z>` vanish by symmetry, so the last
  term is zero. Substituting reduces the bracket to
  `[1 - kappa^2/(4 Omega^2) - 2 R/h]`.
- **Original literature:** asymmetric drift, Binney & Tremaine 1987
  §4.2.1(d). The "softening of v_R^2 at small radii" (clamping
  `<v_phi>^2 >= 0`) is noted on manual page 17 as Hernquist's
  prescription; we clamp to zero, which is the standard approach.

### Bulge / halo Maxwellian sampling

```
F(v, r) dv = 4 pi v^2 (1/(2 pi sigma^2))^(3/2) exp(-v^2 / (2 sigma^2)) dv
```

- **Source manual:** eq 2.65.
- **Implementation:** equivalently, each Cartesian velocity component
  is sampled i.i.d. from `Normal(0, sigma)`. This produces the 3D
  Maxwellian and is statistically isotropic. Standard probability;
  see Binney & Tremaine 1987 §4.4.

---

## 7. `stability_test.py`

### Disk scale-length fit (exponential)

```
h = -1 / slope of ( ln Sigma_d(R) vs R )
```

- **Derivation:** taking the log of `Sigma_d(R) = Sigma_0 exp(-R/h)`
  gives `ln Sigma_d = ln Sigma_0 - R/h`, a straight line of slope
  -1/h. Standard linear regression.

### Disk scale-height fit (sech^2)

```
z_0 = sqrt(12 * <z^2> / pi^2)
```

- **Derivation:** for the normalised density `rho(z) = (1/(2 z_0)) sech^2(z/z_0)`,
  using the standard integral
  `integral_{-infty}^{infty} u^2 sech^2(u) du = pi^2 / 6`
  (Gradshteyn & Ryzhik 3.527; equivalent to a residue calculation
  or differentiation of `tanh`'s generating function),
  we get
  ```
  <z^2> = (z_0^2 / 2) * (pi^2 / 6) = pi^2 z_0^2 / 12
  ```
  Inverting yields the formula above. **Earlier code had the
  bare-integral value `pi^2 / 6` in place of <z^2>, which omits the
  factor of 1/2 from the density normalisation — fixed 2026-05-13.**

---

## Manual citations -- abbreviated bibliography pointers

These are the manual's bibliography entries (manual pages 25-29) that
underpin formulas used here:

- Aarseth, S. (1963) — softening method (general)
- Bahcall, J. N. & Soneira, R. M. (1980) — sech^2 disk vertical profile
- Barnes, J. E. (1986) — tree algorithm (not used; direct N^2 chosen)
- Binney, J. & Tremaine, S. (1987) — *Galactic Dynamics*, the foundational
  textbook; nearly all kinematic identities trace here
- Bottema, R. (1993) — observed disk velocity-dispersion radial scaling
- Freeman, K. C. (1970) — exponential disk surface density
- Hernquist, L. (1990) — bulge analytic profile (manual eq 2.67-2.69)
- Hernquist, L. (1992, 1993b) — system of units used on manual page 18
- Hernquist, L. (1993a) — N-body realisations of compound galaxies; the
  IC construction recipe (manual page 14 onward)
- Hockney, R. W. & Eastwood, J. W. (1981) — *Computer Simulations Using
  Particles*; leapfrog discretisation (manual eq 2.12)
- Lewis, J. R. & Freeman, K. C. (1989) — disk radial dispersion scaling
- Rybicki, G. B. (1971) — two-body relaxation timescale with softening
- Spitzer, L. (1942) — isothermal slab, sech^2 profile
- Toomre, A. (1964) — Q stability criterion (manual eq 2.47-2.48)

For background on the broader code style (no formula derivations, just
conventions), see
`../../../.claude/projects/-Users-zoutongshen-Library-CloudStorage-Dropbox-Academics-02-Courses-26Spring-CMPH/memory/feedback_code_guidelines.md`
(local-only, Schaller 2026 course guidelines).
