# Two-Dimensional FFT Analysis of Air-Gap Flux Density

## Scope

This note describes a general post-processing method for transient electric-machine
simulations that provide radial and circumferential air-gap flux density. The method
uses a two-dimensional discrete Fourier transform (2-D DFT) to identify coupled time
and spatial harmonics in the magnetic field and in the Maxwell stress.

The same framework supports analysis of torque ripple, radial electromagnetic force,
and unbalanced magnetic pull (UMP). It is independent of a particular motor geometry
or solver, provided that the sampling and coordinate conventions are defined.

## Sampled field data

Let the sampled field be `Bq(tj, theta_k)`, where `q` is either `r` or `theta`:

- `Br` is the radial flux density in tesla.
- `Btheta` is the circumferential flux density in tesla.
- `tj = t0 + j Delta_t`, for `j = 0, ..., Nt - 1`.
- `theta_k = 2 pi k / Ntheta`, for `k = 0, ..., Ntheta - 1`.

The angular samples must represent one mechanical circumference after any sector
expansion. The first and last samples of a periodic time interval must not be included
twice in the FFT input. Uniform sampling in both dimensions is assumed.

## Two-dimensional DFT

For either field component, the normalized coefficient at time order `n` and spatial
order `m` is

$$
\hat{B}_q[n,m] = \frac{1}{N_t N_\theta}
\sum_{j=0}^{N_t-1}\sum_{k=0}^{N_\theta-1}
B_q(t_j,\theta_k)
\exp\left[-\mathrm{i}2\pi\left(\frac{nj}{N_t}+\frac{mk}{N_\theta}\right)\right].
$$

The time-order axis is the DFT-bin order. Its corresponding frequency is
`fn = n / (Nt Delta_t)`. The spatial order `m` is the number of sinusoidal cycles
around the mechanical circumference. Positive and negative orders are a paired complex
representation of the same real field; their phase is meaningful when tracing modal
interactions.

## Sector expansion and symmetry

A model that covers only one sector must be expanded to a full mechanical circumference
before the spatial FFT. Let the sector angle be `alpha = 2 pi / P`, and let `ell` be
the repeated-sector index. The field expansion is

$$
B_q(t, \theta + \ell\alpha) = s^\ell B_q(t, \theta),
$$

where `s = +1` for periodic symmetry and `s = -1` for anti-periodic symmetry.

The anti-periodic sign reversal must be applied to both `Br` and `Btheta` in each
alternating sector. Omitting it changes the permitted spatial orders and can produce
incorrect force and torque harmonics.

## Maxwell stress

Using the permeability of free space `mu0`, the radial and circumferential Maxwell
stress components evaluated in the air gap are

$$
F_r(t,\theta) = \frac{B_r^2(t,\theta)-B_\theta^2(t,\theta)}{2\mu_0},
\qquad
F_\theta(t,\theta) = \frac{B_r(t,\theta)B_\theta(t,\theta)}{\mu_0}.
$$

A 2-D DFT is then applied to `Fr` and `Ftheta` in the same way as to the flux-density
components. Because the stress is quadratic in flux density, interactions of two field
modes create sum and difference orders:

$$
(n,m) = (n_1 \mathbin{\pm} n_2,\; m_1 \mathbin{\pm} m_2).
$$

This relationship is useful when linking a magnetic sideband to a force or torque
component.

## Integrated force and torque

For air-gap radius `R` and active stack length `L`, the Cartesian radial force is
obtained from the radial stress:

$$
F_x(t) = R L \int_0^{2\pi} F_r(t,\theta)\cos\theta\,d\theta,
\qquad
F_y(t) = R L \int_0^{2\pi} F_r(t,\theta)\sin\theta\,d\theta.
$$

The UMP magnitude is `sqrt(Fx^2 + Fy^2)`. The spatial first-order component of `Fr`
therefore determines the mean force direction and magnitude. A static eccentricity,
for example, commonly appears as an increase in the `m = 1` radial-stress component.

The air-gap torque is calculated from the circumferential stress:

$$
T(t) = -R^2 L \int_0^{2\pi} F_\theta(t,\theta)\,d\theta.
$$

The sign follows the selected rotor-positive convention and should be checked against
the electromagnetic solver output. Only the spatial zeroth-order component of
`Ftheta` contributes to the integrated torque; its nonzero time orders describe torque
ripple.

## Reading mode maps

A sparse 2-D mode map plots time order on the horizontal axis, spatial order on the
vertical axis, and coefficient magnitude by marker size or color. It is normally most
readable to display positive spatial orders only, while retaining the complex spectra
for quantitative calculations.

Useful selection rules are:

- Inspect `Fr(n = 0, m = 1)` and its nearby time orders for UMP-related behavior.
- Inspect `Ftheta(n, m = 0)` for torque ripple orders.
- Compare `Br` and `Btheta` modes with the stress modes using the sum-and-difference
  relationship to identify likely sources of a force component.
- Compare reconstructed torque with an independently integrated solver torque before
  interpreting small harmonics.

## Practical checks

1. Confirm units: flux density in tesla, radius and stack length in metres, and stress
   in pascals.
2. Expand sector data with the declared periodic or anti-periodic symmetry.
3. Remove a duplicated periodic endpoint before applying the time FFT.
4. Use the same angular grid and full-circumference normalization for field, stress,
   force, and torque calculations.
5. Retain complex FFT coefficients for phase-sensitive comparisons; use magnitudes only
   for display and ranking.

The repository user guide describes the supported JSON and CSV inputs and the generated
analysis files: [Air-gap force 2D FFT guide](README_airgap_force_fft.md).