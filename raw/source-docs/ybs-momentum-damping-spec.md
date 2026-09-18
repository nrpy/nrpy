# YBS-MOM: Hyperbolic-Character Momentum-Constraint Cleaning for BSSN

## Status and naming

This document specifies a proposed numerical adjustment to vacuum BSSN. It
replaces the parabolic momentum-gradient feedback of Yo, Lin, and Cao with an
auxiliary-field relaxation system whose added principal block has finite
characteristic speeds. The construction is a new synthesis, not an equation
claimed to appear in the cited papers.

The filename and label **YBS-MOM** do not identify this proposal with the
historical Yo--Baumgarte--Shapiro (YBS) adjustment. Three terms must remain
distinct:

1. Yo, Lin, and Cao's momentum adjustment adds
   \(h f(\alpha)\mathcal M_{\langle i,j\rangle}\) to the conformal
   trace-free extrinsic-curvature equation. Since \(\mathcal M_i\) already
   contains a derivative of \(\bar A_{ij}\), this is parabolic at principal
   order.
2. The historical YBS adjustment instead adds the Gamma constraint to the
   conformal-connection evolution equation. It is not the momentum-gradient
   term and is not the cleaner proposed here.
3. Etienne's coarse-grid-adjusted Hamiltonian-constraint damping (CAHD)
   remains parabolic. Only its level-local AMR/CFL coefficient scaling is
   transferred below; CAHD is not a prior hyperbolization.

## Conventions

Use geometric units \(G=c=1\) and the ADM sign convention

\[
K_{ij}=-\frac{1}{2\alpha}
       (\partial_t-\mathcal L_\beta)\gamma_{ij}.
\]

Write

\[
\gamma_{ij}=e^{4\phi}\bar\gamma_{ij},\qquad
K_{ij}=e^{4\phi}
\left(\bar A_{ij}+\frac13\bar\gamma_{ij}K\right),
\]

with

\[
\det\bar\gamma=\det\hat\gamma,\qquad
\bar\gamma^{ij}\bar A_{ij}=0.
\]

Here \(\hat\gamma_{ij}\) is the reference metric; in Cartesian BSSN,
\(\det\hat\gamma=1\). Let \(\bar D_i\) be the covariant derivative compatible
with \(\bar\gamma_{ij}\), and raise barred indices with
\(\bar\gamma^{ij}\). The lower-index vacuum momentum constraint is

\[
\boxed{
\mathcal M_i
=\bar D_j\bar A^j{}_i
 +6\bar A^j{}_i\,\partial_j\phi
 -\frac23\partial_iK
=0 .}
\tag{1}
\]

An implementation using \(W=e^{-2\phi}\) may replace
\(6\partial_j\phi\) by \(-3\partial_jW/W\). A matter implementation must use
its own consistently normalized, matter-complete momentum residual. An
upper-index, physically raised, or reference-metric-rescaled diagnostic cannot
be substituted into the equations below without converting it to Eq. (1) or
rederiving the coefficients.

For a covector \(v_i\), define

\[
\mathcal S_{ij}[v]
\equiv\bar D_{\langle i}v_{j\rangle}
=\frac12(\bar D_iv_j+\bar D_jv_i)
 -\frac13\bar\gamma_{ij}\bar D_kv^k .
\tag{2}
\]

Angle brackets therefore denote the symmetric trace-free (STF) projection in
three spatial dimensions. Define the advective tensor derivative

\[
\mathcal D_t\equiv\partial_t-\mathcal L_\beta .
\tag{3}
\]

The new field below is a weight-zero covector, so

\[
(\mathcal L_\beta Q)_i
=\beta^k\partial_kQ_i+Q_k\partial_i\beta^k.
\tag{4}
\]

Codes whose stored \(\bar A_{ij}\) is a tensor density or a rescaled
reference-metric variable must convert the added tensor consistently with the
host BSSN right-hand side.

## Proposed system

Introduce three evolved components of an auxiliary covector \(Q_i\). On an
active AMR level \(n\), evolve

\[
\boxed{
\begin{aligned}
\mathcal D_t\bar A_{ij}
 &= (\mathcal D_t\bar A_{ij})_{\rm BSSN}
   +\ell_{M,n}\,\mathcal S_{ij}[Q],\\[1mm]
\tau_{M,n}\mathcal D_tQ_i
 &=\mathcal M_i-Q_i .
\end{aligned}}
\tag{5}
\]

Require \(\ell_{M,n}>0\) and \(\tau_{M,n}>0\) whenever the cleaner is
enabled. The plus sign is tied to the momentum-constraint convention in
Eq. (1). The first addition is explicitly STF, although the host code's usual
finite-precision trace projection remains necessary.

The baseline contains no lapse window: \(f(\alpha)=1\). A window that reaches
zero can make a one-sided first-order coupling non-diagonalizable. Any future
window must multiply the two coupling directions consistently and requires a
new variable-coefficient principal and CFL analysis.

If the cleaner is disabled, the code must omit Eq. (5) and \(Q_i\), or freeze
the entire auxiliary block behind an explicit feature switch. It must not
evolve a one-way \(Q_i\) equation at zero strength.

### Constraint surface and initialization

The augmented exact constraint surface consists of all base BSSN constraints,
including

\[
\mathcal H=0,\quad \mathcal M_i=0,\quad \mathcal C^i=0,\quad
\det\bar\gamma-\det\hat\gamma=0,\quad
\bar\gamma^{ij}\bar A_{ij}=0,
\]

together with \(Q_i=0\). Both new right-hand sides then vanish. Thus, if the
base continuum equations preserve their full constraint surface, every exact
BSSN solution embedded with \(Q_i=0\) remains an exact solution of the
augmented system.

Initialize \(Q_i=0\). This avoids turning differentiation noise in numerical
initial data into an immediate curvature correction. Preserve \(Q_i\) in
checkpoint/restart data and prolong it when creating new AMR points.

## Frozen-coefficient principal analysis

This section analyzes only the added block. Freeze the metric, shift, and
coefficients at a point, use a local conformal orthonormal frame, discard base
BSSN and lower-order couplings, and write

\[
\mathcal M_i\simeq\partial^j\bar A_{ij}.
\]

The cleaning contribution obeys

\[
\mathcal D_t\mathcal M_i
=\ell_M\left(
 \frac12\Delta Q_i
 +\frac16\partial_i\partial^jQ_j
\right).
\tag{6}
\]

Using the second equation in (5) gives the momentum-constraint telegraph
equation

\[
\boxed{
\tau_M\mathcal D_t^2\mathcal M_i
 +\mathcal D_t\mathcal M_i
=\ell_M\left(
 \frac12\Delta\mathcal M_i
 +\frac16\partial_i\partial^j\mathcal M_j
\right).}
\tag{7}
\]

For a Fourier direction \(s_i=k_i/|k|\), decompose
\(\mathcal M_i=\mathcal M_i^T+s_i\mathcal M_L\), with
\(s^i\mathcal M_i^T=0\). The two transverse and one longitudinal families
have speeds relative to \(\mathcal D_t\)

\[
\boxed{
c_T=\sqrt{\frac{\ell_M}{2\tau_M}},\qquad
c_L=\sqrt{\frac{2\ell_M}{3\tau_M}}.}
\tag{8}
\]

For modes proportional to \(e^{\sigma t+i k_jx^j}\), define the
shift-relative growth rate \(\widehat\sigma=\sigma-i\beta^jk_j\). The
dispersion relations are

\[
\tau_M\widehat\sigma_T^2+\widehat\sigma_T
 +\frac{\ell_M}{2}|k|^2=0,
\qquad
\tau_M\widehat\sigma_L^2+\widehat\sigma_L
 +\frac{2\ell_M}{3}|k|^2=0.
\tag{9}
\]

Every nonzero-wavenumber root has negative real part for positive
\(\ell_M,\tau_M\); an underdamped mode has envelope
\(e^{-t/(2\tau_M)}\). A spatially uniform momentum-constraint violation is
not removed because the feedback acts through a gradient.

The isolated first-order principal symbol also has a complete directional
eigenbasis. With transverse basis covectors \(e_A^i\), define
\(\bar A_{sA}=s^ie_A^j\bar A_{ij}\),
\(\bar A_{ss}=s^is^j\bar A_{ij}\), \(Q_A=e_A^iQ_i\), and
\(Q_s=s^iQ_i\). Characteristic combinations may be chosen as

\[
w_A^\pm=\bar A_{sA}
 \pm\sqrt{\frac{\ell_M\tau_M}{2}}\,Q_A,
\qquad
w_L^\pm=\bar A_{ss}
 \pm\sqrt{\frac{2\ell_M\tau_M}{3}}\,Q_s .
\tag{10}
\]

They propagate at \(\pm c_T\) and \(\pm c_L\), respectively, relative to
\(\mathcal D_t\); the two remaining STF tensor combinations have zero speed
within the added block. The shift translates these into the corresponding
coordinate speeds. These are frozen-background cleaner eigenfields, not a
boundary decomposition of the fully coupled BSSN system.

For constant coefficients and vanishing boundary flux, the cleaning-only
energy

\[
E=\frac12\int\left(
 \bar A_{ij}\bar A^{ij}+\ell_M\tau_M Q_iQ^i
\right)d^3x
\]

satisfies

\[
\frac{dE}{dt}=-\ell_M\int Q_iQ^i\,d^3x\le0.
\tag{11}
\]

This energy check fixes the matched signs in Eq. (5). It applies only to the
isolated, frozen, boundary-free cleaning block.

### Quasi-steady limit

When \(\tau_M\mathcal D_t\mathcal M_i\) is small relative to
\(\mathcal M_i\), iteration of the auxiliary equation gives

\[
Q_i=\mathcal M_i-\tau_M\mathcal D_t\mathcal M_i
 +O(\tau_M^2\mathcal D_t^2\mathcal M_i).
\tag{12}
\]

Hence the leading curvature correction is

\[
\mathcal D_t\bar A_{ij}
=(\mathcal D_t\bar A_{ij})_{\rm BSSN}
 +\ell_M\bar D_{\langle i}\mathcal M_{j\rangle},
\tag{13}
\]

with the first omitted term
\(-\ell_M\tau_M\bar D_{\langle i}
\mathcal D_t\mathcal M_{j\rangle}\). Equation (13) is the covariant
principal form of the Yo--Lin--Cao adjustment, with their \(h f(\alpha)\)
replaced by \(\ell_M\). Its transverse and longitudinal diffusion rates are
\(-\ell_M|k|^2/2\) and \(-2\ell_M|k|^2/3\). The positive sign in Eq. (5)
therefore damps rather than antidamps this convention for \(\mathcal M_i\).

## CAHD-inspired AMR scaling

Let \(\Delta s_n\) be the level-local cell scale used by the evolution's
discrete CFL estimate, \(\Delta t_n\) the timestep actually used on level
\(n\), and \(\mathrm{CFL}_0\) the reference Courant factor on the
finest-timestep level. For Cartesian grids,
\(\Delta s_n=\min_i\Delta x_n^i\). For mapped grids, derive the scale from
the same local derivative-symbol radius used by the hyperbolic timestep
estimate.

Transfer Etienne's CAHD scaling directly to the quasi-steady momentum
coefficient:

\[
\boxed{
\ell_{M,n}
=C_M\,\Delta t_n\,\mathrm{CFL}_0
 \left(\frac{\Delta s_n}{\Delta t_n}\right)
=C_M\,\mathrm{CFL}_0\,\Delta s_n,}
\tag{14}
\]

where \(C_M>0\) is a new, dimensionless tunable strength. The unsimplified
form records the CAHD construction; code should evaluate the simplified form.
It is independent of whether the level subcycles, increases on coarser levels,
and vanishes linearly under uniform refinement.

Use the simplest one-parameter relaxation baseline

\[
\boxed{\tau_{M,n}=\Delta s_n.}
\tag{15}
\]

Then the cleaning speeds are level independent:

\[
c_{T,n}=\sqrt{\frac{C_M\mathrm{CFL}_0}{2}},\qquad
c_{L,n}=\sqrt{\frac{2C_M\mathrm{CFL}_0}{3}}.
\tag{16}
\]

Equation (15) fixes a reference coordinate signal speed of one; it is a
resolution-scale design choice, not an empirical calibration. No value of
\(C_M\) is inherited from the different CAHD operator. Determine an allowed
and useful range for the chosen gauge, discretization, AMR hierarchy, and
physical problem.

## Timestep requirements

Hyperbolization removes the independent parabolic restriction
\(\Delta t=O(\Delta s^2/\ell_M)\), but two ordinary method-of-lines
restrictions remain.

First, include the cleaner in the full hyperbolic spectral radius. In a
direction \(s_i\), the cleaner contribution is bounded by
\(|\beta^s|+c_L\), so schematically

\[
\Delta t_n\le C_{\rm hyp}
\frac{\Delta s_n}
{\max\!\left[\rho_{\rm BSSN}(s),\,|\beta^s|+c_{L,n}\right]}.
\tag{17}
\]

The code's multidimensional full-system estimate takes precedence over this
one-dimensional schematic bound.

Second, explicit integration of \(-Q_i/\tau_{M,n}\) requires

\[
\boxed{\frac{\Delta t_n}{\tau_{M,n}}\le R_{\rm neg},}
\tag{18}
\]

where \(R_{\rm neg}\) is the negative-real-axis stability radius of the
actual time integrator, including its chosen safety factor. Under Eq. (15),
this becomes a level-independent bound on \(\Delta t_n/\Delta s_n\) when
levels subcycle at a fixed Courant ratio. An exact split source update or an
IMEX scheme can relax Eq. (18), but not Eq. (17) or the need for stage-consistent
coupling.

## Implementation and validation requirements

At each method-of-lines stage:

1. Fill BSSN and \(Q_i\) ghost zones at the same stage time.
2. Compute the lower-index \(\mathcal M_i\) from the current stage using
   production-compatible derivative operators and the normalization in
   Eq. (1).
3. Compute \(\bar D_iQ_j\), apply Eq. (2), and add
   \(+\ell_{M,n}\mathcal S_{ij}[Q]\) to the existing \(\bar A_{ij}\) right-hand
   side. Keep \(\ell_{M,n}\) outside the derivative.
4. Add \((\mathcal M_i-Q_i)/\tau_{M,n}\), together with the covector Lie
   transport in Eq. (4), to the coordinate-time \(Q_i\) right-hand side.
5. Evolve \(Q_i\) through prolongation, restriction, synchronization,
   checkpoint/restart, parity, and interpatch transformations as a covector.

Do not reset \(Q_i=\mathcal M_i\) at each stage; that restores the parabolic
system. Do not retain the old direct momentum-gradient term concurrently
unless deliberately defining and analyzing a mixed system.

Level jumps in \(\ell_{M,n}\) and \(\tau_{M,n}\) change cleaning impedance and
may reflect cleaning waves. Measure these reflections before considering a
coefficient taper. With subcycling, time-prolong \(Q_i\) at the same order as
the BSSN state. At physical boundaries, derive conditions from the full
coupled characteristic system; the isolated fields in Eq. (10) are useful
diagnostics but are not by themselves a constraint-preserving boundary
condition.

Minimum validation includes:

- zero added right-hand sides on exact constraint-satisfying data with
  \(Q_i=0\);
- STF contraction of the added \(\bar A_{ij}\) right-hand side;
- transverse and longitudinal plane-wave speeds and damping from
  Eqs. (8)--(9);
- the quasi-steady relation \(Q_i\simeq\mathcal M_i\) and diffusion rates in
  Eq. (13);
- uniform-grid convergence and a pulse crossing a 2:1 refinement interface;
- simultaneous checks of momentum, Hamiltonian, Gamma, determinant, and trace
  constraints; and
- convergence of physical observables, not merely smaller constraint norms.

## Dimensions and claim boundary

In geometric units,

\[
[\bar A_{ij}]=L^{-1},\quad
[\mathcal M_i]=[Q_i]=L^{-2},\quad
[\ell_M]=[\tau_M]=L,
\]

while \(C_M\) and \(\mathrm{CFL}_0\) are dimensionless. Thus
\([\ell_M\bar DQ]=L^{-2}=[\mathcal D_t\bar A]\) and
\([(\mathcal M-Q)/\tau_M]=L^{-3}=[\mathcal D_tQ]\).

Equations (5)--(10) establish finite speeds and a complete eigenbasis for the
isolated frozen added block when enabled. They do **not** prove strong
hyperbolicity or well-posedness of the full gauge-fixed BSSN system, nor do
they prove nonlinear stability, monotone finite-domain constraint decay,
AMR-interface stability, or improved waveform accuracy. Those claims require
the full coupled principal symbol, compatible boundary analysis, and numerical
validation.

## References

- H.-J. Yo, C.-Y. Lin, and Z. Cao, “Modifications for numerical stability of
  black hole evolution,” [arXiv:1205.5111](https://arxiv.org/abs/1205.5111),
  [Phys. Rev. D 86, 064027 (2012)](https://doi.org/10.1103/PhysRevD.86.064027).
  Equations (52)--(56) motivate the momentum-gradient adjustment; Eq. (47)
  identifies the distinct YBS Gamma-constraint adjustment.
- H.-J. Yo, T. W. Baumgarte, and S. L. Shapiro, “Improved numerical stability
  of stationary black hole evolution calculations,”
  [arXiv:gr-qc/0209066](https://arxiv.org/abs/gr-qc/0209066),
  [Phys. Rev. D 66, 084026 (2002)](https://doi.org/10.1103/PhysRevD.66.084026).
- Z. B. Etienne, “Improved Moving-Puncture Techniques for Compact Binary
  Simulations,” [arXiv:2404.01137](https://arxiv.org/abs/2404.01137),
  [Phys. Rev. D 110, 064045 (2024)](https://doi.org/10.1103/PhysRevD.110.064045).
  Equation (26) supplies the level-local AMR/CFL scaling; its CAHD principal
  equation remains parabolic.
