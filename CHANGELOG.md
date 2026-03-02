# DEBUG_L4PTSAL Override Changelog

All changes are to override files in `DEBUG_L4PTSAL/`, never to `src/`.

---

## 2026-03-02

### Depth-Integrated Momentum Budget Output

#### cppdefs.opt
- Added `#define DIAGNOSTICS` -- enables 7-term momentum diagnostic budget.

#### diagnostics.opt (NEW override)
- Set `diag_uv = .true.`, `output_period = 900`, `nrpf = 24`.

#### diagnostics.F (NEW override)
- Added depth-integrated diagnostic arrays `Udiag_int(nx,ny,7)`,
  `Vdiag_int(nx,ny,7)`. Computed as sum over k of Udiag(i,j,k,n) in new
  subroutine `calc_diag_int`. Averaged version in `calc_diag_int_avg`.
  Written to diagnostics NetCDF as `u_pgr_int`, `u_cor_int`, `u_adv_int`,
  `u_dis_int`, `u_hmx_int`, `u_vmx_int`, `u_cpl_int` (and `v_*` equivalents).
  Units: m^2/s^2 (depth-integrated du/dt contribution per term).

#### wec_frc.F
- Added `hvf_int_xi`, `hvf_int_eta` -- depth-integrated horizontal vortex
  force. Accumulated in `compute_wec_rhs_uv_terms` from the Stokes-drift x
  velocity shear terms (cff1). Reset at start of each call (keeps only
  corrector values). Written to WEC output.
- Added `r_D` output (bottom drag coefficient from coupling module).

#### wec_frc.opt
- Added `wrt_hvf_int = .true.`, `wrt_rD = .true.`.

### Diagnostic Term Assignment Fixes (double-counting audit)

#### prsgrd.F
- **Moved DIAGNOSTICS `ipgr` capture** from BEFORE the `#ifdef WEC` block to AFTER
  it. Previously, `ipgr = ru * dxdyi` was captured at the end of the main k-loop
  (after hydrostatic PG only), BEFORE the WEC pressure additions (kvf, sup, calP,
  Kapsrf). Now ipgr correctly includes the total pressure gradient with all WEC
  corrections. Without WEC, behavior is unchanged.

#### compute_horiz_rhs_uv_terms.h (NEW override)
- **Added Stokes-Coriolis to `icori` diagnostic**. Previously, icori captured only
  `UFx` (Coriolis on Eulerian v) for u-momentum and `VFe` (Coriolis on Eulerian u)
  for v-momentum. With WEC, `UFe` (= f * v_stokes) and `VFx` (= f * u_stokes) are
  also added to ru/rv but were missing from icori -- they fell into the `idiss`
  residual instead. Now icori includes both Eulerian and Stokes-Coriolis, guarded
  by `#ifdef WEC`. Without WEC, behavior is unchanged.

### WEC Isolation from Diagnostic Residual (iwec)

#### diagnostics.F
- Increased `nmd` from 7 to 8. Added `iwec = 8` -- captures all WEC explicit
  RHS contributions (horizontal vortex force J + body forces + streaming).
- Added `Udiag = 0; Vdiag = 0` after allocation (ensures iwec=0 when WEC off).
- Added 3D NetCDF variables `u_wec`, `v_wec` and depth-integrated `u_wec_int`,
  `v_wec_int` to output files and writing routines.

#### step3d_uv1.F
- Added `ru_save`, `rv_save` local arrays (`#if defined WEC && defined DIAGNOSTICS`).
- Before `compute_wec_rhs_uv_terms`: save `ru`, call subroutine, then compute
  `ru_save = ru - ru_save` (WEC subroutine contribution only, excluding vertical
  advection which is added later inside the j-loop).
- Inside j-loop after frc2dx: capture `Udiag(:,:,:,iwec)` from `ru_save*dxdyi`
  plus `dm*dn*frc2dx*dxdyi` at k=1 (streaming, only when `!BODY_FRICTION`).
- Modified `idiss = ru*dxdyi - ipgr - icori - iadv - iwec` (was `- iadv` only).
  Now `idiss` is approximately pure numerical dissipation from advection.

**Budget mapping after all fixes:**
```
ipgr   = total PGF including WEC corrections (kvf, sup, calP, Kapsrf)
icori  = Coriolis + Stokes-Coriolis (f*v_euler + f*v_stokes)
iadv   = horizontal advection (centered 4th-order fluxes)
iwec   = vortex force J + streaming + body forces (all WEC explicit RHS)
idiss  ≈ numerical dissipation from advection scheme
ihmix  = horizontal viscous mixing
ivmix  = implicit solver: wind + breaking + bottom drag + vertical diffusion
icoup  = 2D/3D coupling correction
```

**ivmix decomposition (offline from existing outputs):**
- Σ ivmix = sustr + τ_bot  (vertical diffusion is conservative, sums to 0)
- τ_wind = sustr - brk2dx   (SURFACE_BREAK contaminates sustr)
- τ_bot  = -r_D * u_bot     (from r_D and bottom velocity outputs)
- Σ brk2dx already output separately as 2D WEC field

---

## 2026-02-28

### compute_shelter_fac.h
- **Fixed geometric artifacts** in nearshore velocity field caused by
  SHELTER_BRK. Two changes:
  1. **Cosine taper** replaces linear ramp in Step 3 conversion. The linear
     ramp `(1-LI)/(1-Lf0)` had a derivative discontinuity at the Lf0 contour,
     creating sharp gradient lines in shelter_fac that imprinted onto brk.
     Cosine taper `0.5*(1+cos(pi*(LI-Lf0)/(1-Lf0)))` has zero derivative at
     both endpoints (LI=Lf0 and LI=1), eliminating the kink entirely.
  2. **Post-smoothing** (Step 4, N_post=4 passes) of shelter_fac after
     conversion. Uses the same isotropic Laplacian stencil but WITH masking
     (`*umask`, `*vmask` on flux terms). Removes staircase artifacts from the
     grid-aligned discrete coastline mask. Land cells stay at 0.

### wec_frc.F
- Added `pi_sh` to local variable declarations for SHELTER_BRK
  (needed by cosine taper in compute_shelter_fac.h).

---

## 2026-02-26

### cppdefs.opt
- Added `#define SHELTER_BRK` under `#ifdef BREAKING` — coastline shelter
  factor that scales down `brk(i,j)` in concave coastal regions
  (embayments) where WW3's coarser coastline could not resolve wave
  sheltering. Applied after `SMOOTH_BRK`.

### wec_frc.F
- Added `shelter_fac(GLOBAL_2D_ARRAY)` module-level array (SHELTER_BRK
  guarded). Precomputed once during `init_arrays_wec_tile` from `rmask`.
- Widened scratch array guard (`FX_brk`, `FE_brk`, `FE1_brk`) from
  `SMOOTH_BRK` to `SMOOTH_BRK || SHELTER_BRK` (shared for both).
- Added shelter_fac allocation in `init_arrays_wec_tile`.
- Added local variables `iter_s`, `cff_sh`, `cff1_sh`, `Lf0_sh`,
  `cff_sum_sh` in `init_arrays_wec_tile` for shelter computation.
- Added `#include "compute_shelter_fac.h"` at end of
  `init_arrays_wec_tile` (BREAKING && SHELTER_BRK guarded).
- Added `brk = brk * shelter_fac` in `wec_forces_tile` after
  SMOOTH_BRK block (BREAKING && SHELTER_BRK && !BRK0 guarded).

### compute_shelter_fac.h (NEW)
- One-time computation of shelter factor from land mask.
  Algorithm: (1) initialize `shelter_fac = 1 - rmask` (land influence
  field), (2) apply 8 passes of isotropic Laplacian smoother (same
  stencil as smooth_brk.h but WITHOUT masking — land signal must bleed
  into ocean), (3) convert to scaling factor with threshold Lf0=0.5:
  `shelter_fac = max(0, (1-land_influence)/(1-Lf0))` where
  land_influence > Lf0, else shelter_fac=1.
  Physical motivation: WW3 ran at coarser resolution with smoother
  coastline. Fine-scale embayments in ROMS would shelter waves but
  WW3 cannot see them. shelter_fac corrects this mismatch.
  Parameters: N_shelter=8 (σ≈52m at Δx=30m), Lf0=0.5 (straight-coast
  threshold). In embayments: land bleeds from 3 sides → high
  land_influence → shelter_fac < 1. On straight coast: land from
  1 side → land_influence ≈ Lf0 → shelter_fac ≈ 1.

### extract_data.F (NEW override, copied from src/)
- **Fixed `find_local_points` for non-contiguous (scattered) points.**
  The original algorithm assumed points within each MPI subdomain form a
  contiguous block in the input array (true for boundary extraction, false
  for scattered moorings). It scanned forward from the first local point
  and stopped at the first non-local point, silently dropping any local
  points after the gap — or worse, including non-local points whose
  `ip`/`jp` indices fall outside the subdomain's array bounds, causing
  segfaults and MPI deadlocks.
  Replaced with a two-pass algorithm: (1) count all local points,
  (2) collect them with a `global_idx` mapping. Works for both contiguous
  (boundary) and scattered (mooring) point sets.
- **Added `global_idx` array** to `extract_object` type: maps local point
  index to global input array index, used for angle extraction.
- **Fixed angle extraction** in `read_extraction_objects`: uses `global_idx`
  instead of the contiguous `i0:i1` range to extract rotation angles from
  the input array. Old code: `obj%ang = object(i0:i1,3)` (wrong for
  scattered points). New code: loops over `global_idx`.
- **Guarded divide-by-zero in `compute_coef`** (line 418): when all 4
  stencil cells are masked (e.g. u/v stencil for a point near coastline
  like Port Los Angeles), `sum(coef)=0` caused NaN → crash. Now sets
  `coef=0` so interpolation returns 0.0 (effectively skipping that
  variable at that point). Scalar fields (zeta) with at least one ocean
  corner cell still interpolate correctly via the rho-mask stencil.

---

## 2026-02-25

### cppdefs.opt
- Added `#define SMOOTH_BRK` under `#ifdef BREAKING` — Laplacian smoothing
  of brk(i,j) to spread breaking momentum over minimum surfzone width.

### wec_frc.F
- Added `FX_brk`, `FE_brk`, `FE1_brk` scratch arrays (module-level,
  `BREAKING && SMOOTH_BRK` guarded).
- Added allocations in `init_wec_arrays`.
- Added local variables `iter_brk`, `cff_s`, `cff1_s` and hardcoded
  parameter `N_brk_smooth=4` in `wec_forces_tile`.
- Inserted `#include "smooth_brk.h"` after BRK0 block, before FRC0 block,
  guarded by `BREAKING && SMOOTH_BRK && !BRK0`.

### smooth_brk.h (NEW)
- Iterative isotropic Laplacian smoothing of `brk(i,j)` at rho-points.
  Adapted from `lmd_kpp_smooth_hbl.h`. `N_brk_smooth` passes with MPI
  exchange between each pass. Mask-aware (Neumann BC at coastline).
  Conservative. Stencil: isotropic Laplacian (cff=1/12, cff1=3/16).
  Benefits both SURFACE_BREAK (stress) and BODY_BREAK (body force) paths
  since both read brk(i,j). Only brk is smoothed, not ebrk (Akb mixing
  unaffected).

---

## 2026-02-19

### cppdefs.opt
- `WAVE_DIFF_FROM_LM`: changed from `#define` to `#undef`.
  Reason: Akw is hardcoded to 0.0 (wec_frc.F:2080), so the entire
  `llmm → keffkw → Akw` pipeline is dead code. No need to read the
  extra `wwlm` field from WW3 each forcing time step.

### wec_frc.F
- **Re-instated AKB0 guard** (line 2056): changed
  `# if defined LMD_MIXING` → `# if defined LMD_MIXING && !defined AKB0`.
  This restores the debug switch that was present in `src/wec_frc.F:2469`
  but was dropped when the override rewrote this section. Defining `AKB0`
  in `cppdefs.opt` now correctly disables the entire Akb (and Akw) mixing
  block. Updated matching `#endif` comment.

- **New BODY_BREAK flag** for breaking momentum as 3D body force.
  Replaced all `#ifndef SURFACE_BREAK` guards with
  `#if defined SURFACE_BREAK && defined BODY_BREAK` at 7 locations:
  - Declarations (line 131): `brk_r`, `brk3dx`, `brk3de`
  - Allocations (line 819): same arrays
  - `ru(k)` addition (line 1197): `+dm_u * brk3dx`
  - `rv(k)` addition (line 1229): `+dn_v * brk3de`
  - Body force computation (line 2200): `brk_r`, `brk3dx`, `brk3de` calculation
  - MPI exchange (line 2275): `exchange_xxx(brk3dx,brk3de)`
  - Commented-out exchange (line 2296): updated for consistency

  Design mirrors `WAVE_FRICTION` / `BODY_FRICTION` pattern:
  - `SURFACE_BREAK` alone → stress (default, brk2dx via sustr/top cell)
  - `SURFACE_BREAK` + `BODY_BREAK` → body force (brk3dx in ru(k))
  - Neither → no breaking momentum

- **Guarded sustr contamination** to prevent double-counting when `BODY_BREAK`
  is defined. Wrapped lines 1724 and 1749 (`sustr += brk2dx`, `svstr += brk2de`)
  with `#if defined SURFACE_BREAK && !defined BODY_BREAK`. Without this guard,
  defining `BODY_BREAK` would give 2x breaking: once through sustr and once
  through brk3dx in ru(k).

### cppdefs.opt (continued)
- Added `#undef BODY_BREAK` under `SURFACE_BREAK` with comment.
- **Reactivated BODY_FRICTION** streaming body force computation.
  The code was entirely commented out in src/wec_frc.F:2640-2739 and absent
  from the override. Added active implementation using FF_FUNC3 profile shape:
  - Added `frc_r` declaration and allocation under `#ifdef BODY_FRICTION`
  - Vertical profile computation: `frc_r(i,j,k) = fb * frc(i,j)` using
    FF_FUNC3 bottom-decaying exponential shape (lines 2260-2290)
  - Projection to u/v grid: `frc3dx(k)`, `frc3de(k)` (lines 2292-2312)
  - MPI exchange: `exchange_xxx(frc3dx,frc3de)` added
  - **UORBOFF switch** for near-bed orbital amplitude `abot`:
    - `UORBOFF` undef: `abot = wh/(2*sinh(kD))` (monochromatic, Soulsby)
    - `UORBOFF` defined: `abot = Uorb_mag * inv_f` (from WW3 orbital velocity)
  - `abot` → `delta` (WBBL thickness, U2010 Eq. 56) → vertical scale `kfrc`
  - FPE guard: `if (wh > eps)` with `else frc_r=0.0` (from Duck_F_opt pattern).
    Prevents `(abot/Zob)**0.83` trap with `-fpe0` when abot=0.
  - **Fixed fb1 initialization** relative to Duck_F_opt version: `fb1 = 1.0-exp(-2*kfrc*Dstp)`
    (from src/ FF_FUNC3), NOT `fb1 = 0` (which Duck_F_opt had incorrectly copied from
    the breaking profile initialization). With `fb1=0`, the streaming profile telescopes
    to ~0 and the body force is effectively dead.

### Duck_F_opt/wec_frc.F (cross-fix)
- **Fixed fb1 initialization bug** in BODY_FRICTION block (line 2317):
  `fb1 = 0.0` → `fb1 = 1.0-exp(-2.*kfrc*Dstp(i,j))`.
  The previous value was copied from the breaking body force code where fb1=0
  is correct (breaking profile is zero at bottom). For the streaming profile
  (bottom-concentrated, FF_FUNC3), the correct initial value at z=-h is
  `1-exp(-2kD)`. With fb1=0, the recursive integration telescoped to ~0,
  making BODY_FRICTION effectively dead even when defined.

- **Removed `#undef BBL`**. The override developer already merged all `BBL`
  functionality into `UORBOFF`: every `#ifdef BBL` in src/ (uorb/vorb
  declarations, wwuob/wwvob file arrays, allocations, initializations)
  was replaced by `#ifdef UORBOFF` in the override. The src/ guard
  `#if !defined BBL && defined BBL_S95` was also replaced with just
  `#if defined BBL_S95`. The `#undef BBL` line was vestigial — BBL has
  zero active references in the override.

## 2026-02-23

### cppdefs.opt
- **New `BREAKING` master key** for all breaking effects (momentum + mixing).
  Added `#define BREAKING` and nested all breaking sub-flags under `#ifdef BREAKING`:
  `SURFACE_BREAK`, `BODY_BREAK`, `BRK0`, `AKB0`, `Akb_U2010`, `FKV_FUNC4`.
  When `BREAKING` is undefined, ebrk=0, brk=0, Akb block skipped → zero breaking.
  Moved `BRK0` from standalone to nested under BREAKING. Added `AKB0` (previously
  only in wec_frc.F guard, now explicitly undef'd in cppdefs.opt).

### wec_frc.F
- **Guarded ebrk/brk computation with `BREAKING`** (lines 1629, 1641):
  `#ifdef BREAKING` / `#else ebrk=0.0` / `#endif` and same for brk.
  When BREAKING is undefined, both are explicitly set to zero, ensuring
  no breaking signal enters the system.

- **Added `BREAKING` to all breaking CPP guards** — 12 locations:
  - Declarations (line 131): `BREAKING && SURFACE_BREAK && BODY_BREAK`
  - Allocations (line 820): same
  - ru(k) addition (line 1200): same
  - rv(k) addition (line 1232): same
  - BRK0 block (line 1662): `BREAKING && BRK0`
  - sustr contamination (lines 1732, 1759): `BREAKING && SURFACE_BREAK && !BODY_BREAK`
  - Akb block (line 2071): `BREAKING && LMD_MIXING && !AKB0`
  - Body force computation (line 2215): `BREAKING && SURFACE_BREAK && BODY_BREAK`
  - MPI exchange (line 2353): same
  - Commented-out exchange (line 2377): same (consistency)

- **Fixed kbrk bug in body force block** (line 2226): added `kbrk = 1.0/inv_k`
  after `inv_k` computation. Previously, `kbrk` was only set inside the
  `Akb_U2010` block (Akb mixing). The body force block reused the stale
  scalar from the Akb loop's last iteration — wrong for all (i,j) except
  the last. Same fix as Duck_F_opt (from MPI periodic debug session).

- **Guarded frc computation with `WAVE_FRICTION`** (line 1639):
  `#ifdef WAVE_FRICTION` / `#else frc=0.0` / `#endif`.
  When WAVE_FRICTION is undefined, frc is explicitly set to zero, ensuring
  no streaming signal enters the system (frc2dx/frc2de naturally zero from frc=0).
  Mirrors the BREAKING guard for ebrk/brk.

### lmd_kpp.F
- **Guarded Cr(z) streaming term with `WAVE_FRICTION` and `FRC_HBBL0`** (line 269):
  Changed `#ifdef WEC` to `#if defined WEC && defined WAVE_FRICTION && !defined FRC_HBBL0`.
  The C_Ek*f²*|z|/H streaming stabilization term (D. Hypolite addition) in the
  bottom boundary layer detection is now:
  - Active by default when WAVE_FRICTION is defined (current behavior)
  - Disabled when WAVE_FRICTION is undefined (no streaming → no BBL effect)
  - Explicitly zeroed when FRC_HBBL0 is defined (debug switch, like AKB0/BRK0)

### cppdefs.opt (continued)
- Nested streaming sub-flags under `#ifdef WAVE_FRICTION`:
  `BODY_FRICTION` (undef=bottom cell default, defined=body force) and
  `FRC_HBBL0` (undef=active BBL streaming, defined=zero BBL streaming).
  BODY_FRICTION was previously implicit (not listed in cppdefs.opt);
  now explicitly undef'd for clarity.

### lmd_kpp.F (continued)
- **New `WMB0` flag** for wave-enhanced bottom friction velocity (line 465):
  Changed `#if defined WEC && defined BBL_S95` to
  `#if defined WEC && defined BBL_S95 && !defined WMB0`.
  When `WMB0` is defined, the BKPP friction velocity `wmb` falls back to
  the standard current-only formula `vonKar²*cff/log(1+0.5*Hz/Zob)`
  instead of the wave-enhanced `vonKar*sqrt(r_D*cff)`.
  This decouples: BBL_S95 still enhances momentum drag (r_D in
  compute_rd_bott_drag.h), but the enhanced r_D no longer feeds into
  BKPP mixing depth/intensity. Active by default (undef'd).

### wec_frc.F (continued)
- **Fixed missing `WAVE_DIFF_FROM_LM` guard** (line 1470): wrapped
  `llmm(i,j) = wwlm(i,j) * rmask(i,j)` with `#ifdef WAVE_DIFF_FROM_LM`.
  Previously only guarded by `#ifndef MONOWAVE`. When `WAVE_DIFF_FROM_LM`
  is undef'd, `wwlm` is not declared/allocated → compile error.
  Pre-existing bug exposed by our `#undef WAVE_DIFF_FROM_LM` change.

### cppdefs.opt (continued)
- Added `#undef WMB0` nested under `#ifdef BBL_S95` with comment.

### wec_frc.opt
- Added 4 output flags for 3D body force fields (all default `.false.`):
  `wrt_brk3dx`, `wrt_brk3de` (need `BODY_BREAK`),
  `wrt_frc3dx`, `wrt_frc3de` (need `BODY_FRICTION`).

### wec_frc.F (continued)
- **Added 3D body force NetCDF output** for `brk3dx`, `brk3de`, `frc3dx`, `frc3de`:
  - Name/long_name/units strings (after brk2de strings)
  - `_avg` declarations with CPP guards (`BREAKING && SURFACE_BREAK && BODY_BREAK`
    for brk3d*, `WAVE_FRICTION && BODY_FRICTION` for frc3d*)
  - `_avg` allocations with same CPP guards
  - `ncreate` blocks with same CPP guards (xi_u/eta_rho for *x, xi_rho/eta_v for *e)
  - Averaging accumulation in `calc_avg_wec` with same CPP guards
  - History write in `wrt_wec_his` with same CPP guards
  - Average write in `wrt_wec_avg` with same CPP guards

- **Fixed `wst_avg` accumulation bug** (line 2877): `vst_avg`/`vst` → `wst_avg`/`wst`.
  Pre-existing copy-paste error: wst averages were being accumulated from vst data.

- **Fixed `brk2de` avg write bug** (line 3104): `frc2de_name` → `brk2de_name`.
  Pre-existing copy-paste error: brk2de averages were written under the frc2de variable
  name in the NetCDF file.

- **New `FRC0` flag** for zeroing streaming momentum (mirrors `BRK0` for breaking).
  Added `#if defined WAVE_FRICTION && defined FRC0` block after `BRK0` block:
  sets `frc(i,j)=0.D0` over the tile. When defined, streaming momentum (frc2dx/frc2de)
  is zero while streaming Cr(z) in BBL detection and Akw mixing remain active.

### cppdefs.opt (continued)
- Added `#undef FRC0` nested under `#ifdef WAVE_FRICTION` with comment.
