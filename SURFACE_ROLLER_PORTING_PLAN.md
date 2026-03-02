# Surface Roller Porting Plan: WKB → WAVE_OFFLINE

## Context

The surface roller model (U2010 Eqs. 40–45, 48, 59) is fully implemented in the
older WKB online wave model (`ROMS_DHypolite/wkb_wwave.F.old`) but was never
ported to the `WAVE_OFFLINE` pathway used by the current code with WW3 coupling.

The roller physically spreads breaking momentum shoreward: a fraction α_r of
ε^b feeds the roller, which propagates at phase speed c and dissipates as ε^r.
It is ε^r (not ε^b) that drives the nearshore currents and breaking-induced
mixing (Eq. 48, 59). This provides a physically motivated alternative to
`SMOOTH_BRK` (isotropic Laplacian smoothing).

**Existing infrastructure** (already in `DEBUG_L4PTSAL/wec_frc.F`):
- `rdsp(GLOBAL_2D_ARRAY)` — allocated under `#ifdef SURFACE_ROLLER`, currently unused
- `erol(GLOBAL_2D_ARRAY)` — allocated as module-level array, currently unused
- `wvqb` / `wwqb` — Q_b forcing arrays, allocated under `SURFACE_ROLLER`
- `nc_wwqb` — forcing reader for qb, defined in `wec_frc.opt`
- `set_frc_data(nc_wwqb, wwqb, 'r')` — already reads qb when `SURFACE_ROLLER` defined
- `roller_fraction = 0.0` — parameter in `wec_forces_tile`, currently zero

**Forcing file**: Contains `qb` (fraction of breaking waves). Also provides ε^b
(`eb`), wave amplitude, direction, period, and orbital velocities.

**Reference equations** (U2010 = Uchiyama et al. 2010):
- Eq. 40: roller action balance (prognostic PDE)
- Eq. 41: roller dissipation rate
- Eq. 48: modified breaking body force (stress)
- Eq. 50: depth-averaged breaking acceleration
- Eq. 59: modified breaking-induced mixing K^b_v

See `surface_roller.tex` for the full pedagogical derivation.

---

## Equation Summary

### Roller action balance (U2010 Eq. 40)

```
∂A^r/∂t + ∇·(A^r c) = (α_r ε^b − ε^r) / σ
```

where:
- `A^r = E^r/σ` is roller action density [m³/s] (code units: divided by ρ₀)
- `c = σ/k` is phase speed (NOT group speed)
- `α_r` is roller fraction (0 to 1; U2010 baseline: α_r = 1.0)
- `ε^b` is breaking dissipation [m³/s³] (= `ebrk` = `wveb` in code, already /ρ₀)

### Roller dissipation (U2010 Eq. 41)

```
ε^r = g sin(β) E^r / c = g sin(β) σ A^r / (σ/k) = g sin(β) A^r k σ
```

In code units (divided by ρ₀ and σ for action):
```
rdsp = g sin(β) war_code k     (ε^r/ρ₀, [m³/s³])
```
where `war_code = A^r/ρ₀` [m³/s].

### Modified breaking stress (U2010 Eq. 48/50)

```
brk(i,j) = [(1 − α_r) ebrk + rdsp] × kw × inv_f
```

### Modified Akb (U2010 Eq. 59)

```
K^b_v = c_b × ((1 − α_r) ebrk + rdsp)^(1/3) × H_* × D × f^Kv(z)
```

Replace `ebrk` → `(1 − α_r) ebrk + rdsp` in the existing Akb block.

---

## Implementation Strategy

### Approach: Prognostic roller with forward Euler

Solve the roller action balance (Eq. 40) at each baroclinic time step using
1st-order upwind advection (identical scheme to `wkb_wwave.F.old`).

**Why prognostic (not diagnostic)?** The diagnostic formula (Eq. 45 with Q_b)
gives a local estimate of roller energy based on WW3's Q_b field. But Q_b peaks
at the breaking point — the diagnostic roller would NOT be displaced shoreward.
The whole point of the roller is the spatial propagation via ∇·(A^r c), which
only the prognostic PDE captures.

**Why not use Q_b in the prognostic equation?** Q_b is not needed by Eq. 40 —
the source is α_r × ε^b. However, Q_b is useful for:
1. Diagnostic initialization (avoid cold-start transient)
2. Validation (compare prognostic ε^r against α_r × ε^b)
3. Future: modulated α_r (e.g., α_r × Q_b, only create roller where breaking)

### State variable

**Repurpose `erol`** (already allocated as `GLOBAL_2D_ARRAY`) to store the roller
action density `A^r/ρ₀` [m³/s]. This is persistent between time steps.

Currently `erol` is allocated but never used. With the roller, it becomes the
prognostic roller state variable.

### Time stepping

`wec_forces_tile` is called **twice** per baroclinic step (predictor + corrector).
The roller should advance only **once** per step. Use a guard:

```fortran
if (PREDICTOR_2D_STEP(iic)) then
   ! advance roller: erol^{n+1} from erol^n
endif
! always compute rdsp from current erol for brk and Akb
```

Forward Euler is adequate because:
- The roller equilibrium time is O(10 s) (surfzone width / phase speed)
- The baroclinic step is O(3 s)
- 1st-order upwind already introduces numerical diffusion that damps instabilities

### CFL safety

Phase speed c = σ/k. In surfzone (D < 10 m), c < 10 m/s.
With Δx = 30 m, dt = 3 s: CFL = 10×3/30 = 1.0 (marginal).

**Safety measure**: clamp `erol` to non-negative after each step.
If CFL > 1 in deeper water, it doesn't matter because ε^b → 0 there and
the roller decays exponentially regardless. Add optional subcycling if needed.

---

## File-by-File Changes

### 1. `cppdefs.opt` — Add CPP flag

Under `#ifdef BREAKING`, add:
```fortran
#    define SURFACE_ROLLER   /* prognostic surface roller (U2010 Eq.40) */
```

This activates the existing `rdsp`, `wvqb`/`wwqb` declarations and allocations
that are already ifdef'd throughout `wec_frc.F`.

### 2. `wec_frc.opt` — Roller parameters

Add roller configuration parameters:
```fortran
#ifdef SURFACE_ROLLER
      real, parameter :: alpha_r = 1.0     ! roller fraction (U2010: 0-1)
      real, parameter :: sin_beta = 0.1    ! roller slope parameter sin(β)
#endif
```

### 3. `wec_frc.F` — Module-level arrays

**No new arrays needed.** The existing infrastructure is sufficient:
- `erol(GLOBAL_2D_ARRAY)` — repurpose as prognostic A^r/ρ₀
- `rdsp(GLOBAL_2D_ARRAY)` — roller dissipation ε^r/ρ₀ (already under SURFACE_ROLLER)
- `wvqb` / `wwqb` — Q_b from WW3 (already under SURFACE_ROLLER)

Add scratch arrays for advection fluxes:
```fortran
#ifdef SURFACE_ROLLER
      real, allocatable, dimension(:,:) :: FX_rol, FE_rol
#endif
```

Allocate in `init_wec_arrays`:
```fortran
#ifdef SURFACE_ROLLER
      allocate( FX_rol(GLOBAL_2D_ARRAY)); FX_rol=0.
      allocate( FE_rol(GLOBAL_2D_ARRAY)); FE_rol=0.
#endif
```

### 4. `wec_frc.F` — Local variables in `wec_forces_tile`

```fortran
# ifdef SURFACE_ROLLER
      real c_phase, cpx, cpy, cxu, cyv
      real alpha_r_loc, sin_beta_loc, dt_rol
#  include "wec_frc.opt"          ! reads alpha_r, sin_beta
# endif
```

Note: `alpha_r_loc` / `sin_beta_loc` are local copies to avoid namespace
conflicts. Or simply use the parameter names from `wec_frc.opt` directly.

### 5. `wec_frc.F` — Roller time stepping (new include file)

After the ebrk/brk computation (after `MASKING` block at ~line 1700) and
BEFORE the BRK0 block, insert:

```fortran
# if defined BREAKING && defined SURFACE_ROLLER
#  include "step_roller.h"
# endif
```

This must come **after** `ebrk`, `kw`, `fr`, `wdrx`, `wdre` are computed
(they are inputs to the roller), and **before** the `brk(i,j)` line is
modified to use `rdsp`.

### 6. `wec_frc.F` — Modified `brk` computation

Replace (line ~1691):
```fortran
          brk(i,j) =cff*kw(i,j)*inv_f(i,j)       ! 2DH breaking term
```
with:
```fortran
# ifdef SURFACE_ROLLER
          brk(i,j) =((1.0-alpha_r)*cff + rdsp(i,j))
     &                    *kw(i,j)*inv_f(i,j)     ! 2DH breaking + roller
# else
          brk(i,j) =cff*kw(i,j)*inv_f(i,j)       ! 2DH breaking term
# endif
```
where `cff = ebrk(i,j)`.

### 7. `wec_frc.F` — Modified Akb computation

In the Akb block (~line 2071), replace `ebrk(i,j)` with the combined
dissipation. Currently:
```fortran
          cff1    =ebrk(i,j)**c1o3
```
Change to:
```fortran
# ifdef SURFACE_ROLLER
          cff1    =((1.0-alpha_r)*ebrk(i,j)+rdsp(i,j))**c1o3
# else
          cff1    =ebrk(i,j)**c1o3
# endif
```

The FPE guard (`if (wh > eps)`) already protects the `**c1o3` operation.

### 8. `wec_frc.F` — Roller Stokes drift (optional, lower priority)

In the Stokes drift computation, replace `act(i,j)` with
`act(i,j) + erol(i,j)/rho0` to include the roller contribution
(U2010 Eq. 43-44). This is a secondary effect and can be deferred.

### 9. Create `DEBUG_L4PTSAL/step_roller.h` (NEW FILE)

Prognostic roller advection + source/sink, adapted from
`ROMS_DHypolite/wkb_wwave.F.old` lines 309-356.

```fortran
! step_roller.h — Prognostic surface roller (U2010 Eq. 40)
! Advances erol (= A^r/ρ₀) by one baroclinic time step.
! Uses 1st-order upwind advection at phase speed c = σ/k.
! Source: α_r × ε^b/σ.  Sink: g sin(β) × A^r × k / σ.
!
! Called from wec_forces_tile after ebrk/kw/fr/wdrx/wdre computed,
! before brk(i,j) computation.

      if (PREDICTOR_2D_STEP(iic)) then

        dt_rol = dt        ! baroclinic time step

! Compute phase velocity components at rho-points.
! c = σ/k = fr/kw;  cx = c × wdrx;  cy = c × wdre.

        do j=jstr-1,jend+1
          do i=istr-1,iend+1
            c_phase = fr(i,j) / max(kw(i,j), eps)
            cpx_arr(i,j) = c_phase * wdrx(i,j)
            cpy_arr(i,j) = c_phase * wdre(i,j)
          enddo
        enddo

! MPI exchange of erol before advection (fill ghost points).

#  ifdef EXCHANGE
        call exchange_xxx(erol)
#  endif

! Neumann BC at physical boundaries.

#  ifndef EW_PERIODIC
        if (WESTERN_EDGE) then
          do j=jstr-1,jend+1
            erol(istr-1,j) = erol(istr,j)
          enddo
        endif
        if (EASTERN_EDGE) then
          do j=jstr-1,jend+1
            erol(iend+1,j) = erol(iend,j)
          enddo
        endif
#  endif
#  ifndef NS_PERIODIC
        if (SOUTHERN_EDGE) then
          do i=istr-1,iend+1
            erol(i,jstr-1) = erol(i,jstr)
          enddo
        endif
        if (NORTHERN_EDGE) then
          do i=istr-1,iend+1
            erol(i,jend+1) = erol(i,jend)
          enddo
        endif
#  endif

! 1st-order upwind advection fluxes (wkb_wwave.F pattern).
! FX_rol = flux in xi-direction at u-point (i+1/2,j).
! FE_rol = flux in eta-direction at v-point (i,j+1/2).

        do j=jstr,jend
          do i=istr-1,iend
            cxu = 0.5*(cpx_arr(i+1,j) + cpx_arr(i,j))
            FX_rol(i,j) = 0.5*dn_u(i+1,j)*(
     &          (cxu - abs(cxu))*erol(i+1,j)
     &        + (cxu + abs(cxu))*erol(i,j)   )
          enddo
        enddo
        do j=jstr-1,jend
          do i=istr,iend
            cyv = 0.5*(cpy_arr(i,j+1) + cpy_arr(i,j))
            FE_rol(i,j) = 0.5*dm_v(i,j+1)*(
     &          (cyv - abs(cyv))*erol(i,j+1)
     &        + (cyv + abs(cyv))*erol(i,j)   )
          enddo
        enddo

! Update erol: advection + source (α_r ε^b/σ) − sink (g sinβ A^r k/σ)
! Note: ebrk = ε^b/ρ₀, erol = A^r/ρ₀, fr = σ, kw = k.

        do j=jstr,jend
          do i=istr,iend
            erol(i,j) = erol(i,j) - dt_rol * pm(i,j)*pn(i,j)
     &          * ( FX_rol(i,j) - FX_rol(i-1,j)
     &            + FE_rol(i,j) - FE_rol(i,j-1) )
     &        + dt_rol * ( alpha_r * ebrk(i,j)
     &                   - g * sin_beta * erol(i,j) * kw(i,j) )
     &                   / max(fr(i,j), eps)
#  ifdef MASKING
            erol(i,j) = erol(i,j) * rmask(i,j)
#  endif
            erol(i,j) = max(erol(i,j), 0.0)   ! non-negativity
          enddo
        enddo

      endif  ! PREDICTOR_2D_STEP

! Compute roller dissipation rdsp = ε^r/ρ₀ = g sinβ × erol × kw × fr
! (always, for both predictor and corrector).

      do j=jstr-1,jend+1
        do i=istr-1,iend+1
          rdsp(i,j) = g * sin_beta * erol(i,j) * kw(i,j)
     &                * fr(i,j)
#  ifdef MASKING
     &                * rmask(i,j)
#  endif
        enddo
      enddo
```

**Note**: The above is a template. The exact loop bounds, scratch array names
(`cpx_arr`, `cpy_arr`), and `dn_u`/`dm_v` metric terms need to match the
override's conventions. `cpx_arr`/`cpy_arr` could reuse existing scratch
arrays (e.g., `FX_brk`, `FE_brk` if SMOOTH_BRK is disabled, or add new ones).

### 10. `cppdefs.opt` — Disable SMOOTH_BRK when roller is active

```fortran
#  ifdef BREAKING
#    define SURFACE_ROLLER   /* prognostic surface roller (U2010 Eq.40) */
#   ifndef SURFACE_ROLLER
#    define SMOOTH_BRK       /* Laplacian smoothing (surrogate for roller) */
#   endif
```

Or simply `#undef SMOOTH_BRK` when `SURFACE_ROLLER` is defined.
The two should NOT be used simultaneously.

---

## Diagnostic Initialization (Optional)

To avoid the cold-start transient (erol = 0 → takes O(100/c) seconds to
spin up), initialize erol from the diagnostic formula (U2010 Eq. 45):

```fortran
! In init_wec_arrays or first call to wec_forces_tile:
if (FIRST_TIME_STEP) then
  do j=...; do i=...
    ! Diagnostic roller action: A^r/ρ₀ ≈ g D A_R Q_b / (2 L_p σ_p)
    ! Simplified: use ebrk as proxy
    if (ebrk(i,j) > eps .and. kw(i,j) > eps) then
      c_phase = fr(i,j) / kw(i,j)
      erol(i,j) = alpha_r * ebrk(i,j) / (g * sin_beta * kw(i,j))
    endif
  enddo; enddo
endif
```

This gives the steady-state local balance `α_r ε^b = ε^r = g sinβ erol k σ`,
i.e. `erol = α_r ε^b / (g sinβ k σ)`.

---

## Interaction with Existing Flags

| Flag | With roller | Notes |
|------|-------------|-------|
| `SURFACE_ROLLER` | **NEW** master switch | Activates roller PDE + modified brk/Akb |
| `SMOOTH_BRK` | `#undef` | Redundant when roller active |
| `BRK0` | still works | Zeros brk AFTER roller modifies it |
| `AKB0` | still works | Zeros Akb AFTER roller modifies ebrk input |
| `BODY_BREAK` | still works | Body force uses brk(i,j) which includes roller |
| `SURFACE_BREAK` | still works | Stress uses brk2dx/de which derives from brk(i,j) |
| `alpha_r = 0` | degenerates | No roller, erol→0, rdsp→0, brk = ebrk×k/σ (original) |

---

## Testing Plan

### 1. Compilation
- Define `SURFACE_ROLLER`, verify clean compilation with Intel ifx
- Check that `SURFACE_ROLLER` + `BRK0` still works (roller stepping but brk zeroed)

### 2. α_r = 0 regression
- Set `alpha_r = 0.0` → source term is zero → erol stays at 0 → rdsp = 0
- brk and Akb should match the no-roller run exactly

### 3. Sensitivity sweep (matches U2010 Runs 1–7, Table 2)
- α_r = {0, 0.25, 0.5, 0.75, 1.0}
- Compare cross-shore profiles of ε^tot, v̄, u(x,z)
- Expect: increasing α_r shifts peak v̄ shoreward and reduces v̄ gradient
  (U2010 Fig. 6)

### 4. Conservation check
- Print area-integrated (1-α_r)ebrk + rdsp vs area-integrated ebrk
- In steady state, should match (roller conserves total momentum)
- Transient: total = (1-α_r)∫ebrk + ∫rdsp ≤ ∫ebrk (some in roller storage)

### 5. L4PTSAL comparison
- Run A: BRK0 (no breaking) — baseline
- Run B: SMOOTH_BRK N=4 (Laplacian smoothing) — current fix
- Run C: SURFACE_ROLLER α_r=1.0 — physical roller
- Compare nearshore currents, especially in embayments/headlands
- Check CFL: monitor max(c × dt / dx) during run

### 6. Akb check
- Compare K^b_v profiles with and without roller
- With roller, Akb should be slightly broader and shifted shoreward
  (U2010 Fig. 3)

---

## Summary of New/Modified Files

| File | Action | Lines (est.) |
|------|--------|--------------|
| `cppdefs.opt` | Add `SURFACE_ROLLER`, restructure `SMOOTH_BRK` | ~5 |
| `wec_frc.opt` | Add `alpha_r`, `sin_beta` parameters | ~4 |
| `wec_frc.F` | Add scratch arrays, local vars, include, modify brk/Akb | ~30 |
| `step_roller.h` | **NEW**: prognostic roller PDE solver | ~100 |
| `CHANGELOG.md` | Document all changes | ~40 |

Total: ~180 lines of new/modified Fortran.

---

## References

- Uchiyama, Y., J.C. McWilliams, and A.F. Shchepetkin (2010). Wave–current
  interaction in an oceanic circulation model with a vortex-force formalism:
  Application to the surf zone. *Ocean Modelling*, 34, 16–35.
  - Eq. 40: roller action balance
  - Eq. 41: roller dissipation
  - Eq. 48: breaking body force with roller
  - Eq. 50: depth-averaged breaking acceleration with roller
  - Eq. 59: breaking-induced mixing K^b_v with roller
  - Table 1: sin β = 0.1, α_r = 1.0
  - Table 2: Run 1 (α_r = 1, best fit to DUCK94)
  - Fig. 6: effect of α_r on ε^tot and v̄ profiles

- Svendsen, I.A. (1984a). Mass flux and undertow in a surf zone.
  *Coastal Engineering*, 8, 347–365. (Diagnostic roller formula, Eq. 45)

- Nairn, R.B., J.A. Roelvink, and H.N. Southgate (1991). Transition zone
  width and implications for modelling surfzone hydrodynamics. Proc. 22nd
  ICCE, ASCE, 68–81. (Roller model formulation)
