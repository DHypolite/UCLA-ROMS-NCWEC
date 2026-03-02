# Plan: Depth-Integrated Momentum Budget Output

## Context

Diagnose the depth-integrated momentum tendency equation term by term:

```
∂t ∫u dz = ∫B^b dz + ∫B^wd dz + τ_wind + τ^cd_bot + ∫∇⊥ϕ dz + ∫J dz + ∫D⊥ dz
```

where D⊥ = Coriolis + Eulerian advection + lateral diffusion.

This requires outputting each term as a 2D (depth-integrated) field. Some terms are
already available, others need new code.

## Term-by-Term Assessment

| # | Term | Meaning | Code variable | Status |
|---|------|---------|---------------|--------|
| 1 | ∂t ∫u dz | Transport tendency | d(ubar·D)/dt | Offline from ubar, H |
| 2 | ∫B^b dz | Breaking momentum | `brk2dx`, `brk2de` | **Already output** (2D) |
| 3 | ∫B^wd dz | Wave drag/streaming | `frc2dx`, `frc2de` | **Already output** (2D) |
| 4 | τ_wind | Wind stress | `sustr - brk2dx` | **Recoverable** (see below) |
| 5 | τ^cd_bot | Bottom drag stress | `r_D · u(k=1)` | **Needs new output** |
| 6 | ∫∇⊥ϕ dz | Pressure gradient | DIAGNOSTICS `ipgr` | **Needs DIAGNOSTICS** |
| 7 | ∫J dz | Vortex force | horizontal VF | **Needs new accumulator** |
| 8 | ∫D⊥ dz | Cor + adv + diff | DIAGNOSTICS terms | **Needs DIAGNOSTICS** |

### Key findings from code exploration:

1. **Wind stress is recoverable without code changes**: In SURFACE_BREAK mode,
   `sustr = τ_wind + brk2dx` (wec_frc.F:1854). Since `brk2dx` is separately output,
   `τ_wind = sustr_output - brk2dx_output`.

2. **The existing DIAGNOSTICS module** (`diagnostics.F`) provides a 7-term 3D momentum
   budget: ipgr, icori, iadv, idiss, ihmix, ivmix, icoup. Units are `du/dt × dz` [m²/s²]
   (already vertically integrated per cell). Depth integral = Σ_k Udiag(i,j,k,n).

3. **`iadv` is computed independently** via 4th-order advective fluxes
   (`set_diags_u_4th_adv` in diagnostics.F:283-382). It does NOT include WEC vortex
   force — that goes into `idiss`.

4. **`idiss`** is a closure residual: `idiss = total_ru - ipgr - icori - iadv`.
   For the current config (SURFACE_BREAK, bottom-cell streaming), it captures:
   numerical dissipation + horizontal vortex force + frc2dx at k=1.

5. **`ihmix`** is set in `visc3d_S_FS.F` (horizontal viscosity, separate from idiss).

6. **`ivmix`** captures the implicit vertical solver effect:
   `ivmix = (u_new - u_prev)/dt - total_ru`. Depth-integrated, this ≈ sustr - τ_bot.

7. **`r_D`** (bottom drag coefficient) is computed in `compute_rd_bott_drag.h`,
   stored in shared memory, but NOT currently output.

## Mapping: User Terms ↔ Diagnostic Terms

```
∫∇⊥ϕ dz  =  Σ_k ipgr(k)          ← total PGF including WEC corrections (kvf, calP, Kapsrf)
∫D⊥ dz   =  Σ_k (icori+iadv+ihmix)(k)  ← Coriolis + advection + h-mixing (clean, no WEC)
∫J dz    =  NEW custom accumulator  ← horizontal vortex force saved separately
∫B^b dz  =  brk2dx/brk2de          ← already output
∫B^wd dz =  frc2dx/frc2de          ← already output
τ_wind   =  sustr - brk2dx         ← offline subtraction
τ^cd_bot =  r_D · u(k=1)           ← from r_D output + u(k=1) output, offline
```

**Budget closure check from diagnostics:**
```
Σ_k (all 7 terms) = total tendency
Residual = Σ_k (idiss + ivmix + icoup) - ∫J - frc2dx - sustr + τ_bot ≈ numerical_diss
```

## Implementation Plan

### Step 1: Enable DIAGNOSTICS (`cppdefs.opt`)

Add after the WEC block (around line 60):
```fortran
#define DIAGNOSTICS
```

### Step 2: Create `diagnostics.opt` override

Copy `src/diagnostics.opt` → `DEBUG_L4PTSAL/diagnostics.opt`. Changes:
```fortran
logical, parameter, public :: diag_uv  = .true.    ! Enable momentum diagnostics
integer, parameter         :: output_period = 900   ! Match avg output period (seconds)
logical, parameter         :: diag_avg = .true.     ! Time-averaged output
```

### Step 3: Add depth-integration to `diagnostics.F`

Copy `src/diagnostics.F` → `DEBUG_L4PTSAL/diagnostics.F`. Changes:

**a) Module-level arrays** (after line 53, Udiag_avg declaration):
```fortran
real, allocatable, dimension(:,:,:), public :: Udiag_int  ! depth-integrated u diags (nx,ny,nmd)
real, allocatable, dimension(:,:,:), public :: Vdiag_int  ! depth-integrated v diags (nx,ny,nmd)
```

**b) Allocation** (in `init_diagnostics`, after Udiag allocation ~line 186):
```fortran
allocate(Udiag_int(nx,ny,nmd)); Udiag_int = 0
allocate(Vdiag_int(nx,ny,nmd)); Vdiag_int = 0
```

**c) Depth-integration subroutine** (new, ~25 lines):
```fortran
subroutine calc_diag_int
  ! Sum Udiag/Vdiag over k to get depth-integrated budget terms.
  ! Units: m^2/s^2 (du/dt × dz, vertically integrated).
  integer :: i,j,k,n
  Udiag_int = 0
  Vdiag_int = 0
  do n = 1, nmd
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          Udiag_int(i,j,n) = Udiag_int(i,j,n) + Udiag(i,j,k,n)
          Vdiag_int(i,j,n) = Vdiag_int(i,j,n) + Vdiag(i,j,k,n)
        enddo
      enddo
    enddo
  enddo
end subroutine
```

**d) Call depth-integration** in `do_diagnostics` (before writing, ~line 922):
```fortran
call calc_diag_int
```

**e) NetCDF variable creation** (in `create_diagvars`, add 14 new 2D vars):
- `u_pgr_int`, `u_cor_int`, `u_adv_int`, `u_dis_int`, `u_hmx_int`, `u_vmx_int`, `u_cpl_int`
- `v_pgr_int`, `v_cor_int`, `v_adv_int`, `v_dis_int`, `v_hmx_int`, `v_vmx_int`, `v_cpl_int`
- Dimensions: (xi_u/xi_v, eta_u/eta_v, time) — 2D + time

**f) NetCDF writing** (new subroutine `wrt_diag_uv_int`, ~20 lines):
Write Udiag_int(:,:,1:7) and Vdiag_int(:,:,1:7) to the diagnostics file.

**g) If diag_avg**: also accumulate running averages of Udiag_int/Vdiag_int
(in `calc_diag_avg`), write averaged versions.

### Step 4: Add vortex force accumulator in `wec_frc.F`

**a) Module-level 2D arrays** (after shelter_fac declaration, ~line 201):
```fortran
#if defined WEC && defined BREAKING
      real, allocatable, dimension(:,:) :: hvf_int_xi   ! ∫J_xi dz
      real, allocatable, dimension(:,:) :: hvf_int_eta  ! ∫J_eta dz
#endif
```

**b) Allocations** in `init_wec_arrays` (~line 790):
```fortran
allocate(hvf_int_xi(GLOBAL_2D_ARRAY)); hvf_int_xi = 0.
allocate(hvf_int_eta(GLOBAL_2D_ARRAY)); hvf_int_eta = 0.
```

**c) Accumulate in `compute_wec_rhs_uv_terms`** (~line 1270):
The horizontal vortex force is computed as `cff1` (u-component) and a similar
expression for v. Currently `cff1` is added to `ru` along with breaking/streaming
body forces. We need to save the vortex force contribution separately.

In the k-loop, BEFORE adding cff1 to ru:
```fortran
! Save depth-integrated horizontal vortex force
hvf_int_xi(i,j) = hvf_int_xi(i,j) + cff * cff1   ! cff includes Hz*dn volume factor
```
And similarly for the v-component.

**Important**: Reset `hvf_int_xi = 0` at the START of `compute_wec_rhs_uv_terms`
(before the k-loop). The subroutine is called twice per timestep (predictor + corrector),
so we keep only the last call's values (matching the diagnostics timing).

**d) Add to WEC output** (`wrt_wec_his` and `wrt_wec_avg`):
```fortran
if (wrt_hvf_int) call ncwrite(ncid, hvf_int_xi_name, hvf_int_xi, start_2d)
if (wrt_hvf_int) call ncwrite(ncid, hvf_int_eta_name, hvf_int_eta, start_2d)
```

### Step 5: Output r_D for offline bottom drag stress

**a) In `wec_frc.opt`**, add:
```fortran
logical, parameter :: wrt_rD     = .true.   ! bottom drag coefficient r_D [m/s]
logical, parameter :: wrt_hvf_int = .true.  ! depth-integrated horiz. vortex force
```

**b) In `wrt_wec_his`/`wrt_wec_avg`**, add r_D output.
r_D is stored in shared memory (via `ext_copy_prv2shr_2d_tile`). Need to verify
it's accessible from the WEC output subroutine. If not, copy it into a WEC module
array during `wec_forces_tile`.

**c) Offline τ_bot computation** (Python):
```python
tau_bot_u = 0.5 * (r_D[:, 1:] + r_D[:, :-1]) * u[:, 0, :, :]  # k=0 = bottom
tau_bot_v = 0.5 * (r_D[1:, :] + r_D[:-1, :]) * v[:, 0, :, :]
```

### Step 6: NetCDF variable definitions for WEC output

In `def_wec_his` / `def_wec_avg` (or the equivalent section of wec_frc.F):
- `hvf_int_xi`: "Depth-integrated horizontal vortex force, xi-component" [m^2/s^2]
- `hvf_int_eta`: "Depth-integrated horizontal vortex force, eta-component" [m^2/s^2]
- `r_D`: "Bottom drag coefficient" [m/s]

## Files Modified

| File | Action | Scope |
|------|--------|-------|
| `cppdefs.opt` | Edit (existing override) | Add `#define DIAGNOSTICS` |
| `diagnostics.opt` | New override (copy from src/) | Enable `diag_uv`, set periods |
| `diagnostics.F` | New override (copy from src/) | Add depth-integration arrays, subroutine, NetCDF output |
| `wec_frc.F` | Edit (existing override) | Add hvf_int accumulator, r_D output |
| `wec_frc.opt` | Edit (existing override) | Add `wrt_rD`, `wrt_hvf_int` flags |

## Post-Processing Budget (Python pseudocode)

```python
import xarray as xr

# Load outputs
his = xr.open_dataset('*_his.nc')      # history file
dia = xr.open_dataset('*_dia.nc')      # diagnostics file
wec = xr.open_dataset('*_wec_his.nc')  # WEC history file

# Already-available 2D terms:
Bb   = wec.brk2dx          # ∫B^b dz (breaking)
Bwd  = wec.frc2dx          # ∫B^wd dz (streaming)

# Wind stress (recover by subtracting breaking from total sustr):
tau_wind = his.sustr - wec.brk2dx

# Bottom drag stress (offline from r_D and bottom velocity):
rD = wec.r_D
u_bot = his.u.isel(s_rho=0)
tau_bot = -0.5 * (rD + rD.shift(xi_rho=-1)) * u_bot  # negative = opposes motion

# From depth-integrated diagnostics:
PGF  = dia.u_pgr_int       # ∫∇⊥ϕ dz (total PGF incl. WEC corrections)
Dperp = dia.u_cor_int + dia.u_adv_int + dia.u_hmx_int  # ∫D⊥ dz
J    = wec.hvf_int_xi       # ∫J dz (horizontal vortex force)

# Transport tendency (offline):
H = his.zeta + his.h
ubar_H = his.ubar * 0.5*(H + H.shift(xi_rho=-1))
ddt_uH = ubar_H.diff('ocean_time') / ubar_H.ocean_time.diff('ocean_time')

# Budget closure check:
residual = ddt_uH - (Bb + Bwd + tau_wind + tau_bot + PGF + J + Dperp)
# Should be small (numerical dissipation + coupling errors)

# Cross-check with full diagnostics sum:
total_diag = sum(dia[f'u_{t}_int'] for t in ['pgr','cor','adv','dis','hmx','vmx','cpl'])
# total_diag should equal ddt_uH
```

## Verification

1. **Compile test**: Build with `#define DIAGNOSTICS` and `diag_uv = .true.`
2. **Short run**: 1-2 output periods, check that diagnostic file is created with all variables
3. **Budget closure**: Verify `Σ (7 diag terms)_int ≈ ∂t(ubar·D)` (should match to machine precision)
4. **WEC term check**: Verify `Σ_k idiss ≈ ∫J + frc2dx + numerical_diss` (consistency)
5. **Wind stress check**: Verify `sustr_output - brk2dx_output` matches a separate pure-wind run (if available)
6. **r_D sanity**: Verify r_D values are O(0.001-0.01) m/s in deep water, larger near coast with waves

## CHANGELOG Entry

```
### cppdefs.opt
- Added `#define DIAGNOSTICS` — enables 7-term momentum diagnostic budget.

### diagnostics.opt (NEW override)
- Set `diag_uv = .true.`, `output_period = 900`, `diag_avg = .true.`.

### diagnostics.F (NEW override)
- Added depth-integrated diagnostic arrays `Udiag_int(nx,ny,7)`,
  `Vdiag_int(nx,ny,7)`. Computed as Σ_k Udiag(i,j,k,n) in new subroutine
  `calc_diag_int`. Written to diagnostics NetCDF as `u_pgr_int`, `u_cor_int`,
  `u_adv_int`, `u_dis_int`, `u_hmx_int`, `u_vmx_int`, `u_cpl_int` (and v_*).
  Units: m²/s² (depth-integrated du/dt contribution per term).

### wec_frc.F
- Added `hvf_int_xi`, `hvf_int_eta` — depth-integrated horizontal vortex force.
  Accumulated in `compute_wec_rhs_uv_terms` from the Stokes-drift × velocity
  shear terms (cff1). Reset at start of each call. Written to WEC output.
- Added `r_D` output (bottom drag coefficient from shared memory).

### wec_frc.opt
- Added `wrt_hvf_int = .true.`, `wrt_rD = .true.`.
```
