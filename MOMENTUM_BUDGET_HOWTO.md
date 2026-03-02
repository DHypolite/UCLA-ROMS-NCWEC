# Momentum Budget Diagnostics — What to Turn On

## Override files required (in `DEBUG_L4PTSAL/`)

| File | Status | Purpose |
|------|--------|---------|
| `cppdefs.opt` | existing override | CPP flags |
| `diagnostics.F` | **new override** | depth-integrated arrays, iwec index, NetCDF output |
| `diagnostics.opt` | **new override** | enable momentum diagnostics, set output period |
| `prsgrd.F` | existing override | ipgr capture moved after WEC block |
| `compute_horiz_rhs_uv_terms.h` | **new override** | Stokes-Coriolis in icori |
| `step3d_uv1.F` | existing override | iwec capture, idiss modified |
| `wec_frc.F` | existing override | hvf_int, r_D output |
| `wec_frc.opt` | existing override | output flags for WEC 2D fields |

## Options to turn on

### cppdefs.opt
```fortran
#define DIAGNOSTICS        ! <-- REQUIRED: enables 8-term momentum budget
```
Already set. Everything else (WEC, BREAKING, WAVE_FRICTION, etc.) unchanged.

### diagnostics.opt
```fortran
diag_uv      = .true.      ! <-- REQUIRED: enables momentum diagnostics
wrt_diag_3d  = .false.     ! Write 3D per-level budget (16 vars × nz, large files)
wrt_diag_int = .true.      ! Write 2D depth-integrated budget (16 vars, compact)
diag_avg     = .true.      ! time-averaged output (set .false. for snapshots)
output_period = 900        ! averaging period in seconds (match wec_frc avg)
nrpf = 24                 ! records per file
```

### wec_frc.opt — 2D fields for offline budget decomposition
```fortran
wrt_brk2dx  = .true.      ! breaking stress xi   (= ∫B^b dz)
wrt_brk2de  = .true.      ! breaking stress eta
wrt_frc2dx  = .true.      ! streaming stress xi  (= ∫B^wd dz)
wrt_frc2de  = .true.      ! streaming stress eta
wrt_hvf_int = .true.      ! depth-int vortex force (= ∫J dz)
wrt_rD      = .true.      ! bottom drag coeff r_D [m/s]
```
All already set to `.true.`.

## Output files and variables

### Diagnostics file (`*_dia.nc`)
8-term 3D budget (units: m^2/s^2 = du/dt × dz):
- `u_pgr`, `v_pgr` — pressure gradient (incl. WEC: kvf, sup, calP, Kapsrf)
- `u_cor`, `v_cor` — Coriolis + Stokes-Coriolis
- `u_adv`, `v_adv` — horizontal advection (centered)
- `u_wec`, `v_wec` — **WEC explicit RHS** (vortex force + body forces + streaming)
- `u_dis`, `v_dis` — numerical dissipation (small residual)
- `u_hmx`, `v_hmx` — horizontal viscous mixing
- `u_vmx`, `v_vmx` — implicit vertical solver (wind + breaking + bottom drag + Akv)
- `u_cpl`, `v_cpl` — 2D/3D coupling correction

Depth-integrated versions (`*_int` suffix): same 8 terms summed over k.

### WEC file (`*_wec_his.nc` or `*_wec_avg.nc`)
- `brk2dx`, `brk2de` — breaking momentum stress [m^2/s^2]
- `frc2dx`, `frc2de` — streaming momentum stress [m^2/s^2]
- `hvf_int_xi`, `hvf_int_eta` — depth-integrated vortex force [m^2/s^2]
- `r_D` — bottom drag coefficient [m/s]

### History file (`*_his.nc`)
- `sustr`, `svstr` — total surface stress (wind + breaking) [m^2/s^2]

## Depth-integrated budget equation

```
∂t(ubar·H) = Σ ipgr + Σ icori + Σ iadv + Σ iwec + Σ idiss + Σ ihmix + Σ ivmix + Σ icoup
```

### Mapping to physical terms

```
∫∇⊥ϕ dz   = Σ ipgr          (total PGF incl. WEC pressure corrections)
∫D⊥ dz    = Σ (icori + iadv + ihmix)  (Coriolis + advection + h-mixing)
∫J dz     = hvf_int          (from WEC output, cross-check: ≈ Σ iwec - frc2dx)
∫B^b dz   = brk2dx           (from WEC output)
∫B^wd dz  = frc2dx           (from WEC output)
τ_wind    = sustr - brk2dx   (offline subtraction)
τ_bot     = -r_D · u_bot     (offline from r_D and bottom velocity)
num_diss  = Σ idiss           (should be small)
```

### Closure checks
1. `Σ (all 8 diag terms) ≈ ∂t(ubar·H)` — should match to machine precision
2. `Σ iwec ≈ hvf_int + frc2dx + body_forces` — WEC cross-check
3. `Σ ivmix ≈ sustr + τ_bot` — vertical diffusion is conservative (sums to 0)
4. `Σ idiss ≈ 0` — should be small (pure numerical dissipation)
