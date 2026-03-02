# WEC CPP Options Tree (Code-Verified)

All flags verified against `DEBUG_L4PTSAL/wec_frc.F` (override) and `src/wec_frc.F`.

## Master Flag

```
WEC                          Master WEC switch
```

Everything below is inside `#ifdef WEC`.

---

## 1. Wave Forcing Source (mutually exclusive)

```
WAVE_OFFLINE                 Read wave fields from file (WW3 or other wave model)
│
├── MONOWAVE                 Stokes drift formulation choice
│   ├── defined  → Monochromatic (U2010 Eq. 13-14: cosh 2k(z+h)/sinh 2kD)
│   │              Computes ust2d, vst2d from stk*wdrx/wdre
│   │              Uses kw for pressure terms (calP) and Akw
│   │
│   └── undef    → Broadband two-scale (R2021 Eq. 15-17)
│       │          Reads ust2d/vst2d from WW3 (ustr2d/vstr2d)
│       │          Computes keff iteratively
│       │          Normalizes by fac2d
│       │          Uses keff for pressure terms (calP)
│       │          Allocates: keff, fac2d, inv_exef
│       │
│       └── WAVE_DIFF_FROM_LM   Mean wavelength input (broadband only)
│           ├── defined  → Read mean wavelength (llmm) from WW3 for Akw wavenumber
│           │              keffkw = 2*pi/llmm (wec_frc.F:2069, inside #ifndef MONOWAVE)
│           └── undef    → Use online computation
│           Only meaningful with !MONOWAVE — monochromatic uses kw directly
│           **CURRENTLY UNDEF'd**: Akw is hardcoded to 0.0 (line 2080), so the
│           entire llmm → keffkw → Akw pipeline is dead code. No need to read
│           the extra wwlm field from WW3. Can re-enable if Akw is restored.
│           Affects: wec_frc.F (6 locations)
│
│   Affects: wec_frc.F (15+ locations), prsgrd.F, compute_prsgrd_wec_terms.h
│
├── SUP_OFF                  Wave setup input source
│   ├── defined  → Read sup (setup) from WW3 file
│   └── undef    → Compute sup online from wave fields
│   Affects: wec_frc.F (9 locations), compute_prsgrd_wec_terms.h
│
├── UORBOFF                  Orbital velocity input source
│   ├── defined  → Read uorb/vorb from WW3 file (more accurate)
│   └── undef    → Compute Uorb_mag from Soulsby formula: 0.5*fr*wh/sinh(kD)
│   Either way → Uorb_mag feeds tauw_orb → r_D (BBL_S95 implicit drag)
│   Affects: wec_frc.F (10 locations)
│
└── SURFACE_ROLLER           Surface roller model
    ├── defined  → Enable roller energy pathway (additional breaking cascade)
    │              Reads roller data from file (WAVE_OFFLINE+SURFACE_ROLLER)
    │              Adds roller contributions to momentum and mixing
    └── undef    → No roller
    Affects: wec_frc.F (7 locations in override, 12+ in src)

ANA_WEC_FRC                  Analytical wave forcing (alternative to WAVE_OFFLINE)
                             Mutually exclusive with WAVE_OFFLINE in practice
                             Affects: wec_frc.F (6 locations)

WKB_WWAVE                   WKB online wave model (alternative to WAVE_OFFLINE)
                             Mostly commented out in current code
```

---

## 2. Wave Breaking

```
BREAKING                     Master switch for ALL breaking effects (momentum + mixing)
│                            When undef'd: ebrk=0, brk=0, Akb block skipped → zero breaking
│                            All sub-flags below only meaningful when BREAKING is defined
│
├── 2a. Breaking Momentum
│   │
│   ├── SURFACE_BREAK        Enables breaking momentum injection
│   │   ├── BODY_BREAK (sub-flag, controls representation)
│   │   │   ├── undef (default) → Surface stress: brk2dx added to sustr (top cell)
│   │   │   │                      Guard: BREAKING && SURFACE_BREAK && !BODY_BREAK
│   │   │   └── defined          → 3D body force: brk3dx(k) in ru(k), brk3de(k) in rv(k)
│   │   │                          Guard: BREAKING && SURFACE_BREAK && BODY_BREAK
│   │   │                          Allocates: brk_r, brk3dx, brk3de (3D arrays)
│   │   │
│   │   │   Design mirrors WAVE_FRICTION / BODY_FRICTION pattern:
│   │   │   SURFACE_BREAK alone      → stress (default)
│   │   │   SURFACE_BREAK+BODY_BREAK → body force
│   │   │   Neither SURFACE_BREAK    → no breaking momentum (even with BREAKING)
│   │   │
│   │   └── Affects: wec_frc.F (7 locations: declarations, allocations, body force
│   │                computation, ru/rv addition, sustr contamination, MPI exchange)
│   │
│   ├── BRK0                 Debug: zeros brk entirely (momentum only, NOT Akb)
│   │   Guard: BREAKING && BRK0
│   │   defined  → brk = 0 (no breaking momentum at all)
│   │
│   ├── SMOOTH_BRK           Laplacian smoothing of brk(i,j)
│   │   Guard: BREAKING && SMOOTH_BRK && !BRK0
│   │   Spreads brk isotropically via N_brk_smooth=4 Laplacian passes.
│   │   Represents unresolved cross-shore extent of surfzone forcing.
│   │   Applied BEFORE SHELTER_BRK. Only brk is smoothed, not ebrk (Akb unaffected).
│   │   Affects: wec_frc.F (FX_brk/FE_brk/FE1_brk scratch, smooth_brk.h include)
│   │
│   └── SHELTER_BRK          Coastline shelter: scale brk by embayment concavity
│       Guard: BREAKING && SHELTER_BRK && !BRK0
│       Precomputed shelter_fac from rmask via smoothed land_influence.
│       Corrects for WW3 coastline being coarser than ROMS (unresolved sheltering).
│       Applied AFTER SMOOTH_BRK: brk = brk * shelter_fac.
│       Parameters: N_shelter=8, Lf0=0.5 (hardcoded in compute_shelter_fac.h)
│       Affects: wec_frc.F (shelter_fac array, compute_shelter_fac.h include)
│
├── 2b. Breaking-Induced Mixing (Akb)
│   │   Guard: BREAKING && LMD_MIXING && !AKB0
│   │
│   ├── Akb_U2010            Akb recipe selection
│   │   ├── defined  → U2010 simple scaling: Akb = c_b * (ebrk^(1/3)) * wh * F(z)
│   │   │              Key parameter: bconst = c_b = 0.03 (calibrated, was 0.1)
│   │   └── undef   → Romero 2021 TKE-based (the #else branch)
│   │                  FKV_FUNC4 defined in cppdefs.opt but NEVER REFERENCED (GHOST)
│   │
│   ├── AKB0                 Debug: zeros Akb entirely
│   │   Guard: skips Akb block when defined (in the master #if)
│   │
│   └── FKV_FUNC4            GHOST — defined but zero references
│
└── Source: ebrk/brk computed from wveb (combined S_ds + S_db from WW3)
    Guard: #ifdef BREAKING (else ebrk=0, brk=0)
```

---

## 3. Wave-Current Bottom Interaction

### 3a. Bottom Streaming (WAVE_FRICTION)

```
WAVE_FRICTION                Master switch for ALL streaming effects (momentum + BBL Cr(z))
│                            When undef'd: frc=0, frc2dx/frc2de=0, Cr(z) term skipped
│                            All sub-flags below only meaningful when WAVE_FRICTION is defined
│
├── FRC0                    Debug: zeros frc entirely (momentum only, NOT Cr(z))
│                            Guard: WAVE_FRICTION && FRC0 (wec_frc.F, after BRK0 block)
│                            undef (default) → frc computed from wved
│                            defined          → frc=0 → frc2dx/frc2de=0
│                            Mirrors BRK0 for breaking.
│
├── Momentum pathway:
│   ├── BODY_FRICTION (sub-flag, controls vertical distribution)
│   │   ├── undef (default) → Bottom cell only: frc2dx added to ru(1)
│   │   │                      In pre_step3d4S.F:283 and step3d_uv1.F:141
│   │   │                      WBBL-unresolved approximation (U2010 §3.3)
│   │   └── defined          → Distributed body force: frc3dx(k) in ru(k) for all k
│   │                          In wec_frc.F:1204,1236 (compute_wec_rhs_uv_terms)
│   │                          Uses f^wd(z) vertical profile (U2010 Eq. 54-55)
│   │
│   │   Mutually exclusive #ifdef guards, no double-counting.
│   │   Affects: wec_frc.F (5 locations), pre_step3d4S.F, step3d_uv1.F
│   │
│   └── Source: frc = wved * kw * inv_f (wec_frc.F:1639)
│       Guard: #ifdef WAVE_FRICTION (else frc=0)
│       frc2dx = 0.5*(frc(i-1)*wdrx(i-1) + frc(i)*wdrx(i))
│
├── BBL detection pathway:
│   └── FRC_HBBL0            Debug: zeros streaming contribution in Cr(z) for BBL detection
│       Guard: WEC && WAVE_FRICTION && !FRC_HBBL0 (in lmd_kpp.F:269)
│       undef (default) → C_Ek*f²*|z|/H term active in Cr(z) (D. Hypolite addition)
│       defined          → streaming stabilization term zeroed
│
└── PHYSICS: B^wd = (ε^wd/ρ₀σ) k f^wd(z) — body force (U2010 Eq. 54)
    Source: wved = ε^wd (wave frictional dissipation from wave model)
    Direction: wave propagation direction (projected via wdrx/wdre)

Bottom streaming profile shapes (src/ only, commented out):
  FF_FUNC1                   Bottom streaming vertical distribution
                             Defined but commented out in src/wec_frc.F:1809
```

### 3b. Wave-Enhanced Bottom Drag (Implicit)

```
BBL_S95                      Soulsby 1995 wave-enhanced bottom drag coefficient
│   defined  → Modifies r_D drag coefficient using Uorb_mag and tauw_orb
│              r_D enters IMPLICITLY in tridiagonal solver (not additive force)
│              Computed in compute_rd_bott_drag.h:49-55 (#ifdef BBL_S95)
│
│   Two downstream effects of wave-enhanced r_D:
│   1. Momentum drag: r_D in tridiagonal solver (always active with BBL_S95)
│   2. BKPP mixing: wmb = vonKar*sqrt(r_D*|U|) → deeper/stronger bottom BL
│
│   ├── WMB0               Debug: disable wave-enhanced wmb in BKPP mixing
│   │   Guard: WEC && BBL_S95 && !WMB0 (in lmd_kpp.F:465)
│   │   undef (default) → wmb uses wave-enhanced r_D (consistent with momentum)
│   │   defined          → wmb falls back to current-only log-layer formula
│   │                      Decouples: r_D still enhances momentum drag,
│   │                      but BKPP mixing depth/intensity unaffected by waves
│   │
│   Input chain:
│   UORBOFF → Uorb_mag → tauw_orb → cffwec2 enhancement factor → r_D
│   ├── UORBOFF defined  → Uorb_mag = sqrt(uorb²+vorb²) from WW3 file
│   └── UORBOFF undef    → Uorb_mag = 0.5*fr*wh/sinh(kD) (Soulsby formula)
│
│   Affects: wec_frc.F (6 locations), lmd_kpp.F (2), compute_rd_bott_drag.h

Alternative BBL formulations (src/ only, mostly commented out):
  BBL_F00                    Feddersen 2000 bottom drag
  BBL_QUAD2D                 Quadratic 2D bottom drag

BBL                          REMOVED from cppdefs.opt
│   Was the original ROMS BBL module flag (Styles & Glenn BBL96).
│   In src/: guarded uorb/vorb declarations and BBL_S95 block.
│   Override merged all BBL orbital velocity functionality into UORBOFF.
│   Override removed BBL guard from BBL_S95 block.
│   No active references remain — removed from cppdefs.opt.
```

---

## 4. Wave-Turbulence Interaction

```
KPPLT                        Langmuir turbulence parameterization for KPP
                             GHOST FLAG — undef'd in cppdefs.opt, ZERO references
                             in any .F file (neither override nor src/)
                             Possibly planned but never implemented
```

---

## 5. Internal Flags (inside src/wec_frc.F, not in cppdefs.opt)

```
KH_SOULSBY / KH_HUNT        Wavenumber solver choice (src/:1786-1787)
                             KH_HUNT active, KH_SOULSBY undef'd

FB_WSCALE                    Breaking velocity scale (src/:1802, commented out)
KVF0                         Zero Kv override (src/:1804, commented out, no #ifdef refs)
```

---

## 6. Ghost/Dead Flags in cppdefs.opt

| Flag | Status | Notes |
|------|--------|-------|
| `FKV_FUNC4` | GHOST | Defined but never `#ifdef`'d anywhere — R2021 path is #else of Akb_U2010 |
| `KPPLT` | GHOST | Undef'd, zero code references |
| `LOG_BDRAG` | GHOST | Undef'd, zero code references |

---

## Current cppdefs.opt Configuration

```fortran
#define WAVE_OFFLINE
#define BREAKING            ← master switch for ALL breaking effects
  #define SURFACE_BREAK     ← breaking momentum as surface stress
  #undef  BODY_BREAK        ← stress (default); define for body force
  #undef  BRK0              ← breaking momentum NOT zeroed
  #undef  AKB0              ← breaking mixing NOT zeroed
  #undef  Akb_U2010         ← using Romero 2021 Akb
  #define FKV_FUNC4         ← GHOST (no effect)
#define WAVE_FRICTION       ← bottom streaming ON
  #undef  FRC0              ← streaming momentum NOT zeroed
  #undef  BODY_FRICTION     ← bottom cell (default); define for body force
  #undef  FRC_HBBL0         ← Cr(z) streaming active (default); define to zero
#define BBL_S95             ← wave-enhanced drag ON
  #undef  WMB0              ← wave-enhanced wmb in BKPP active (default); define to disable
#define UORBOFF             ← read Uorb from WW3
#define SUP_OFF             ← read setup from WW3
#undef  WAVE_DIFF_FROM_LM   ← Akw disabled, llmm input not needed
#undef  KPPLT               ← GHOST
#undef  MONOWAVE            ← broadband (R2021) Stokes
#undef  LOG_BDRAG           ← GHOST
#undef  ANA_WEC_FRC         ← not analytical
#undef  WKB_WWAVE           ← not online wave model
```

---

## Verification: Streaming Code Paths

The streaming momentum representation is **correctly implemented** (no double-counting):

| Guard | Path | Where applied |
|-------|------|---------------|
| `WAVE_FRICTION && !BODY_FRICTION` | Bottom cell body force: `frc2dx → ru(1)` | `pre_step3d4S.F:283`, `step3d_uv1.F:141` |
| `WAVE_FRICTION && BODY_FRICTION` | Distributed body force: `frc3dx(k) → ru(k)` | `wec_frc.F:1201,1233` (`compute_wec_rhs_uv_terms`) |

Both paths use `ru` (body force). Neither enters through `bustr`/`bvstr`.
The WBBL thickness (δ_w ~ few cm, U2010 Eq. 56) is unresolvable on the model grid.
The `!BODY_FRICTION` path lumps the depth-integrated body force into the bottom cell.
U2010 Table 2 Run 25 tests f^wd = "S" (Eq. 57: τ^wd_bot = ε^wd k/σ) but this stress
formulation is not implemented as a true bustr boundary stress in the code.
Mutually exclusive guards, single pathway at each stage.

---

## Streaming vs Enhanced Bottom Drag: Two Distinct Phenomena

**These are NOT the same effect.** They represent different physical mechanisms and
are not double-counting, despite both involving wave-bottom interaction.

### Streaming (WAVE_FRICTION → frc2dx)

- **What**: Wave momentum transferred to mean flow via WBBL Reynolds stress divergence
- **Paper**: U2010 Eq. 54: B^wd = (ε^wd/ρ₀σ) k f^wd(z)
- **Code**: `frc = wved * kw * inv_f` (wec_frc.F:1632)
- **Source**: `wved` = ε^wd (rate of wave energy dissipation at bottom)
- **Direction**: Wave propagation direction (projected via wdrx/wdre)
- **Effect**: Waves GIVE momentum to the mean current
- **Enters**: `ru(1)` or `ru(k)` — explicit body force

### Wave-Enhanced Bottom Drag (BBL_S95 → r_D)

- **What**: Wave orbital motion enhances effective bottom roughness for mean current
- **Paper**: U2010 Eq. 60: τ^cd_bot = τ_c [1 + 1.2(|τ_w|/(|τ_w|+|τ_c|))^3.2]
- **Code**: `tauw_orb` (wec_frc.F:1777) → `r_D` (compute_rd_bott_drag.h:49-54)
- **Source**: `Uorb_mag` (wave orbital velocity, from file or Soulsby formula)
- **Direction**: Opposing the mean current direction
- **Effect**: Bottom TAKES more momentum from the mean current
- **Enters**: `r_D` coefficient in tridiagonal solver — implicit

### Why they are independent

| | Streaming (B^wd) | Enhanced drag (τ^cd_bot) |
|---|---|---|
| **Energy source** | Wave energy (ε^wd) | Mean current kinetic energy |
| **Direction of force** | Wave propagation direction | Opposing current direction |
| **Code input** | `wved` (dissipation rate) | `Uorb_mag` (orbital velocity) |
| **Code output** | `frc2dx → ru(1)` | `tauw_orb → r_D → tridiag solver` |
| **Paper equation** | Eq. 54 | Eq. 60 |
| **Would exist without mean current?** | Yes (streaming occurs even in still water) | No (drag requires a current to oppose) |

Streaming transfers wave momentum to the current. Enhanced drag increases the rate
at which the bottom dissipates current momentum. These can coexist and partly
oppose each other: streaming pushes flow in the wave direction, enhanced drag
decelerates whatever current exists.

---

## Breaking Type Distinction: None (Verified)

**Romero 2021 Section 2.3 (p.4):**
> "We do not distinguish between mixing due to white-capping in deep water
> and depth-induced breaking. Our approach is suitable for both deep and
> shallow wave breaking, giving depth-averaged diffusivities"

**Code verification (wec_frc.F):**
The code makes NO distinction between depth-induced and deep-water (whitecapping) breaking.
A single input field `wweb` is read from WW3 (line 1278), which is the combined total
breaking dissipation (WW3 sums S_ds + S_db internally before export).

| Stage | Code | Uses |
|-------|------|------|
| Input | `wveb = wweb * rmask` (line 1461) | Single combined field |
| ε_b | `ebrk = wramp2 * wveb` (line 1626) | Single scalar per grid point |
| Momentum | `brk = ebrk * kw * inv_f` (line 1634) → `brk2dx` (line 1715) | Combined ebrk |
| Mixing (U2010) | `vs = ebrk^(1/3)` (line 2113) → Akb profile | Combined ebrk |
| Mixing (R2021) | `Akb = ebrk/hs * kb0(k)` (line 2156) | Combined ebrk |

No CPP flags, no conditional logic, and no separate arrays exist for different breaking types.
The only difference between the two Akb recipes (U2010 vs R2021) is the vertical profile
shape — both use the same combined ebrk as input.

**NOTE**: Comment at line 2091 says "Eddy viscosity due to **depth-induced** wave breaking"
which is misleading — the code handles ALL breaking types (deep-water + depth-induced)
through this same ebrk pathway. Consider correcting to "Eddy viscosity due to wave breaking".

---

## Summary: 3 Non-Conservative WEC Effects

| Effect | Energy source | Momentum pathway | Mixing pathway |
|--------|-------------|------------------|----------------|
| **Breaking** | ε^b (`ebrk`) | stress (`sustr`) OR body force (`ru(k)`) | Akb → Akv + Akt |
| **Streaming** | ε^wd (`wved`) | body force only: `ru(1)` or `ru(k)` | BBL Cr(z) stabilization |
| **Bottom drag** | BBL_S95 | implicit coefficient `r_D` | Enhanced wmb → deeper BKPP (WMB0 to disable) |

---

## Proposed Clean-Up Summary

1. **Remove ghost flags**: `FKV_FUNC4`, `KPPLT`, `LOG_BDRAG`
2. **Remove BBL from cppdefs.opt**: footgun with no WEC purpose (override already ignores it)
3. **Fix breaking momentum**: Replace fragile `#ifndef SURFACE_BREAK` with explicit
   `BODYFORCE_BREAK` (see `BREAKING_FIX_PLAN.md`)
4. **Keep UORBOFF**: Valid input source switch (WW3 file vs Soulsby formula for Uorb_mag)
5. **Keep internal src/ flags** (`KH_HUNT`, `FKV_FUNC2`, `FB_FUNC3`): stable recipe choices
   that don't need user-facing CPP keys
