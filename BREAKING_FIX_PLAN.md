# Plan: Fix Wave Breaking Double-Counting in DEBUG_L4PTSAL

## Current State: Not Actively Broken

The current `DEBUG_L4PTSAL` config (`SURFACE_BREAK` defined, `BRK0` undef) does **not**
double-count. The override developer made a coordinated set of changes:
1. Added `sustr += brk2dx` unconditionally in `wec_frc.F` (line 1724)
2. Removed the explicit `brk2dx` block from the `pre_step3d4S.F` override
3. Defined `SURFACE_BREAK` to disable the body force path (`brk3dx`)

Result: the **only** breaking pathway is via contaminated `sustr`, delivering exactly 1×
at every stage (predictor, corrector, barotropic). This works but is **fragile** —
double-counting would reappear if `SURFACE_BREAK` is undef'ed (body force + sustr = 2×)
or if the source `pre_step3d4S.F` is used (explicit + sustr = 2×).

## Goal

Replace the fragile sustr-contamination approach with a clean, explicit, CPP-guarded design
that is safe under all flag combinations. See `WAVE_BREAKING_ANALYSIS.md` §4 for full analysis.

## Design: Three Modes via CPP Flags

| Mode | CPP flags | Behavior |
|------|-----------|----------|
| **No breaking** (default) | neither `SURFACE_BREAK` nor `BODYFORCE_BREAK` | No breaking momentum — safe default, equivalent to `BRK0` |
| **Surface stress** | `#define SURFACE_BREAK` | `brk2dx` added explicitly to top cell in predictor + corrector + barotropic |
| **Body force** | `#define BODYFORCE_BREAK` | `brk3dx(k)` added to `ru(k)` via `compute_wec_rhs_uv_terms` |

`BRK0` can be removed (or kept as a separate safety switch for debugging).
Defining both `SURFACE_BREAK` and `BODYFORCE_BREAK` simultaneously must be prohibited.

## Files to Modify (all overrides, never src/)

### 1. `wec_frc.F` — Remove sustr contamination, guard body force with new key

**a) Remove unconditional `sustr += brk2dx` (lines 1721–1724, 1746–1749):**
```fortran
! REMOVE these lines:
          !ADD brk2dx to sustr ...
          sustr(i,j) = sustr(i,j) + brk2dx(i,j)    ! line 1724
          ...
          svstr(i,j) = svstr(i,j) + brk2de(i,j)    ! line 1749
```

**b) Replace `#ifndef SURFACE_BREAK` with `#ifdef BODYFORCE_BREAK` at all locations:**

| Location | Lines | What |
|----------|-------|------|
| Array declarations | 131–134 | `brk_r`, `brk3dx`, `brk3de` allocation guards |
| Array allocations | 819–825 | Same arrays |
| Body force computation | 2200–2251 | `brk_r(k) = fb * brk`, `brk3dx`, `brk3de` |
| Addition to `ru(k)` | 1197–1200 | `ru += cff * dm_u * brk3dx` |
| Addition to `rv(k)` | 1229–1232 | `rv += cff * dn_v * brk3de` |
| MPI exchange | 2275–2277 | `call exchange_xxx(brk3dx,brk3de)` |

Change: `# ifndef SURFACE_BREAK` → `# ifdef BODYFORCE_BREAK`

**c) Add compile-time guard against defining both:**
```fortran
# if defined SURFACE_BREAK && defined BODYFORCE_BREAK
        ERROR: Cannot define both SURFACE_BREAK and BODYFORCE_BREAK
# endif
```

### 2. `pre_step3d4S.F` — Add explicit SURFACE_BREAK brk2dx to predictor

The current override removed the `#ifdef SURFACE_BREAK` block. Restore it so that
when `SURFACE_BREAK` is defined, breaking enters the predictor explicitly (sustr is now clean):

**u-component (after line 406, inside top cell):**
```fortran
          DC(i,N)=cff*( DC(i,N)
     &                                  +dtau*sustr(i,j)
#    if defined WEC && defined SURFACE_BREAK
     &                            + dtau*brk2dx(i,j)
#    endif
     &                                                 )
```

**v-component (after line 469, inside top cell):**
```fortran
            DC(i,N)=cff*( DC(i,N)
     &                                    +dtau*svstr(i,j)
#    if defined WEC && defined SURFACE_BREAK
     &                            + dtau*brk2de(i,j)
#    endif
     &                                                   )
```

### 3. `step3d_uv1.F` — Add explicit SURFACE_BREAK brk2dx to corrector

The corrector needs `brk2dx`/`brk2de` in both the top cell and barotropic accumulator.
Changes required in **both** the `IMPLICIT_BOTTOM_DRAG` path (active, lines 158–315)
and the `#else` path (currently dead code, lines 317–378).

**IMPLICIT_BOTTOM_DRAG path (ACTIVE):**

u top cell (line 178–180):
```fortran
          DC(i,N)=cff*( u(i,j,N,nnew) +DC(i,0)*ru(i,j,N)
     &                                     +dt*sustr(i,j)
#if defined WEC && defined SURFACE_BREAK
     &                                     +dt*brk2dx(i,j)
#endif
     &                                                  )
```

u barotropic (lines 212–214):
```fortran
          rufrc(i,j)=ru(i,j,1) +dm_u(i,j)*dn_u(i,j)*(
     &                                      sustr(i,j)
#if defined WEC && defined SURFACE_BREAK
     &                                     +brk2dx(i,j)
#endif
     &            -0.5*(r_D(i-1,j)+r_D(i,j))*DC(i,1) )
```

v top cell (lines 255–257):
```fortran
            DC(i,N)=cff*( v(i,j,N,nnew) +DC(i,0)*rv(i,j,N)
     &                                       +dt*svstr(i,j)
#if defined WEC && defined SURFACE_BREAK
     &                                       +dt*brk2de(i,j)
#endif
     &                                                    )
```

v barotropic (lines 287–289):
```fortran
            rvfrc(i,j)=rv(i,j,1) +dm_v(i,j)*dn_v(i,j)*(
     &                                        svstr(i,j)
#if defined WEC && defined SURFACE_BREAK
     &                                       +brk2de(i,j)
#endif
     &              -0.5*(r_D(i,j-1)+r_D(i,j))*DC(i,1) )
```

**Non-IMPLICIT path (dead code, lines 317–378) — same pattern:**

u top cell (line 322–323): add `+dt*brk2dx`
u rufrc (lines 326–328): add `+brk2dx`
v top cell (line 352–353): add `+dt*brk2de`
v rvfrc (lines 356–359): add `+brk2de`
(Same `#if defined WEC && defined SURFACE_BREAK` guards)

### 4. `cppdefs.opt` — Add BODYFORCE_BREAK, clean up

```fortran
#    define SURFACE_BREAK      /* Breaking as surface stress */
#    undef  BODYFORCE_BREAK    /* Breaking as 3D body force (mutually exclusive with SURFACE_BREAK) */
                               /* Neither defined = no breaking momentum (safe default) */
#    undef  BRK0               /* Legacy: zeros brk — redundant if neither BREAK flag defined */
```

## Verification Matrix

After the fix, each mode should produce exactly 1× breaking at every stage:

| Config | Predictor | Corrector | Barotropic |
|--------|-----------|-----------|------------|
| `SURFACE_BREAK` | `sustr + brk2dx` = 1× | `sustr + brk2dx` = 1× | `sustr + brk2dx` = 1× |
| `BODYFORCE_BREAK` | via `ru(k)` = 1× | via `ru(k)` = 1× | via `Σru(k)` = 1× |
| Neither (default) | 0× | 0× | 0× |

For `BODYFORCE_BREAK`: `compute_wec_rhs_uv_terms` is called in both predictor
(`pre_step3d4S.F:172`) and corrector (`step3d_uv1.F:121`), adding `brk3dx` to `ru(k)`.
Barotropic `rufrc` accumulates `Σru(k)` which naturally includes the body force.
No changes needed to predictor or corrector for this path — it works by construction
(same pattern as friction).

## Files Touched Summary

| File | Already an override? | Action |
|------|---------------------|--------|
| `wec_frc.F` | Yes | Remove sustr contamination; `#ifndef SURFACE_BREAK` → `#ifdef BODYFORCE_BREAK` |
| `pre_step3d4S.F` | Yes | Add `#ifdef SURFACE_BREAK` brk2dx to predictor top cell |
| `step3d_uv1.F` | Yes | Add `#ifdef SURFACE_BREAK` brk2dx to corrector top cell + barotropic |
| `cppdefs.opt` | Yes | Add `BODYFORCE_BREAK`, update comments |

No files need to be copied from `src/` — all required overrides already exist.

---

## Scan for Other Dissipation Double-Counting

Checked all wave dissipation effects for the same class of bug (force entering through
two simultaneous pathways at the same time-stepping stage):

| Effect | Input | Pathways | Double-counted? |
|--------|-------|----------|-----------------|
| **Breaking momentum** | `wveb` → `ebrk` → `brk` | `sustr` contamination + body force `ru(k)` | **Yes** (latent — safe only because SURFACE_BREAK defined and pre_step3d4S override) |
| **Breaking mixing (Akb)** | `wveb` → `ebrk` | `Akv = Kv` then `+= Akb` in `lmd_kpp.F:515-525` | **No** — Akv is reset to Kv before each addition. Old buggy averaging (`0.5*Akv + 0.5*Kv`, lines 513-514) correctly commented out. |
| **Friction (frc2dx)** | `wved` → `frc` | `ru(1)` bottom cell only. NEVER enters `bustr`/`sustr`. | **No** — one pathway, applied at both stages. Mutually exclusive with `BODY_FRICTION` 3D path. |
| **Enhanced bottom drag (r_D)** | `Uorb_mag`, `tauw_orb` via BBL_S95 | Coefficient in tridiagonal solver + rufrc | **No** — just enhances a single drag coefficient, not an additive force. |
| **Akw (wave diffusivity)** | — | Hardcoded to `0.0` (`wec_frc.F:2080`) | **N/A** — disabled. |
| **ebrk → both brk and Akb** | `wveb` | Momentum via `brk`, mixing via `Akb` | **No** — different physical effects (momentum transfer vs turbulent mixing) from the same energy source. Physically correct. |

**Friction is clean.** Two representations exist:
- 2D: `frc2dx → ru(1)` guarded by `WAVE_FRICTION && !BODY_FRICTION`
- 3D: `frc3dx → ru(k)` guarded by `WAVE_FRICTION && BODY_FRICTION`

These are **mutually exclusive** — unlike the old breaking design where body force and
sustr contamination were both active simultaneously. This is the correct pattern.

**Bottom drag vs streaming are independent physics:**
- `r_D` (BBL_S95): wave orbital motion enhances the drag coefficient on the Eulerian current
- `frc` (bottom streaming): orbital motion creates a net Stokes-layer streaming force

These are additive by construction (different mechanisms), not duplicate accounting.

**Conclusion: breaking momentum is the only double-counting issue.** All other dissipation
effects are correctly implemented with single-pathway, CPP-guarded designs.
