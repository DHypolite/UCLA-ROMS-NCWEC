#!/usr/bin/env python
"""
compare_budgets.py — Depth-integrated momentum budget comparison (WEC vs No-WEC)

Compares 8 depth-integrated momentum diagnostic terms between WEC and
no-WEC simulations.  Produces 3 figures:
  1. Budget bar-chart comparison (spatial means)
  2. Breaking-pressure-gradient imbalance maps  (WEC only)
  3. Sub-balance / closure checks               (WEC only)
"""

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

# ══════════════════════════════ User config ══════════════════════════
grid_file = 'grid.nc'
dia_wec   = 'wec_dia.nc'
dia_nowec = 'nowec_dia.nc'
wec_file  = 'wec_wec_his.nc'       # brk2dx, frc2dx, hvf_int, r_D

imin, imax = 50, 80                # nearshore region (rho-grid indices, inclusive)
jmin, jmax = 20, 60
trec = -1                          # time record (-1 = last)

# Coordinate mode: 'grid' or 'rotated'
coord_mode = 'grid'

# Only used when coord_mode='rotated':
# shore_angle = angle of ALONG-SHORE direction in degrees, CCW from XI-AXIS
# Cross-shore = 90 deg CCW from along-shore (points offshore / left of along-shore)
shore_angle = 0.0

# ══════════════════════════════ Constants ════════════════════════════
DIA_TERMS = ['pgr', 'cor', 'adv', 'dis', 'hmx', 'vmx', 'cpl', 'wec']

TERM_LABELS = {
    'pgr': 'PGF',  'cor': 'Cor',  'adv': 'Adv',  'dis': 'Diss',
    'hmx': 'Hmix', 'vmx': 'Vmix', 'cpl': 'Cpl',  'wec': 'WEC',
}

# ══════════════════════════════ Helpers ══════════════════════════════

def u_to_rho(u, nxi_rho):
    """Interpolate u-grid (eta_rho, xi_u) -> rho-grid (eta_rho, xi_rho).

    Interior: average adjacent u-points.  Boundaries: copy nearest.
    """
    out = np.empty((u.shape[0], nxi_rho))
    out[:, 0]    = u[:, 0]
    out[:, -1]   = u[:, -1]
    out[:, 1:-1] = 0.5 * (u[:, :-1] + u[:, 1:])
    return out


def v_to_rho(v, neta_rho):
    """Interpolate v-grid (eta_v, xi_rho) -> rho-grid (eta_rho, xi_rho).

    Interior: average adjacent v-points.  Boundaries: copy nearest.
    """
    out = np.empty((neta_rho, v.shape[1]))
    out[0, :]    = v[0, :]
    out[-1, :]   = v[-1, :]
    out[1:-1, :] = 0.5 * (v[:-1, :] + v[1:, :])
    return out


def any_to_rho(data, grid_type, neta_rho, nxi_rho):
    """Interpolate from any staggered grid to rho-grid."""
    if grid_type == 'u':
        return u_to_rho(data, nxi_rho)
    elif grid_type == 'v':
        return v_to_rho(data, neta_rho)
    return data  # already rho


def detect_grid(var):
    """Return 'u', 'v', or 'rho' from a NetCDF variable's dimension names."""
    dims = var.dimensions
    if 'xi_u' in dims:
        return 'u'
    if 'eta_v' in dims:
        return 'v'
    return 'rho'


def read_2d(ds, vname, trec):
    """Read a single 2D field at time record *trec*, return (data, grid_type)."""
    var  = ds[vname]
    data = var[trec, :, :]
    return np.asarray(data), detect_grid(var)


def rotate(u_rho, v_rho, theta_deg):
    """Rotate (xi, eta) components to (cross-shore, along-shore).

    theta = shore_angle (deg, CCW from xi-axis to along-shore direction).
      along  =  u*cos(th) + v*sin(th)
      cross  = -u*sin(th) + v*cos(th)
    Returns (cross, along).
    """
    th = np.deg2rad(theta_deg)
    c, s = np.cos(th), np.sin(th)
    along =  u_rho * c + v_rho * s
    cross = -u_rho * s + v_rho * c
    return cross, along


def safe_pctile(a, pct):
    """np.nanpercentile that returns a small default if all NaN."""
    v = a[np.isfinite(a)]
    if v.size == 0:
        return 1e-10
    out = np.nanpercentile(np.abs(v), pct)
    return out if out > 0 else 1e-10


# ══════════════════════════════ Load grid ════════════════════════════
print('Loading grid ...')
grd = nc.Dataset(grid_file, 'r')
rmask = np.asarray(grd['mask_rho'][:])
h     = np.asarray(grd['h'][:])
neta_rho, nxi_rho = rmask.shape
grd.close()

# Rho-grid subset (inclusive bounds)
js, je = jmin, jmax + 1
is_, ie = imin, imax + 1

mask_sub = rmask[js:je, is_:ie]
mask_nan = np.where(mask_sub > 0.5, 1.0, np.nan)   # 1 over ocean, NaN over land

print(f'  rho-grid shape : ({neta_rho}, {nxi_rho})')
print(f'  subset         : j=[{jmin}:{jmax}], i=[{imin}:{imax}]')

# ══════════════════════════════ Load diagnostics ════════════════════

def load_dia(fname):
    """Load all 8 depth-integrated budget terms, interpolate to rho-grid,
    subset, and mask.  Returns dicts u[term], v[term]."""
    ds = nc.Dataset(fname, 'r')
    u, v = {}, {}
    for t in DIA_TERMS:
        # u-terms (xi_u, eta_rho)
        uname = f'u_{t}_int'
        udata, ugrid = read_2d(ds, uname, trec)
        u[t] = any_to_rho(udata, ugrid, neta_rho, nxi_rho)[js:je, is_:ie] * mask_nan

        # v-terms (xi_rho, eta_v)
        vname = f'v_{t}_int'
        vdata, vgrid = read_2d(ds, vname, trec)
        v[t] = any_to_rho(vdata, vgrid, neta_rho, nxi_rho)[js:je, is_:ie] * mask_nan
    ds.close()
    return u, v


print('Loading WEC diagnostics ...')
u_w, v_w = load_dia(dia_wec)

print('Loading no-WEC diagnostics ...')
u_n, v_n = load_dia(dia_nowec)

# ══════════════════════════════ Load WEC 2D fields ══════════════════
print('Loading WEC history fields ...')
wds = nc.Dataset(wec_file, 'r')

def load_wec_field(vname):
    data, grid = read_2d(wds, vname, trec)
    return any_to_rho(data, grid, neta_rho, nxi_rho)[js:je, is_:ie] * mask_nan

brk_xi  = load_wec_field('brk2dx')
brk_eta = load_wec_field('brk2de')
frc_xi  = load_wec_field('frc2dx')
frc_eta = load_wec_field('frc2de')
hvf_xi  = load_wec_field('hvf_int_xi')
hvf_eta = load_wec_field('hvf_int_eta')
r_D     = load_wec_field('r_D')

wds.close()

# ══════════════════════════════ Rotation ════════════════════════════
if coord_mode == 'rotated':
    comp_labels = ('cross-shore', 'along-shore')
    for t in DIA_TERMS:
        u_w[t], v_w[t] = rotate(u_w[t], v_w[t], shore_angle)
        u_n[t], v_n[t] = rotate(u_n[t], v_n[t], shore_angle)
    brk_xi,  brk_eta  = rotate(brk_xi,  brk_eta,  shore_angle)
    frc_xi,  frc_eta  = rotate(frc_xi,  frc_eta,  shore_angle)
    hvf_xi,  hvf_eta  = rotate(hvf_xi,  hvf_eta,  shore_angle)
else:
    comp_labels = ('xi-component', 'eta-component')

# After rotation:  u_* / *_xi = first component,  v_* / *_eta = second

# ══════════════════════════════ Summary statistics ══════════════════
def smean(f):
    return np.nanmean(f)

print('\n─── Spatial-mean budget [m²/s²] ───')
hdr = f'{"term":>6s}  {"u_WEC":>11s}  {"u_noWEC":>11s}  {"v_WEC":>11s}  {"v_noWEC":>11s}'
print(hdr)
for t in DIA_TERMS:
    print(f'{TERM_LABELS[t]:>6s}  {smean(u_w[t]):+11.4e}  {smean(u_n[t]):+11.4e}'
          f'  {smean(v_w[t]):+11.4e}  {smean(v_n[t]):+11.4e}')
sum_uw = smean(sum(u_w[t] for t in DIA_TERMS))
sum_un = smean(sum(u_n[t] for t in DIA_TERMS))
sum_vw = smean(sum(v_w[t] for t in DIA_TERMS))
sum_vn = smean(sum(v_n[t] for t in DIA_TERMS))
print(f'{"Sum":>6s}  {sum_uw:+11.4e}  {sum_un:+11.4e}  {sum_vw:+11.4e}  {sum_vn:+11.4e}')

print(f'\nbrk_xi  mean = {smean(brk_xi):+.4e}')
print(f'frc_xi  mean = {smean(frc_xi):+.4e}')
print(f'hvf_xi  mean = {smean(hvf_xi):+.4e}')
print(f'wec_int mean = {smean(u_w["wec"]):+.4e}')
print(f'hvf+frc mean = {smean(hvf_xi + frc_xi):+.4e}')

# ══════════════════════════════ Figure 1 ════════════════════════════
print('\nPlotting Figure 1: Budget bar charts ...')
fig1, axes1 = plt.subplots(1, 2, figsize=(14, 5), sharey=True)

for ax, (dw, dn), label in zip(
        axes1,
        [(u_w, u_n), (v_w, v_n)],
        comp_labels):

    means_w = [smean(dw[t]) for t in DIA_TERMS]
    means_n = [smean(dn[t]) for t in DIA_TERMS]
    means_w.append(np.nanmean(sum(dw[t] for t in DIA_TERMS)))
    means_n.append(np.nanmean(sum(dn[t] for t in DIA_TERMS)))

    xlabels = [TERM_LABELS[t] for t in DIA_TERMS] + ['Sum']
    x = np.arange(len(xlabels))
    w = 0.35

    ax.bar(x - w/2, means_w, w, label='WEC',    color='tab:blue')
    ax.bar(x + w/2, means_n, w, label='No WEC', color='tab:orange')
    ax.set_xticks(x)
    ax.set_xticklabels(xlabels, rotation=45, ha='right')
    ax.set_ylabel(r'$\overline{F}$  [m$^2$/s$^2$]')
    ax.set_title(f'Depth-integrated budget — {label}')
    ax.axhline(0, color='k', lw=0.5)
    ax.legend(fontsize=9)

fig1.tight_layout()
fig1.savefig('fig1_budget_bars.png', dpi=150)
print('  -> fig1_budget_bars.png')

# ══════════════════════════════ Figure 2 ════════════════════════════
print('Plotting Figure 2: Breaking-PGF balance ...')
fig2, axes2 = plt.subplots(2, 2, figsize=(12, 10))

# (a) brk + pgr imbalance — first component
imbal_u = brk_xi + u_w['pgr']
vlim = safe_pctile(imbal_u, 95)
im = axes2[0, 0].pcolormesh(imbal_u, cmap='RdBu_r', vmin=-vlim, vmax=vlim,
                             shading='auto')
axes2[0, 0].set_title(f'brk + pgr  ({comp_labels[0]})')
axes2[0, 0].set_aspect('equal')
fig2.colorbar(im, ax=axes2[0, 0], shrink=0.8, label=r'm$^2$/s$^2$')

# (b) brk + pgr imbalance — second component
imbal_v = brk_eta + v_w['pgr']
vlim2 = safe_pctile(imbal_v, 95)
im = axes2[0, 1].pcolormesh(imbal_v, cmap='RdBu_r', vmin=-vlim2, vmax=vlim2,
                             shading='auto')
axes2[0, 1].set_title(f'brk + pgr  ({comp_labels[1]})')
axes2[0, 1].set_aspect('equal')
fig2.colorbar(im, ax=axes2[0, 1], shrink=0.8, label=r'm$^2$/s$^2$')

# (c) brk vs -pgr scatter (first component)
b_flat = brk_xi.ravel()
p_flat = u_w['pgr'].ravel()
ok = np.isfinite(b_flat) & np.isfinite(p_flat)
if ok.any():
    axes2[1, 0].scatter(b_flat[ok], -p_flat[ok], s=2, alpha=0.4)
    lim = max(np.abs(b_flat[ok]).max(), np.abs(p_flat[ok]).max()) * 1.1
    if lim > 0:
        axes2[1, 0].plot([-lim, lim], [-lim, lim], 'k--', lw=0.8, label='1:1')
        axes2[1, 0].set_xlim(-lim, lim); axes2[1, 0].set_ylim(-lim, lim)
axes2[1, 0].set_xlabel(f'brk  ({comp_labels[0]})')
axes2[1, 0].set_ylabel(f'-pgr  ({comp_labels[0]})')
axes2[1, 0].set_title('brk  vs  -pgr  (should be 1:1)')
axes2[1, 0].set_aspect('equal')
axes2[1, 0].legend(fontsize=9)

# (d) fractional imbalance |brk+pgr| / |brk|
with np.errstate(divide='ignore', invalid='ignore'):
    frac = np.abs(imbal_u) / np.abs(brk_xi)
    frac[np.abs(brk_xi) < 1e-12] = np.nan
im = axes2[1, 1].pcolormesh(frac, cmap='hot_r', vmin=0, vmax=2, shading='auto')
axes2[1, 1].set_title(f'|brk+pgr| / |brk|  ({comp_labels[0]})')
axes2[1, 1].set_aspect('equal')
fig2.colorbar(im, ax=axes2[1, 1], shrink=0.8)

fig2.tight_layout()
fig2.savefig('fig2_breaking_balance.png', dpi=150)
print('  -> fig2_breaking_balance.png')

# ══════════════════════════════ Figure 3 ════════════════════════════
print('Plotting Figure 3: Sub-balance checks ...')
fig3, axes3 = plt.subplots(2, 2, figsize=(12, 10))

# (a) u_wec vs hvf + frc scatter (first component)
wec_flat = u_w['wec'].ravel()
sub_flat = (hvf_xi + frc_xi).ravel()
ok = np.isfinite(wec_flat) & np.isfinite(sub_flat)
if ok.any():
    axes3[0, 0].scatter(sub_flat[ok], wec_flat[ok], s=2, alpha=0.4)
    lim = max(np.abs(wec_flat[ok]).max(), np.abs(sub_flat[ok]).max()) * 1.1
    if lim > 0:
        axes3[0, 0].plot([-lim, lim], [-lim, lim], 'k--', lw=0.8, label='1:1')
        axes3[0, 0].set_xlim(-lim, lim); axes3[0, 0].set_ylim(-lim, lim)
axes3[0, 0].set_xlabel(f'hvf + frc  ({comp_labels[0]})')
axes3[0, 0].set_ylabel(f'wec_int  ({comp_labels[0]})')
axes3[0, 0].set_title(r'$u_\mathrm{wec}$  vs  hvf$_\xi$ + frc$_\xi$')
axes3[0, 0].set_aspect('equal')
axes3[0, 0].legend(fontsize=9)

# (b) v_wec vs hvf + frc scatter (second component)
wec_flat = v_w['wec'].ravel()
sub_flat = (hvf_eta + frc_eta).ravel()
ok = np.isfinite(wec_flat) & np.isfinite(sub_flat)
if ok.any():
    axes3[0, 1].scatter(sub_flat[ok], wec_flat[ok], s=2, alpha=0.4)
    lim = max(np.abs(wec_flat[ok]).max(), np.abs(sub_flat[ok]).max()) * 1.1
    if lim > 0:
        axes3[0, 1].plot([-lim, lim], [-lim, lim], 'k--', lw=0.8, label='1:1')
        axes3[0, 1].set_xlim(-lim, lim); axes3[0, 1].set_ylim(-lim, lim)
axes3[0, 1].set_xlabel(f'hvf + frc  ({comp_labels[1]})')
axes3[0, 1].set_ylabel(f'wec_int  ({comp_labels[1]})')
axes3[0, 1].set_title(r'$v_\mathrm{wec}$  vs  hvf$_\eta$ + frc$_\eta$')
axes3[0, 1].set_aspect('equal')
axes3[0, 1].legend(fontsize=9)

# (c) Sum of all 8 terms — spatial map (= tendency, should be small)
tend_u = sum(u_w[t] for t in DIA_TERMS)
vlim = safe_pctile(tend_u, 95)
im = axes3[1, 0].pcolormesh(tend_u, cmap='RdBu_r', vmin=-vlim, vmax=vlim,
                             shading='auto')
axes3[1, 0].set_title(f'Sum of all terms  ({comp_labels[0]})')
axes3[1, 0].set_aspect('equal')
fig3.colorbar(im, ax=axes3[1, 0], shrink=0.8, label=r'm$^2$/s$^2$')

# (d) Numerical dissipation — spatial map (should be small)
vlim = safe_pctile(u_w['dis'], 95)
im = axes3[1, 1].pcolormesh(u_w['dis'], cmap='RdBu_r', vmin=-vlim, vmax=vlim,
                             shading='auto')
axes3[1, 1].set_title(f'Numerical dissipation  ({comp_labels[0]})')
axes3[1, 1].set_aspect('equal')
fig3.colorbar(im, ax=axes3[1, 1], shrink=0.8, label=r'm$^2$/s$^2$')

fig3.tight_layout()
fig3.savefig('fig3_sub_balance.png', dpi=150)
print('  -> fig3_sub_balance.png')

print('\nDone.')
plt.show()
