"""TMR analysis figures.

This script produces the figures used in the paper:

  * Figure 1 (``tmr_analysis.png``) — per-input breakdown of the conditional
    failure probability for the target TMR {4PG1, 4PG2, 4PG3}, computed
    analytically with the conditional model (Model A):

        P_TMR(L) = N(L) / (d_i*d_j + d_i*d_k + d_j*d_k)
        N(L)     = d_i*d_j*P_i*P_j + d_i*d_k*P_i*P_k + d_j*d_k*P_j*P_k

  * Figure 2a (``tmr_top_comparison.png``) — bar chart comparing the best
    heterogeneous TMR configurations against all homogeneous ones.
  * Figure 2b (``tmr_distribution.png``) — distribution of the failure
    probability over all triples.

Figures 2a/2b are built from the *numerical* (analytical) Model A metric,
evaluated for every triple directly from the per-input truth-table fault data
(the same data as the heatmap), so no empirical fault injection is required.
The combinatorial sweep now allows genome repetition, so the 165 triples span
three regimes: homogeneous {A,A,A}, partially degenerate {A,A,B} and fully
degenerate {A,B,C}.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch

# ---- data: per-input single-fault census, RECONCILED with the dual-injection
# experiment. D[g] is the number of active gates of genome g; N[L][i] is the
# number of active gates whose individual stuck-at fault PRESERVES the output of
# genome i for input L. These counts were recomputed directly from the evolved
# genotypes (PG_collection) using the same fault model as
# calculateFitnessDualClusterInjectionRedundant, so the Model A metric below
# matches the empirical dual-injection EP exactly (the previous heatmap counts
# were off by one preserved gate per cell, which lowered Model A systematically).
genomes = ['4PG1', '4PG2', '4PG3', '4PG4', '4PG5', '4PG6', '4PG7', '4PG8', '4PG9']
D = {'4PG1': 23, '4PG2': 23, '4PG3': 21, '4PG4': 14, '4PG5': 20,
     '4PG6': 15, '4PG7': 15, '4PG8': 16, '4PG9': 12}
lines = ['0000', '0001', '0010', '0011', '0100', '0101', '0110', '0111',
         '1000', '1001', '1010', '1011', '1100', '1101', '1110', '1111']
N = {
    '0000': [15, 14, 13, 8, 12, 7, 8, 8, 6],   '0001': [16, 15, 14, 8, 13, 7, 9, 10, 6],
    '0010': [15, 16, 14, 8, 14, 7, 9, 10, 6],  '0011': [17, 16, 14, 8, 14, 7, 8, 9, 6],
    '0100': [12, 16, 12, 7, 13, 7, 8, 10, 6],  '0101': [18, 16, 14, 7, 14, 9, 8, 9, 6],
    '0110': [14, 16, 15, 7, 14, 7, 8, 8, 6],   '0111': [21, 16, 12, 7, 16, 9, 8, 10, 6],
    '1000': [14, 17, 12, 8, 11, 7, 9, 10, 6],  '1001': [17, 13, 14, 8, 13, 8, 8, 8, 6],
    '1010': [15, 16, 15, 8, 12, 7, 8, 10, 6],  '1011': [15, 15, 12, 8, 16, 8, 9, 10, 6],
    '1100': [13, 16, 15, 7, 12, 7, 8, 10, 6],  '1101': [19, 15, 14, 7, 13, 8, 8, 10, 6],
    '1110': [13, 19, 14, 7, 15, 7, 8, 10, 6],  '1111': [18, 15, 16, 7, 14, 8, 8, 8, 6],
}


def pf(L, i):
    """Single-module failure probability of genome ``i`` for input ``L``."""
    return (D[genomes[i]] - N[L][i]) / D[genomes[i]]


# ----- target TMR: {4PG1, 4PG2, 4PG3} ----
i, j, k = 0, 1, 2
gi, gj, gk = genomes[i], genomes[j], genomes[k]
di, dj, dk = D[gi], D[gj], D[gk]
W_cond = di * dj + di * dk + dj * dk          # Model A denominator

scenarios_cols = [f'2 faults ∈ {gi}', f'2 faults ∈ {gj}', f'2 faults ∈ {gk}',
                  f'{gi}+{gj}',  f'{gi}+{gk}',  f'{gj}+{gk}',
                  'Max pair',     'Mean of 3 pairs',
                  r'Weighted cond. (Model A)  ★']

mat = np.zeros((len(lines), len(scenarios_cols)))
for li, L in enumerate(lines):
    p1, p2, p3 = pf(L, i), pf(L, j), pf(L, k)
    s12, s13, s23 = p1 * p2, p1 * p3, p2 * p3
    upper = max(s12, s13, s23)
    pair_avg = (s12 + s13 + s23) / 3
    cond_weighted = (di * dj * s12 + di * dk * s13 + dj * dk * s23) / W_cond   # Model A
    mat[li] = [0, 0, 0, s12, s13, s23, upper, pair_avg, cond_weighted]

# ==================================================================
# FIGURE 1 — per-input analysis for TMR {4PG1, 4PG2, 4PG3}
# ==================================================================
fig = plt.figure(figsize=(18, 11))
gs = GridSpec(2, 2, height_ratios=[1, 0.55], width_ratios=[2.2, 1],
              hspace=0.32, wspace=0.18, figure=fig)

ax1 = fig.add_subplot(gs[0, :])
cmap = mcolors.LinearSegmentedColormap.from_list(
    'fail', ['#1a9850', '#fee08b', '#d73027'])
im = ax1.imshow(mat, aspect='auto', cmap=cmap, vmin=0, vmax=mat.max())
ax1.set_xticks(range(len(scenarios_cols)))
ax1.set_xticklabels(scenarios_cols, rotation=22, ha='right')
ax1.set_yticks(range(len(lines)))
ax1.set_yticklabels(lines)
ax1.set_ylabel('Truth-table input')
ax1.set_title(f'Failure probability of TMR{{{gi}, {gj}, {gk}}} per input and scenario\n'
              'Conditional model (Model A) — matches the experimental metric',
              fontsize=13)
# separators: 3 scenarios (=0) | 3 pairs | 3 aggregations
for x in [2.5, 5.5, 6.5, 7.5]:
    ax1.axvline(x, color='white', linewidth=2)
# highlight the recommended column
ax1.add_patch(plt.Rectangle((7.5, -0.5), 1, len(lines),
                            fill=False, edgecolor='blue', linewidth=2.5))
for li in range(len(lines)):
    for ci in range(len(scenarios_cols)):
        v = mat[li, ci]
        color = 'white' if v > 0.10 else 'black'
        weight = 'bold' if ci == 8 else 'normal'
        ax1.text(ci, li, f'{v:.3f}', ha='center', va='center',
                 fontsize=8, color=color, fontweight=weight)
plt.colorbar(im, ax=ax1, fraction=0.025, pad=0.01,
             label='Failure probability')

# curve: per-input aggregations
ax2 = fig.add_subplot(gs[1, 0])
ax2.plot(lines, mat[:, 6], 'o--', label='Max pair (pessimistic bound)',
         color='#d73027', linewidth=1.5, alpha=0.7)
ax2.plot(lines, mat[:, 7], 's--', label='Mean of 3 pairs (assumes equal d)',
         color='#fc8d59', linewidth=1.5, alpha=0.7)
ax2.plot(lines, mat[:, 8], '^-', label='Weighted cond. — Model A  ★',
         color='#1a9850', linewidth=2.5)
ax2.set_xlabel('Truth-table input')
ax2.set_ylabel(r'$P(\mathrm{TMR\ fails}\mid 2\ \mathrm{faults\ in\ distinct\ modules})$')
ax2.set_title('Comparison of the three per-input aggregations')
ax2.legend(loc='upper right', fontsize=9)
ax2.grid(True, alpha=0.3)
plt.setp(ax2.get_xticklabels(), rotation=45, ha='right', fontsize=8)
ax2.set_xticks(range(len(lines)))
ax2.set_xticklabels(lines)

# bars: global aggregations (mean over the 16 inputs)
ax3 = fig.add_subplot(gs[1, 1])
labels_g = ['Max pair', 'Mean of 3 pairs', 'Weighted cond.\n(Model A) ★']
vals_g = [mat[:, 6].mean(), mat[:, 7].mean(), mat[:, 8].mean()]
colors_g = ['#d73027', '#fc8d59', '#1a9850']
edges = ['black', 'black', 'blue']
widths = [1, 1, 2.5]
bars = ax3.bar(labels_g, vals_g, color=colors_g,
               edgecolor=edges, linewidth=widths)
for bar, v in zip(bars, vals_g):
    ax3.text(bar.get_x() + bar.get_width() / 2, v + 0.003,
             f'{v:.4f}\n({100 * v:.2f}%)', ha='center', fontsize=10)
ax3.set_ylabel(r'$\bar{P}(\mathrm{TMR\ fails})$  mean over 16 inputs')
ax3.set_title(f'TMR{{{gi}, {gj}, {gk}}} — global aggregate')
ax3.set_ylim(0, max(vals_g) * 1.30)
ax3.grid(True, axis='y', alpha=0.3)

plt.savefig('tmr_analysis.png', dpi=150, bbox_inches='tight')
plt.savefig('tmr_analysis.eps', format='eps', bbox_inches='tight')
plt.close()
print("Saved tmr_analysis.png / tmr_analysis.eps")

# ==================================================================
# NUMERICAL COMBINATORIAL ANALYSIS (Model A) over ALL triples
# Repetition of genomes is now allowed, so the sweep uses
# combinations_with_replacement and spans three regimes:
#   homogeneous {A,A,A}, partially degenerate {A,A,B}, fully degenerate {A,B,C}.
# The conditional failure probability of each triple is computed analytically
# from the per-input single-module failure probabilities pf(L, i), i.e. from the
# same truth-table fault data as the heatmap -- no empirical injection required.
# ==================================================================
from itertools import combinations_with_replacement


def cond_A(triple):
    """Model A conditional failure probability of a triple (a, b, c), averaged
    over the 16 truth-table inputs. Valid for any triple, including repeated
    genomes; for a homogeneous triple it reduces to E_L[pf(L, a)^2]."""
    a, b, c = triple
    dA, dB, dC = D[genomes[a]], D[genomes[b]], D[genomes[c]]
    W = dA * dB + dA * dC + dB * dC
    vals = [(dA * dB * pf(L, a) * pf(L, b)
             + dA * dC * pf(L, a) * pf(L, c)
             + dB * dC * pf(L, b) * pf(L, c)) / W for L in lines]
    return float(np.mean(vals))


def regime(triple):
    a, b, c = triple
    if a == b == c:
        return 'homogeneous'
    if a == b or b == c or a == c:
        return 'partial'
    return 'distinct'


def tri_name(triple, sep='\n'):
    return sep.join(genomes[g] for g in triple)


# Full combinatorial sweep WITH repetition: 165 triples = C(9+3-1, 3).
records = [(tr, cond_A(tr), regime(tr))
           for tr in combinations_with_replacement(range(9), 3)]

COLOR = {'homogeneous': '#d73027', 'partial': '#f4a259', 'distinct': '#1a9850'}
LABEL = {'homogeneous': 'Homogeneous {A,A,A} (classical TMR)',
         'partial': 'Partially degenerate {A,A,B}',
         'distinct': 'Fully degenerate {A,B,C}'}

# ==================================================================
# FIGURE 2a — overall ranking of the TOP_N most robust configurations.
# Every triple competes in a single ranking, coloured by degeneracy regime,
# so the impact of admitting the partially degenerate {A,A,B} triples is visible.
# ==================================================================
TOP_N = 10
ranking = sorted(records, key=lambda r: r[1])
top = ranking[:TOP_N]
counts = {k: sum(1 for tr, v, kd in top if kd == k) for k in COLOR}

fig2a, ax = plt.subplots(figsize=(9, 6))
names = [tri_name(tr) for tr, v, kd in top]
vals = [v for tr, v, kd in top]
colors = [COLOR[kd] for tr, v, kd in top]
bars = ax.bar(range(len(names)), vals, color=colors, edgecolor='black')
ax.set_xticks(range(len(names)))
ax.set_xticklabels(names, fontsize=8)
ax.set_ylabel(r'$\bar{P}^{(A)}_{\mathrm{TMR}}$  — mean conditional failure probability (Model A)')
ax.set_title(f'{TOP_N} most robust TMR configurations (numerical Model A)\n'
             f"{counts['distinct']} fully degenerate, {counts['partial']} partially "
             f"degenerate, {counts['homogeneous']} homogeneous")
ax.grid(True, axis='y', alpha=0.3)
ax.set_ylim(0, max(vals) * 1.35)
for bar, v in zip(bars, vals):
    ax.text(bar.get_x() + bar.get_width() / 2, v + max(vals) * 0.01, f'{v:.4f}',
            ha='center', fontsize=7)
ax.legend(handles=[Patch(color=COLOR[k], label=LABEL[k])
                   for k in ['distinct', 'partial', 'homogeneous']], loc='upper right')

plt.tight_layout()
plt.savefig('tmr_top_comparison.png', dpi=150, bbox_inches='tight')
plt.savefig('tmr_top_comparison.eps', format='eps', bbox_inches='tight')
plt.close()
print("Saved tmr_top_comparison.png / tmr_top_comparison.eps")

# ==================================================================
# FIGURE 2b — distribution of all 165 triples, stacked by regime.
# ==================================================================
vals_by_kind = {k: [v for tr, v, kd in records if kd == k] for k in COLOR}

fig2b, ax = plt.subplots(figsize=(8, 6))
ax.hist([vals_by_kind['distinct'], vals_by_kind['partial'], vals_by_kind['homogeneous']],
        bins=18, stacked=True,
        color=[COLOR['distinct'], COLOR['partial'], COLOR['homogeneous']],
        edgecolor='black',
        label=[f"{LABEL['distinct']} (n={len(vals_by_kind['distinct'])})",
               f"{LABEL['partial']} (n={len(vals_by_kind['partial'])})",
               f"{LABEL['homogeneous']} (n={len(vals_by_kind['homogeneous'])})"])
best = ranking[0]
ax.axvline(best[1], color='blue', linewidth=2,
           label=f'best overall = {best[1]:.4f} ({tri_name(best[0], "+")})')
ax.set_xlabel(r'$\bar{P}^{(A)}_{\mathrm{TMR}}$  — conditional failure probability (Model A)')
ax.set_ylabel('Number of triples')
ax.set_title('Distribution of the 165 TMR triples by degeneracy regime')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('tmr_distribution.png', dpi=150, bbox_inches='tight')
plt.savefig('tmr_distribution.eps', format='eps', bbox_inches='tight')
plt.close()
print("Saved tmr_distribution.png / tmr_distribution.eps")

# ==================================================================
# Updated CSV (Model A) for the target TMR {4PG1, 4PG2, 4PG3}
# ==================================================================
import csv
with open('tmr_PG1_PG2_PG3.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['# TMR{4PG1, 4PG2, 4PG3} -- conditional model'])
    w.writerow(['# Denominator of the "Weighted cond." column:',
                f'{di * dj} + {di * dk} + {dj * dk} = {W_cond}'])
    w.writerow([])
    w.writerow(['Input'] + scenarios_cols)
    for li, L in enumerate(lines):
        w.writerow([L] + [f'{v:.5f}' for v in mat[li]])
    w.writerow([])
    w.writerow(['Mean over 16 inputs'] +
               [f'{mat[:, c].mean():.5f}' for c in range(len(scenarios_cols))])
print("Saved tmr_PG1_PG2_PG3.csv")

# ==================================================================
# Console summary -- effect of admitting partially degenerate triples
# ==================================================================
best_distinct = min((r for r in records if r[2] == 'distinct'), key=lambda r: r[1])
best_partial = min((r for r in records if r[2] == 'partial'), key=lambda r: r[1])
best_hom = min((r for r in records if r[2] == 'homogeneous'), key=lambda r: r[1])
worst = max(records, key=lambda r: r[1])

print()
print("=" * 72)
print("Numerical Model A sweep -- 165 triples (combinations_with_replacement)")
print("=" * 72)
print(f"  Best overall            : {tri_name(best[0], ',')}  ->  {best[1]:.5f}  [{best[2]}]")
print(f"  Best fully degenerate   : {tri_name(best_distinct[0], ',')}  ->  {best_distinct[1]:.5f}")
print(f"  Best partially degener. : {tri_name(best_partial[0], ',')}  ->  {best_partial[1]:.5f}")
print(f"  Best homogeneous        : {tri_name(best_hom[0], ',')}      ->  {best_hom[1]:.5f}")
print(f"  Worst overall           : {tri_name(worst[0], ',')}      ->  {worst[1]:.5f}  [{worst[2]}]")
print(f"  Top-{TOP_N} composition     : "
      f"{counts['distinct']} distinct / {counts['partial']} partial / {counts['homogeneous']} homogeneous")
gain = (best_hom[1] - best[1]) / best_hom[1] * 100
print(f"  Gain of best overall vs. best homogeneous: {gain:.2f}%")
delta = (best_partial[1] - best_distinct[1]) / best_distinct[1] * 100
print(f"  Best partial vs. best distinct: {delta:+.2f}% "
      f"(does admitting {{A,A,B}} change the optimum?)")
