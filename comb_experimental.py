"""Top-10 most robust TMR triples — EMPIRICAL dual fault-injection data.

Companion to comb.py (which uses the analytical Model A metric). This script
builds the same top-10 bar chart, but from the empirical results of
``calculateFitnessDualClusterInjectionRedundant`` (two independent faults
injected into two distinct modules). The conditional failure probability is
P_fail = EP/100 = (100 - fitness)/100. Genome repetition is allowed, so the 165
triples span homogeneous {A,A,A}, partially degenerate {A,A,B} and fully
degenerate {A,B,C}.
"""
import re
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

genomes = ['4PG1', '4PG2', '4PG3', '4PG4', '4PG5', '4PG6', '4PG7', '4PG8', '4PG9']

RAW_RESULTS = """
Fitness for ('4PG1', '4PG1', '4PG1'): 89.01228733459358 EP: 10.99
Fitness for ('4PG1', '4PG1', '4PG2'): 89.79206049149339 EP: 10.21
Fitness for ('4PG1', '4PG1', '4PG3'): 89.06354515050168 EP: 10.94
Fitness for ('4PG1', '4PG1', '4PG4'): 87.0311167945439 EP: 12.97
Fitness for ('4PG1', '4PG1', '4PG5'): 89.32022084195998 EP: 10.68
Fitness for ('4PG1', '4PG1', '4PG6'): 86.09515996718622 EP: 13.90
Fitness for ('4PG1', '4PG1', '4PG7'): 87.23338802296965 EP: 12.77
Fitness for ('4PG1', '4PG1', '4PG8'): 87.86561264822134 EP: 12.13
Fitness for ('4PG1', '4PG1', '4PG9'): 86.57493061979649 EP: 13.43
Fitness for ('4PG1', '4PG2', '4PG2'): 89.98503465658474 EP: 10.01
Fitness for ('4PG1', '4PG2', '4PG3'): 89.46070234113712 EP: 10.54
Fitness for ('4PG1', '4PG2', '4PG4'): 87.53196930946292 EP: 12.47
Fitness for ('4PG1', '4PG2', '4PG5'): 89.78174603174604 EP: 10.22
Fitness for ('4PG1', '4PG2', '4PG6'): 86.68990976210009 EP: 13.31
Fitness for ('4PG1', '4PG2', '4PG7'): 87.68970467596391 EP: 12.31
Fitness for ('4PG1', '4PG2', '4PG8'): 88.25098814229248 EP: 11.75
Fitness for ('4PG1', '4PG2', '4PG9'): 87.11262719703979 EP: 12.89
Fitness for ('4PG1', '4PG3', '4PG3'): 88.6727078891258 EP: 11.33
Fitness for ('4PG1', '4PG3', '4PG4'): 86.65263876251137 EP: 13.35
Fitness for ('4PG1', '4PG3', '4PG5'): 89.1507703595011 EP: 10.85
Fitness for ('4PG1', '4PG3', '4PG6'): 85.77209098862643 EP: 14.23
Fitness for ('4PG1', '4PG3', '4PG7'): 86.86023622047244 EP: 13.14
Fitness for ('4PG1', '4PG3', '4PG8'): 87.53159224936815 EP: 12.47
Fitness for ('4PG1', '4PG3', '4PG9'): 86.18323442136499 EP: 13.82
Fitness for ('4PG1', '4PG4', '4PG4'): 83.75 EP: 16.25
Fitness for ('4PG1', '4PG4', '4PG5'): 87.06450094161958 EP: 12.94
Fitness for ('4PG1', '4PG4', '4PG6'): 82.74657924743444 EP: 17.25
Fitness for ('4PG1', '4PG4', '4PG7'): 84.06499429874572 EP: 15.94
Fitness for ('4PG1', '4PG4', '4PG8'): 84.92888402625822 EP: 15.07
Fitness for ('4PG1', '4PG4', '4PG9'): 83.09399477806788 EP: 16.91
Fitness for ('4PG1', '4PG5', '4PG5'): 89.3371212121212 EP: 10.66
Fitness for ('4PG1', '4PG5', '4PG6'): 86.13122171945702 EP: 13.87
Fitness for ('4PG1', '4PG5', '4PG7'): 87.24547511312217 EP: 12.75
Fitness for ('4PG1', '4PG5', '4PG8'): 87.8810975609756 EP: 12.12
Fitness for ('4PG1', '4PG5', '4PG9'): 86.59707991803278 EP: 13.40
Fitness for ('4PG1', '4PG6', '4PG6'): 81.62568306010928 EP: 18.37
Fitness for ('4PG1', '4PG6', '4PG7'): 83.06010928961749 EP: 16.94
Fitness for ('4PG1', '4PG6', '4PG8'): 83.94543546694648 EP: 16.05
Fitness for ('4PG1', '4PG6', '4PG9'): 81.99906367041199 EP: 18.00
Fitness for ('4PG1', '4PG7', '4PG7'): 84.34426229508196 EP: 15.66
Fitness for ('4PG1', '4PG7', '4PG8'): 85.1718258132214 EP: 14.83
Fitness for ('4PG1', '4PG7', '4PG9'): 83.42696629213484 EP: 16.57
Fitness for ('4PG1', '4PG8', '4PG8'): 85.88709677419355 EP: 14.11
Fitness for ('4PG1', '4PG8', '4PG9'): 84.33761961722487 EP: 15.66
Fitness for ('4PG1', '4PG9', '4PG9'): 82.32758620689656 EP: 17.67
Fitness for ('4PG2', '4PG2', '4PG2'): 89.59120982986768 EP: 10.41
Fitness for ('4PG2', '4PG2', '4PG3'): 89.23494983277592 EP: 10.77
Fitness for ('4PG2', '4PG2', '4PG4'): 87.23891730605285 EP: 12.76
Fitness for ('4PG2', '4PG2', '4PG5'): 89.60058661145618 EP: 10.40
Fitness for ('4PG2', '4PG2', '4PG6'): 86.5207136997539 EP: 13.48
Fitness for ('4PG2', '4PG2', '4PG7'): 87.38207547169812 EP: 12.62
Fitness for ('4PG2', '4PG2', '4PG8'): 87.9001976284585 EP: 12.10
Fitness for ('4PG2', '4PG2', '4PG9'): 86.78885291396855 EP: 13.21
Fitness for ('4PG2', '4PG3', '4PG3'): 88.63717128642502 EP: 11.36
Fitness for ('4PG2', '4PG3', '4PG4'): 86.60145586897178 EP: 13.40
Fitness for ('4PG2', '4PG3', '4PG5'): 89.16911225238444 EP: 10.83
Fitness for ('4PG2', '4PG3', '4PG6'): 85.84317585301837 EP: 14.16
Fitness for ('4PG2', '4PG3', '4PG7'): 86.7836832895888 EP: 13.22
Fitness for ('4PG2', '4PG3', '4PG8'): 87.39995787700084 EP: 12.60
Fitness for ('4PG2', '4PG3', '4PG9'): 86.12141444114738 EP: 13.88
Fitness for ('4PG2', '4PG4', '4PG4'): 83.67559523809524 EP: 16.32
Fitness for ('4PG2', '4PG4', '4PG5'): 87.08215630885122 EP: 12.92
Fitness for ('4PG2', '4PG4', '4PG6'): 82.83209806157355 EP: 17.17
Fitness for ('4PG2', '4PG4', '4PG7'): 83.95809578107183 EP: 16.04
Fitness for ('4PG2', '4PG4', '4PG8'): 84.75109409190372 EP: 15.25
Fitness for ('4PG2', '4PG4', '4PG9'): 83.00424281984334 EP: 17.00
Fitness for ('4PG2', '4PG5', '4PG5'): 89.4128787878788 EP: 10.59
Fitness for ('4PG2', '4PG5', '4PG6'): 86.27262443438914 EP: 13.73
Fitness for ('4PG2', '4PG5', '4PG7'): 87.2341628959276 EP: 12.77
Fitness for ('4PG2', '4PG5', '4PG8'): 87.81032229965156 EP: 12.19
Fitness for ('4PG2', '4PG5', '4PG9'): 86.60988729508196 EP: 13.39
Fitness for ('4PG2', '4PG6', '4PG6'): 81.85792349726776 EP: 18.14
Fitness for ('4PG2', '4PG6', '4PG7'): 83.10792349726775 EP: 16.89
Fitness for ('4PG2', '4PG6', '4PG8'): 83.91920251836306 EP: 16.08
Fitness for ('4PG2', '4PG6', '4PG9'): 82.0848938826467 EP: 17.92
Fitness for ('4PG2', '4PG7', '4PG7'): 84.20765027322405 EP: 15.79
Fitness for ('4PG2', '4PG7', '4PG8'): 84.9685204616999 EP: 15.03
Fitness for ('4PG2', '4PG7', '4PG9'): 83.30212234706616 EP: 16.70
Fitness for ('4PG2', '4PG8', '4PG8'): 85.62247983870968 EP: 14.38
Fitness for ('4PG2', '4PG8', '4PG9'): 84.13576555023924 EP: 15.86
Fitness for ('4PG2', '4PG9', '4PG9'): 82.21982758620689 EP: 17.78
Fitness for ('4PG3', '4PG3', '4PG3'): 87.75510204081633 EP: 12.24
Fitness for ('4PG3', '4PG3', '4PG4'): 85.61710398445092 EP: 14.38
Fitness for ('4PG3', '4PG3', '4PG5'): 88.46604215456675 EP: 11.53
Fitness for ('4PG3', '4PG3', '4PG6'): 84.82726423902895 EP: 15.17
Fitness for ('4PG3', '4PG3', '4PG7'): 85.85434173669468 EP: 14.15
Fitness for ('4PG3', '4PG3', '4PG8'): 86.59029649595688 EP: 13.41
Fitness for ('4PG3', '4PG3', '4PG9'): 85.07936507936508 EP: 14.92
Fitness for ('4PG3', '4PG4', '4PG4'): 82.58928571428571 EP: 17.41
Fitness for ('4PG3', '4PG4', '4PG5'): 86.32419517102615 EP: 13.68
Fitness for ('4PG3', '4PG4', '4PG6'): 81.7002442002442 EP: 18.30
Fitness for ('4PG3', '4PG4', '4PG7'): 82.93650793650794 EP: 17.06
Fitness for ('4PG3', '4PG4', '4PG8'): 83.8773419203747 EP: 16.12
Fitness for ('4PG3', '4PG4', '4PG9'): 81.86274509803921 EP: 18.14
Fitness for ('4PG3', '4PG5', '4PG5'): 88.88104838709677 EP: 11.12
Fitness for ('4PG3', '4PG5', '4PG6'): 85.47101449275362 EP: 14.53
Fitness for ('4PG3', '4PG5', '4PG7'): 86.52173913043478 EP: 13.48
Fitness for ('4PG3', '4PG5', '4PG8'): 87.20957249070632 EP: 12.79
Fitness for ('4PG3', '4PG5', '4PG9'): 85.8141447368421 EP: 14.19
Fitness for ('4PG3', '4PG6', '4PG6'): 80.67251461988305 EP: 19.33
Fitness for ('4PG3', '4PG6', '4PG7'): 82.03947368421053 EP: 17.96
Fitness for ('4PG3', '4PG6', '4PG8'): 82.996632996633 EP: 17.00
Fitness for ('4PG3', '4PG6', '4PG9'): 80.89022757697457 EP: 19.11
Fitness for ('4PG3', '4PG7', '4PG7'): 83.24561403508773 EP: 16.75
Fitness for ('4PG3', '4PG7', '4PG8'): 84.14702581369248 EP: 15.85
Fitness for ('4PG3', '4PG7', '4PG9'): 82.2289156626506 EP: 17.77
Fitness for ('4PG3', '4PG8', '4PG8'): 84.92726293103449 EP: 15.07
Fitness for ('4PG3', '4PG8', '4PG9'): 83.22115384615385 EP: 16.78
Fitness for ('4PG3', '4PG9', '4PG9'): 81.01851851851852 EP: 18.98
Fitness for ('4PG4', '4PG4', '4PG4'): 78.31632653061224 EP: 21.68
Fitness for ('4PG4', '4PG4', '4PG5'): 83.25066137566137 EP: 16.75
Fitness for ('4PG4', '4PG4', '4PG6'): 77.31331168831169 EP: 22.69
Fitness for ('4PG4', '4PG4', '4PG7'): 78.81493506493507 EP: 21.19
Fitness for ('4PG4', '4PG4', '4PG8'): 80.02717391304348 EP: 19.97
Fitness for ('4PG4', '4PG4', '4PG9'): 77.34962406015038 EP: 22.65
Fitness for ('4PG4', '4PG5', '4PG5'): 86.640625 EP: 13.36
Fitness for ('4PG4', '4PG5', '4PG6'): 82.2863924050633 EP: 17.71
Fitness for ('4PG4', '4PG5', '4PG7'): 83.55221518987341 EP: 16.45
Fitness for ('4PG4', '4PG5', '4PG8'): 84.43567961165049 EP: 15.56
Fitness for ('4PG4', '4PG5', '4PG9'): 82.54905523255815 EP: 17.45
Fitness for ('4PG4', '4PG6', '4PG6'): 76.12403100775194 EP: 23.88
Fitness for ('4PG4', '4PG6', '4PG7'): 77.80038759689923 EP: 22.20
Fitness for ('4PG4', '4PG6', '4PG8'): 79.03375370919882 EP: 20.97
Fitness for ('4PG4', '4PG6', '4PG9'): 76.23207885304659 EP: 23.77
Fitness for ('4PG4', '4PG7', '4PG7'): 79.26356589147287 EP: 20.74
Fitness for ('4PG4', '4PG7', '4PG8'): 80.42470326409496 EP: 19.58
Fitness for ('4PG4', '4PG7', '4PG9'): 77.8673835125448 EP: 22.13
Fitness for ('4PG4', '4PG8', '4PG8'): 81.42755681818183 EP: 18.57
Fitness for ('4PG4', '4PG8', '4PG9'): 79.14169520547945 EP: 20.86
Fitness for ('4PG4', '4PG9', '4PG9'): 76.25 EP: 23.75
Fitness for ('4PG5', '4PG5', '4PG5'): 88.96875 EP: 11.03
Fitness for ('4PG5', '4PG5', '4PG6'): 85.7375 EP: 14.26
Fitness for ('4PG5', '4PG5', '4PG7'): 86.8125 EP: 13.19
Fitness for ('4PG5', '4PG5', '4PG8'): 87.4639423076923 EP: 12.54
Fitness for ('4PG5', '4PG5', '4PG9'): 86.1221590909091 EP: 13.88
Fitness for ('4PG5', '4PG6', '4PG6'): 81.18181818181817 EP: 18.82
Fitness for ('4PG5', '4PG6', '4PG7'): 82.58333333333333 EP: 17.42
Fitness for ('4PG5', '4PG6', '4PG8'): 83.48837209302326 EP: 16.51
Fitness for ('4PG5', '4PG6', '4PG9'): 81.49305555555556 EP: 18.51
Fitness for ('4PG5', '4PG7', '4PG7'): 83.81818181818181 EP: 16.18
Fitness for ('4PG5', '4PG7', '4PG8'): 84.66569767441861 EP: 15.33
Fitness for ('4PG5', '4PG7', '4PG9'): 82.86458333333333 EP: 17.14
Fitness for ('4PG5', '4PG8', '4PG8'): 85.39341517857143 EP: 14.61
Fitness for ('4PG5', '4PG8', '4PG9'): 83.79321808510637 EP: 16.21
Fitness for ('4PG5', '4PG9', '4PG9'): 81.73076923076923 EP: 18.27
Fitness for ('4PG6', '4PG6', '4PG6'): 74.77777777777777 EP: 25.22
Fitness for ('4PG6', '4PG6', '4PG7'): 76.61111111111111 EP: 23.39
Fitness for ('4PG6', '4PG6', '4PG8'): 77.87234042553192 EP: 22.13
Fitness for ('4PG6', '4PG6', '4PG9'): 74.91452991452991 EP: 25.09
Fitness for ('4PG6', '4PG7', '4PG7'): 78.24074074074075 EP: 21.76
Fitness for ('4PG6', '4PG7', '4PG8'): 79.42375886524823 EP: 20.58
Fitness for ('4PG6', '4PG7', '4PG9'): 76.741452991453 EP: 23.26
Fitness for ('4PG6', '4PG8', '4PG8'): 80.45176630434783 EP: 19.55
Fitness for ('4PG6', '4PG8', '4PG9'): 78.04330065359477 EP: 21.96
Fitness for ('4PG6', '4PG9', '4PG9'): 75.0 EP: 25.00
Fitness for ('4PG7', '4PG7', '4PG7'): 79.66666666666666 EP: 20.33
Fitness for ('4PG7', '4PG7', '4PG8'): 80.78014184397163 EP: 19.22
Fitness for ('4PG7', '4PG7', '4PG9'): 78.33333333333333 EP: 21.67
Fitness for ('4PG7', '4PG8', '4PG8'): 81.74252717391305 EP: 18.26
Fitness for ('4PG7', '4PG8', '4PG9'): 79.5547385620915 EP: 20.45
Fitness for ('4PG7', '4PG9', '4PG9'): 76.78571428571429 EP: 23.21
Fitness for ('4PG8', '4PG8', '4PG8'): 82.568359375 EP: 17.43
Fitness for ('4PG8', '4PG8', '4PG9'): 80.60546875 EP: 19.39
Fitness for ('4PG8', '4PG9', '4PG9'): 78.125 EP: 21.88
Fitness for ('4PG9', '4PG9', '4PG9'): 75.0 EP: 25.00
"""


def regime(triple):
    a, b, c = triple
    if a == b == c:
        return 'homogeneous'
    if a == b or b == c or a == c:
        return 'partial'
    return 'distinct'


# Parse: (triple, P_fail, regime). P_fail = EP/100 = (100 - fitness)/100.
records = []
pattern = re.compile(r"\('(4PG\d)', '(4PG\d)', '(4PG\d)'\):\s*([\d.]+)\s*EP:")
for a, b, c, fit in pattern.findall(RAW_RESULTS):
    tr = (a, b, c)
    records.append((tr, (100.0 - float(fit)) / 100.0, regime(tr)))

COLOR = {'homogeneous': '#d73027', 'partial': '#f4a259', 'distinct': '#1a9850'}
LABEL = {'homogeneous': 'Homogeneous {A,A,A} (classical TMR)',
         'partial': 'Partially degenerate {A,A,B}',
         'distinct': 'Fully degenerate {A,B,C}'}

# ---- Top-10 ranking, coloured by degeneracy regime ----
TOP_N = 10
ranking = sorted(records, key=lambda r: r[1])
top = ranking[:TOP_N]
counts = {k: sum(1 for tr, v, kd in top if kd == k) for k in COLOR}

fig, ax = plt.subplots(figsize=(9, 6))
names = ['\n'.join(tr) for tr, v, kd in top]
vals = [v for tr, v, kd in top]
colors = [COLOR[kd] for tr, v, kd in top]
bars = ax.bar(range(len(names)), vals, color=colors, edgecolor='black')
ax.set_xticks(range(len(names)))
ax.set_xticklabels(names, fontsize=8)
ax.set_ylabel(r'$\bar{P}_{\mathrm{TMR}}$  — mean conditional failure probability (dual injection)')
ax.set_title(f'{TOP_N} most robust TMR configurations (empirical dual injection)\n'
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
plt.savefig('tmr_top_comparison_experimental.png', dpi=150, bbox_inches='tight')
plt.savefig('tmr_top_comparison_experimental.eps', format='eps', bbox_inches='tight')
plt.close()
print("Saved tmr_top_comparison_experimental.png / .eps")
print(f"Top-{TOP_N}: {counts['distinct']} distinct / {counts['partial']} partial / "
      f"{counts['homogeneous']} homogeneous")
for tr, v, kd in top:
    print(f"  {tr}  P_fail={v:.4f}  [{kd}]")
