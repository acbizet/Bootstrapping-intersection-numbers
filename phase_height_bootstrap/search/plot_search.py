"""Redraw the finite local-search diagnostic from stored results."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT=Path(__file__).resolve().parent
r=json.loads((OUT/'results.json').read_text());cases=r['cases']
fig,axs=plt.subplots(1,2,figsize=(11.5,4.1),gridspec_kw={'width_ratios':[1.5,1]},constrained_layout=True)
near=[v for v in cases if v['Q']<5900]
for marker,key,label,col in [('o','unconstrained_height_margin','Local peak minus resonant height','#2c6e9f'),('s','phase_margin','Same-phase local maximum minus resonant height','#a14580')]:
 ys=[v['unconstrained_height_margin'] if key!='phase_margin' else v['winning_phase_candidate']['height_margin'] for v in near]
 axs[0].plot([v['Q'] for v in near],ys,marker=marker,color=col,label=label)
axs[0].set_ylim(-80,20);axs[0].axhline(0,color='black',lw=.9);axs[0].set_xlabel(r'$Q$ at $P=10^7$, $n=26$');axs[0].set_ylabel('Height allowance');axs[0].set_title('Phase removes local height allowance');axs[0].legend(fontsize=8,loc='lower left')
for v in near:
 axs[0].annotate(f"{v['target_phase']:+.3f}",(v['Q'],v['winning_phase_candidate']['height_margin']),xytext=(0,-14),textcoords='offset points',ha='center',fontsize=8)
axs[0].text(.98,.98,'Labels: target phase (radians)',transform=axs[0].transAxes,ha='right',va='top',fontsize=8)
v=cases[-1];peak=v['unconstrained_height_margin'];cand=v['winning_phase_candidate'];ys=[peak,cand['height_margin']]
axs[1].bar([0,1],ys,color=['#2c6e9f','#a14580'],width=.65)
for j,y in enumerate(ys):axs[1].text(j,y+5,f'{y:.3f}',ha='center',fontsize=10)
axs[1].set_xticks([0,1],['Local peak','Same-phase\nlocal maximum']);axs[1].set_ylabel('Height allowance');axs[1].set_title(r'$Q=5940$: phase permits a point');axs[1].set_ylim(0,max(ys)*1.15);axs[1].axhline(0,color='black',lw=.8)
for ax in axs:ax.grid(axis='y',alpha=.2);ax.set_axisbelow(True)
fig.suptitle('Finite multistart optimization on the original coefficient torus\nThese local maxima are not certified global contour bounds',fontsize=11)
for ext in ['png','pdf']:fig.savefig(OUT/f'phase_height_search.{ext}',dpi=180)
