"""Plot the exact global reduction and the resulting height exclusions."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from global_height import Profile

HERE=Path(__file__).resolve().parent
rows=json.loads((HERE/'global_height_results.json').read_text())['rows']
fig,axes=plt.subplots(1,2,figsize=(11,4.2),constrained_layout=True)
colors=['#4b6cb7','#287c8e','#d17c16','#8a579c','#47834a']
row=rows[2]
p=Profile(row['source_omega'],row['source_eta'],4096)
u=np.geomspace(.00015,.027,500)
d=np.array([p.deficit(np.sin(v/2)**2) for v in u])
axes[0].plot(u,d,color='black',lw=2,label=r'$D(u)$, $Q=5883$')
for rr,c in zip(rows[:4],colors[:4]):
    axes[0].axhline(rr['budget_B0_minus_hres'],c=c,lw=1.4,
                    label=fr'$M$, $Q={int(rr["Q"])}$')
axes[0].scatter([row['minimizing_abs_u']],
                [p.deficit(row['minimizing_s'])],c='black',s=25,zorder=4)
axes[0].set(xscale='log',xlabel=r'$|u|$',ylabel=r'height deficit $D(u)$',
            ylim=(0,35),title='Whole torus reduced to one real variable')
axes[0].legend(fontsize=8,loc='upper left',ncol=2)
gaps=[r['maximum_height_minus_hres_upper'] for r in rows]
axes[1].axhline(0,c='black',lw=1)
axes[1].bar([str(int(r['Q'])) for r in rows],gaps,color=colors)
axes[1].set(yscale='symlog',xlabel=r'$Q$ at $P=10^7$, $n=26$',
            ylabel='Global maximum height minus resonant height',
            title='A negative value excludes an intersection')
axes[1].set_ylim(-25,600)
for i,g in enumerate(gaps):
    axes[1].annotate(f'{g:.6g}',(i,g),xytext=(0,5 if g>0 else -13),
                     textcoords='offset points',ha='center',fontsize=9)
for ax in axes:
    ax.grid(alpha=.2)
fig.savefig(HERE/'global_height_bootstrap.png',dpi=180)
fig.savefig(HERE/'global_height_bootstrap.pdf')
plt.close(fig)
