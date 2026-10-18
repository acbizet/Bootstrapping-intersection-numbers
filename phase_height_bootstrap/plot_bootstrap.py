"""Overview of global height and phase exclusions; no flow illustration."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent


def main():
    d=json.loads((HERE/'summary.json').read_text());rows=d['rows']
    fig,axs=plt.subplots(1,2,figsize=(12.5,4.6))
    ix=np.arange(len(rows));values=[r['global_height_margin_upper'] for r in rows]
    colors=['#b54a4a' if v<0 else '#4277ab' for v in values]
    axs[0].bar(ix,values,color=colors,width=.65)
    axs[0].axhline(0,color='black',lw=.8)
    axs[0].set_xticks(ix,[str(r['Q']) for r in rows])
    axs[0].set_yscale('symlog',linthresh=2)
    axs[0].set_ylabel(r'Global contour maximum $-\,\mathrm{Re}\,S_{\rm res}$')
    axs[0].set_xlabel(r'Physical charge $Q$, at $P=10^7$, resonance $n=26$')
    axs[0].set_title('The sharpened height bound excludes two targets')
    for x,v in zip(ix,values):
        axs[0].annotate(f'{v:+.3f}',(x,v),xytext=(0,5 if v>0 else -13),
                        textcoords='offset points',ha='center',fontsize=9)
    axs[0].set_ylim(-35,850)
    r=next(r for r in rows if r['Q']==5883)
    bound=r['global_phase_absolute_bound'];phi=r['saddle_phase']
    axs[1].axvspan(-bound,bound,color='#69a681',alpha=.28,
                  label='All contour phases at sufficient height')
    axs[1].axvline(phi,color='#b54a4a',lw=2,label='Required saddle phase')
    axs[1].plot([phi],[0],marker='o',color='#b54a4a')
    axs[1].text(0,.34,rf'$|\Psi|\leq {bound:.6f}$',ha='center',fontsize=12)
    axs[1].text(phi,-.25,rf'$\Psi_{{\rm res}}={phi:.6f}$',ha='right',fontsize=11)
    axs[1].set_xlim(-np.pi,np.pi);axs[1].set_ylim(-.7,.7);axs[1].set_yticks([])
    axs[1].set_xticks([-np.pi,-2,-1,0,1,2,np.pi],
                     [r'$-\pi$','−2','−1','0','1','2',r'$\pi$'])
    axs[1].set_xlabel('Action phase in radians, modulo $2\pi$')
    axs[1].set_title(r'$Q=5883$: a global phase gap excludes attachment')
    axs[1].legend(loc='lower left',fontsize=8)
    fig.suptitle('Finite-charge contour bootstrap: analytical global bounds, numerical inputs',fontsize=12)
    fig.tight_layout()
    fig.savefig(HERE/'bootstrap_overview.png',dpi=190)
    fig.savefig(HERE/'bootstrap_overview.pdf')
    plt.close(fig)


if __name__=='__main__':main()
