"""Plot matched-control roughness estimates from the archived follow-up analysis.

This script reads existing results and produces artwork only. Ratio intervals in
panel B are the archived Monte Carlo propagated intervals. Panel A displays the
archived standard errors times 1.96; g=0 is the same mean divided by itself.
"""
from pathlib import Path
import json
import numpy as np
from orx_figstyle import use_style, figure_grid, panel_labels, save, WIDE, PALETTE

HERE = Path(__file__).resolve().parent
ana = json.loads((HERE / 'modal_hmc_followup_analysis.json').read_text())
use_style()
fig, (a,b) = figure_grid(1,2,width=WIDE,ratio=.47,gridspec_kw={'width_ratios':[1,1.12]})
colors=[PALETTE[k] for k in ('blue','orange','green','purple')]
for n,col,m in zip([8,16,24,32],colors,['o','s','^','D']):
    rr=sorted((r for r in ana['crosscheck_gN_grid'] if r['N']==n),key=lambda r:r['g_N'])
    g=np.array([r['g_N'] for r in rr])
    y=np.array([r['smoothing_vs_measured_g0_pct'] for r in rr])
    err=np.array([r['smoothing_vs_measured_g0_se_pct'] for r in rr])*1.96
    err[g==0]=0  # correlated identity, not a second independent observation
    a.errorbar(g,y,yerr=err,fmt=m+'-',color=col,ms=4,lw=1.1,capsize=2,label=f'{n}')
a.set(xlabel=r'Coupling $g_N$',ylabel='Roughness reduction (%)',xlim=(-.08,2.12),ylim=(-1.5,53),xticks=[0,.5,1,1.5,2],yticks=[0,10,20,30,40,50])
a.legend(title=r'Matrix rank $N$',title_fontsize=8,ncol=2,loc='lower right',fontsize=8)
a.set_title('Measured zero-coupling reference',loc='left',x=.08,pad=4)
rows=ana['primary_4rank_table']
vals=[r['smoothing_vs_measured_g0'] for r in rows]
vals += [ana['lattice_12x8_test']['smoothing_12x8_vs_measured_g0'],ana['extended_n32_smoothing_crosscheck']]
labels=[r'$N=8$, $8\times4$',r'$N=16$, $8\times4$',r'$N=24$, $8\times4$',r'$N=32$, $12\times4$',r'$N=24$, $12\times8$',r'$N=32$, extended']
for i,(v,label,col) in enumerate(zip(vals,labels,colors+[PALETTE['black'],PALETTE['red']])):
    mean=v['smoothing_pct']
    b.errorbar(mean,i,xerr=[[mean-v['ci95_mc_lo']],[v['ci95_mc_hi']-mean]],fmt='o' if i<4 else 'D',color=col,ms=4.5,capsize=3,lw=1.1)
b.axhline(3.5,color='.8',lw=.6)
b.set(yticks=range(6),yticklabels=labels,xlim=(47.45,48.2),xticks=[47.5,47.7,47.9,48.1],xlabel='Roughness reduction at $g_N=2$ (%)',ylim=(5.6,-.6))
b.set_title('Size and volume checks',loc='left',x=.08,pad=4)
panel_labels([a,b])
save(fig,str(HERE.parent/'figs/fig_followup_roughness'))
