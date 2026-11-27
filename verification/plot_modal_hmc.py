"""Observable-specific follow-up diagnostics and finite-resolution benchmarks.

Uses only the existing production summaries and follow-up analysis. No simulation
is launched and no new numerical result file is written.
"""
from pathlib import Path
import json
from orx_figstyle import use_style, figure_grid, panel_labels, save, WIDE, PALETTE

HERE=Path(__file__).resolve().parent
ana=json.loads((HERE/'modal_hmc_followup_analysis.json').read_text())
runs=json.loads((HERE/'modal_hmc_final_results.json').read_text())['runs']
by_tag={r['tag']:r for r in runs}
use_style()
fig,axes=figure_grid(2,2,width=WIDE,ratio=.91)
a,b,c,d=axes.ravel()
blue=PALETTE['blue'];orange=PALETTE['orange'];black=PALETTE['black']
rr=[ana[k] for k in ['extended_n32_batch1','extended_n32_batch2','extended_n48_batch1','extended_n48_batch2']]
rr.append(ana['lattice_12x8_test']['lat_12x8_g2'])
labels=[r'$N=32$, global',r'$N=32$, +local',r'$N=48$, global',r'$N=48$, +local',r'$N=24$, $12\times8$']
for field,dy,col,m,label in [('Q2',-.11,blue,'o','Roughness'),('rho2_eig_mean',.11,orange,'s','Local area')]:
    a.scatter([r['diagnostics'][field]['rhat_rank'] for r in rr],[i+dy for i in range(5)],c=col,marker=m,s=22,label=label,zorder=3)
a.axvline(1.01,color='.45',ls='--',lw=.8)
a.set(yticks=range(5),yticklabels=labels,ylim=(4.6,-.6),xlim=(.99,1.365),xticks=[1,1.1,1.2,1.3],xlabel=r'Rank-normalized split-$\hat R$')
a.legend(loc='upper right',fontsize=8)
a.set_title('Agreement between chains',loc='left',x=.08,pad=4)
for field,dy,col,m,label in [('ess_bulk',-.11,blue,'o','Bulk'),('ess_tail',.11,orange,'s','Tail')]:
    b.scatter([r['diagnostics']['rho2_eig_mean'][field] for r in rr],[i+dy for i in range(5)],c=col,marker=m,s=22,label=label,zorder=3)
b.axvline(400,color='.45',ls='--',lw=.8)
b.set_xscale('log')
b.set(yticks=range(5),yticklabels=labels,ylim=(4.6,-.6),xlim=(10,16000),xticks=[10,100,1000,10000],xlabel='Local-area effective sample size')
b.legend(loc='lower right',fontsize=8)
b.set_title('Information after autocorrelation',loc='left',x=.08,pad=4)
for n in [8,12,16,24,32,48]:
    if n==48:
        vals=ana['extended_n48_batch2']['diagnostics'];col=blue;mfc=blue
    else:
        vals=by_tag[f'final_IR_mu005_N{n}']['observables'];col='.5';mfc='white'
    lo,hi=vals['rho2_eig_min']['mean'],vals['rho2_eig_max']['mean']
    c.plot([n,n],[lo,hi],color=col,lw=1.1)
    c.plot([n,n],[lo,hi],'o',mfc=mfc,mec=col,ms=4)
    c.plot(n,vals['rho2_eig_mean']['mean'],'D',color=col,ms=3.5)
c.axhline(1,color='.6',ls=':',lw=.7)
c.set(xlabel=r'Matrix rank $N$',ylabel='Normalized local area eigenvalue',xticks=[8,16,24,32,48],ylim=(0,2.9),xlim=(4,52))
c.text(.39,.92,'Blue: extended run\nGrey: archived samples',transform=c.transAxes,fontsize=8,va='top')
c.set_title('Mean extrema at one lattice site',loc='left',x=.08,pad=4)
for tag,col,m,lab in [('final_IR_mu005_N48',blue,'o',r'$N=48$, scaled coupling'),('final_UV_gN2.00_N32',orange,'^',r'$N=32$, $g_N=2$')]:
    mm=by_tag[tag]['spin2_modes']
    d.errorbar([q['j'] for q in mm],[q['Z_ratio_0'] for q in mm],yerr=[q['G_meas_0_se']/q['G_free_0'] for q in mm],fmt=m+'-',color=col,ms=3,lw=1,elinewidth=.6,label=lab)
d.axhline(1,color='.5',ls=':',lw=.7)
d.set(xlabel=r'Angular momentum $j$',ylabel=r'$G(j)/G_{\mathrm{free}}(j)$',xticks=[2,10,20,30,40,47],xlim=(1,49),ylim=(.94,1.04))
d.legend(loc='lower right',fontsize=8)
d.set_title('Finite-resolution propagators',loc='left',x=.08,pad=4)
panel_labels(axes)
save(fig,str(HERE.parent/'figs/fig6_modal_b200_hmc'))
