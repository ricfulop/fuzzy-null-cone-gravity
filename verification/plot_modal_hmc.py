"""Regenerate the reviewed HMC figure from the completed production ensembles."""
from pathlib import Path
import json
import numpy as np
from orx_figstyle import use_style, figure_grid, panel_labels, save, WIDE, PALETTE

HERE=Path(__file__).resolve().parent
runs=json.loads((HERE/'modal_hmc_final_results.json').read_text())['runs']
by_tag={r['tag']:r for r in runs}
use_style()
fig,axes=figure_grid(2,2,width=WIDE,ratio=.78)
a,b,c,d=axes.ravel()
colors=[PALETTE[k] for k in ('blue','orange','purple','black')]
summary={'roughness':[], 'area_diagnostics':[]}
# Spin-2 propagator. These errors are archived unblocked Monte Carlo estimates.
for tag,col,m,lab in [('final_IR_mu005_N48',colors[0],'o',r'$N=48$, scaled coupling'),
                      ('final_UV_gN2.00_N32',colors[1],'^',r'$N=32$, $g_N=2$')]:
    modes=by_tag[tag]['spin2_modes']
    a.errorbar([q['j'] for q in modes],[q['Z_ratio_0'] for q in modes],
        yerr=[q['G_meas_0_se']/q['G_free_0'] for q in modes],fmt=m+'-',color=col,ms=3,lw=1,elinewidth=.6,label=lab)
a.axhline(1,color='0.5',ls=':',lw=.6)
a.set(xlabel=r'Angular momentum $j$',ylabel=r'$G(j)/G_{\mathrm{free}}(j)$',xlim=(1,49),xticks=[2,10,20,30,40,47],ylim=(.94,1.04))
a.legend(loc='lower right',fontsize=7)
# Compare measured coupled ensembles at the same rank and lattice.
for n,col,m in zip([8,16,24,32],colors,['o','s','^','D']):
    base=by_tag[f'final_UV_gN0.50_N{n}']['observables']
    q0=base['Q2_anisotropy']['mean'];se0=base['Q2_anisotropy']['se']
    yy=[];ee=[]
    for g in [.5,1,1.5,2]:
        r=by_tag[f'final_UV_gN{g:.2f}_N{n}'];q=r['observables']['Q2_anisotropy'];free=r['observables']['Q2_free']
        ratio=q['mean']/q0
        se=0 if g==.5 else ratio*np.hypot(q['se']/q['mean'],se0/q0)
        yy.append(100*(1-ratio));ee.append(100*se)
        if g==2:
            summary['roughness'].append(dict(N=n,ratio_to_free=q['mean']/free,se_ratio_to_free=q['se']/free,
                reduction_from_g05=yy[-1],se_reduction_from_g05=ee[-1],rhat=q['rhat'],Nt=r['params']['Nt'],Nr=r['params']['Nr']))
    b.errorbar([.5,1,1.5,2],yy,yerr=ee,marker=m,ls='-',color=col,lw=1.1,ms=4,capsize=2,label=f'{n}')
b.set(xlabel=r'Coupling $g_N$',ylabel=r'Reduction from $g_N=0.5$ (%)',xticks=[.5,1,1.5,2],ylim=(-2,47))
b.legend(title='Rank',title_fontsize=7,ncol=2,fontsize=7,loc='lower right')
# Mean local extreme eigenvalues, rather than absolute sampled bounds.
for n in [8,12,16,24,32,48]:
    r=by_tag[f'final_IR_mu005_N{n}'];o=r['observables'];p=r['params']
    assert p['n_dof']==p['Nt']*p['Nr']*(3*n*n-8)
    lo,hi=o['rho2_eig_min']['mean'],o['rho2_eig_max']['mean']
    c.plot([n,n],[lo,hi],color=colors[0],lw=1.3)
    c.scatter([n,n],[lo,hi],color=colors[0],s=15)
    c.scatter(n,o['rho2_eig_mean']['mean'],color=colors[3],s=18,marker='D')
    summary['area_diagnostics'].append(dict(N=n,mean_min=lo,mean_max=hi,
        mean_area=o['rho2_eig_mean']['mean'],rhat_area=o['rho2_eig_mean']['rhat'],rhat_monopole=o['phi0_tadpole']['rhat'],
        mass2_fit=o['m2_grav_ir'],dof_per_chain=p['n_dof'],chains=p['n_chains']))
c.axhline(1,color='0.5',ls=':',lw=.6)
c.set(xlabel=r'Matrix rank $N$',ylabel='Normalized local area eigenvalue',xlim=(5,51),ylim=(0,2.9),xticks=[8,16,24,32,48])
c.text(.45,.88,'Mean extrema',transform=c.transAxes,fontsize=7)
c.text(.45,.79,'Diamonds mark mean area',transform=c.transAxes,fontsize=7)
# Convergence is observable-specific.
for prefix,col,m,field,lab in [
    ('final_IR_mu005_N',colors[0],'o','rho2_eig_mean','Area, scaled coupling'),
    ('final_UV_gN2.00_N',colors[1],'s','rho2_eig_mean','Area, fixed coupling'),
    ('final_IR_mu005_N',colors[2],'^','phi0_tadpole','Monopole, scaled coupling')]:
    rr=sorted((r for r in runs if r['tag'].startswith(prefix)),key=lambda r:r['params']['N'])
    d.plot([r['params']['N'] for r in rr],[r['observables'][field]['rhat'] for r in rr],m+'-',color=col,lw=1.1,ms=4,label=lab)
d.set(xlabel=r'Matrix rank $N$',ylabel=r'Chain diagnostic $\hat{R}$',xlim=(5,51),ylim=(.98,1.5),xticks=[8,16,24,32,48])
d.legend(loc='upper left',fontsize=7)
panel_labels(axes)
save(fig,str(HERE.parent/'figs/fig6_modal_b200_hmc'))
(HERE/'reviewed_hmc_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
