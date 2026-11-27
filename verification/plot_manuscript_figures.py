#!/usr/bin/env python3
"""Rebuild the manuscript artwork without changing simulation or result files.

General manuscript layout at 178 mm; these dimensions are not a claim of journal
acceptance. PDF and SVG retain vector lines/text. PNG is a 300 dpi review copy.

The six main figures explain geometry, the exact B=4 state count, the archived
matched-control roughness result, observable-specific convergence, two-loop
cutoff dependence, and the known Schwarzschild emission benchmark. Supplementary
panels retain archived comparisons but identify illustrative/conditional results.
No observations are generated, dropped, or averaged for visual convenience.
The corrected pair-state sum is evaluated in memory using 2*omega <= B; the
archived single-particle-cutoff JSON is deliberately left unchanged.
"""
from pathlib import Path
import json
import math
import runpy
import sys
import warnings
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.ticker import ScalarFormatter
from scipy.integrate import simpson
from scipy.special import exp1

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
OUT=ROOT/'figs'
sys.path.insert(0,str(HERE))
BLUE='#0072B2'; RUST='#B44F00'; TEAL='#007D67'; PURPLE='#8064A2'
INK='#20262B'; GRAY='#65717A'; LIGHT='#EAF2F5'; GRID='#D6DEE3'
COLORS=[BLUE,RUST,TEAL,PURPLE]
WIDTH=178/25.4
STEMS=[]

def read(name):
    return json.loads((HERE/name).read_text())

def style():
    plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Liberation Sans','DejaVu Sans'],
        'font.size':9,'axes.labelsize':9,'axes.titlesize':10,'legend.fontsize':8,
        'xtick.labelsize':8,'ytick.labelsize':8,'mathtext.fontset':'stixsans',
        'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none',
        'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,
        'axes.grid':False,'lines.linewidth':1.5,'lines.markersize':4.5,
        'xtick.direction':'out','ytick.direction':'out','xtick.major.size':3,
        'ytick.major.size':3,'legend.frameon':False,'figure.facecolor':'white',
        'axes.facecolor':'white','text.color':INK,'axes.labelcolor':INK,
        'xtick.color':INK,'ytick.color':INK,'savefig.transparent':False})


def panel(ax,letter,title):
    ax.set_title(title,loc='left',pad=12,fontweight='medium',fontsize=10)
    ax.text(-.12,1.06,letter,transform=ax.transAxes,weight='bold',fontsize=12,
            ha='right',va='bottom',clip_on=False)


def save(fig,stem):
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    bounds=fig.bbox
    for t in fig.findobj(matplotlib.text.Text):
        if not t.get_visible() or not t.get_text().strip(): continue
        b=t.get_window_extent(renderer)
        # Ignore invisible off-scale tick labels; test other text against page.
        if t in [q for ax in fig.axes for axis in (ax.xaxis,ax.yaxis) for q in axis.get_ticklabels()]:
            continue
        if b.x0 < -2 or b.y0 < -2 or b.x1>bounds.x1+2 or b.y1>bounds.y1+2:
            warnings.warn(f'{stem}: text outside canvas: {t.get_text()!r}')
    for ext in ['pdf','svg','png']:
        fig.savefig(OUT/f'{stem}.{ext}',dpi=300,facecolor='white',bbox_inches=None)
    STEMS.append(stem)
    plt.close(fig)


def geometry():
    fig=plt.figure(figsize=(WIDTH,3.05))
    # Manual positions preserve the drawing geometry and final text size.
    a=fig.add_axes([.025,.17,.29,.66]);a.axis('off')
    a.set(xlim=(-1.36,1.38),ylim=(-1.30,1.38),aspect='equal')
    a.fill([0,1,0,-1],[-1,0,1,0],fc=LIGHT,zorder=0)
    a.plot([0,1,0,-1,0],[-1,0,1,0,-1],color=BLUE,lw=1.6)
    th=np.linspace(0,np.pi,160)
    a.plot(np.cos(th),.16*np.sin(th),color=GRAY,ls=(0,(3,2)),lw=.8)
    a.plot(np.cos(th),-.16*np.sin(th),color=BLUE,lw=1.1)
    a.scatter([0,0],[-1,1],s=18,color=INK,zorder=4)
    a.text(0,1.15,'Later event',ha='center',fontsize=8.5)
    a.text(0,-1.20,'Earlier event',ha='center',fontsize=8.5)
    a.text(0,.04,'Boundary sphere',ha='center',fontsize=8,
           bbox={'facecolor':LIGHT,'edgecolor':'none','pad':1})
    a.annotate('Lightlike\nboundary',xy=(-.62,.38),xytext=(-1.33,.76),
               fontsize=8,ha='left',arrowprops={'arrowstyle':'-','color':GRAY,'lw':.65})
    a.annotate('',xy=(1.22,.83),xytext=(1.22,.25),
               arrowprops={'arrowstyle':'->','color':GRAY,'lw':.8})
    a.text(1.22,.15,'Time',ha='center',fontsize=8)
    b=fig.add_axes([.365,.47,.275,.31]);b.axis('off');b.set(xlim=(0,1),ylim=(-.5,.5))
    u=np.linspace(0,1,300)
    b.axvspan(.19,.83,color=LIGHT)
    for h in [-1,-.5,0,.5,1]: b.plot(u,h*(.43-.29*u**1.5),color=TEAL,lw=1)
    b.text(.5,-.52,'Neighboring light rays',ha='center',fontsize=8)
    t=fig.add_axes([.365,.24,.275,.16]);t.set(xlim=(0,1),ylim=(0,1.12))
    z=(u-.19)/.64;inside=(z>0)&(z<1);weight=np.zeros_like(z)
    weight[inside]=np.exp(1-1/(4*z[inside]*(1-z[inside])))
    t.fill_between(u,weight,color=LIGHT);t.plot(u,weight,color=TEAL,lw=1.4)
    t.set(xticks=[],yticks=[],xlabel='Position along a light ray');t.spines['left'].set_visible(False)
    t.text(.50,1.17,'Smooth averaging weight',ha='center',fontsize=8)
    c=fig.add_axes([.71,.20,.27,.61]);c.axis('off')
    for y,h,txt in [(.58,.39,'Finite angular resolution\nTotal excitation budget\nFinite corner sectors'),(.21,.20,'Finite basis of states')]:
        box=patches.FancyBboxPatch((.01,y),.96,h,boxstyle='round,pad=.025,rounding_size=.025',
             transform=c.transAxes,fc=LIGHT if y>.5 else 'white',ec=BLUE,lw=.8,clip_on=False)
        c.add_patch(box);c.text(.49,y+h/2,txt,ha='center',va='center',fontsize=8.5,linespacing=1.45,transform=c.transAxes)
    for lo,hi in [(.44,.55),(.045,.175)]:
        c.annotate('',xy=(.49,lo),xytext=(.49,hi),xycoords='axes fraction',
                   arrowprops={'arrowstyle':'->','color':INK,'lw':.9})
    c.text(.49,-.035,'Finite closed traces',ha='center',va='center',fontsize=9,transform=c.transAxes)
    for x,l,title in [(.025,'A','Causal diamond'),(.365,'B','Clock prescription'),(.71,'C','State-space assumptions')]:
        fig.text(x,.94,l,fontsize=12,weight='bold',va='top')
        fig.text(x+.035,.94,title,fontsize=9.5,va='top')
    fig.text(.50,.045,'Schematic geometry  •  Finiteness is conditional on the stated assumptions',ha='center',fontsize=8,color=GRAY)
    save(fig,'review_01_causal_diamond')


def state_count():
    fig, (a,b)=plt.subplots(1,2,figsize=(WIDTH,3.35),gridspec_kw={'width_ratios':[1.05,1]})
    fig.subplots_adjust(left=.09,right=.985,bottom=.20,top=.81,wspace=.38)
    k=np.arange(1,9)
    for j,col,m in [(2,BLUE,'o'),(3,RUST,'s'),(4,TEAL,'^')]:
        w=k+(j-1)*(j+2)/(2*k)
        a.plot(k,w,color=col,marker=m,ms=4,mfc='white',lw=1.1,label=f'$j={j}$')
        ok=w<=4
        a.plot(k[ok],w[ok],ls='none',marker=m,ms=6,color=col)
    a.axhspan(0,4,color=LIGHT,zorder=-2);a.axhline(4,color=GRAY,ls='--',lw=.8)
    a.text(7.9,4.18,'Budget $B=4$',ha='right',fontsize=8)
    a.set(xlim=(.65,8.35),ylim=(0,11),xticks=[1,2,3,4,6,8],yticks=[0,2,4,6,8,10],
          xlabel='Radial mode $k$',ylabel=r'One-graviton cost $\omega_{k,j}$')
    a.legend(loc='upper right',ncol=3,handlelength=.8,columnspacing=.7)
    panel(a,'A','Which modes fit the budget?')
    for row,(kk,pol) in enumerate([(1,'+'),(1,'−'),(2,'+'),(2,'−'),(3,'+'),(3,'−')]):
        b.scatter(range(-2,3),[row]*5,s=38,color=BLUE if pol=='+' else TEAL,
                  marker='o' if pol=='+' else 's',zorder=3)
    b.set(xlim=(-2.6,2.6),ylim=(6.25,-.85),xticks=range(-2,3),yticks=range(6),
          yticklabels=['$k=1$, +','$k=1$, −','$k=2$, +','$k=2$, −','$k=3$, +','$k=3$, −'],
          xlabel='Magnetic quantum number $m$')
    b.tick_params(axis='y',length=0);b.spines[['left','top','right']].set_visible(False)
    b.text(0,6.0,'30 one-graviton states + 1 vacuum = 31',ha='center',fontsize=8.7)
    panel(b,'B','Count the states at fixed corner data')
    fig.text(.54,.035,'Each marker is one state. Two gravitons cost at least 6, so no pair fits.',ha='center',fontsize=8,color=GRAY)
    assert min(k+(2-1)*(2+2)/(2*k))*2>4
    save(fig,'review_02_state_count')


def roughness():
    d=read('modal_hmc_followup_analysis.json')
    fig,(a,b)=plt.subplots(1,2,figsize=(WIDTH,3.65),gridspec_kw={'width_ratios':[.92,1.05]})
    fig.subplots_adjust(left=.095,right=.99,bottom=.19,top=.83,wspace=.78)
    for n,col,m in zip([8,16,24,32],COLORS,['o','s','^','D']):
        rows=sorted((r for r in d['crosscheck_gN_grid'] if r['N']==n),key=lambda r:r['g_N'])
        x=np.array([r['g_N'] for r in rows]);y=np.array([r['Q2_mean'] for r in rows])
        e=1.96*np.array([r['Q2_se'] for r in rows])
        a.errorbar(x,y,yerr=e,marker=m,color=col,lw=1.25,ms=4,capsize=2,label=f'$N={n}$')
    a.set(xlim=(-.06,2.09),ylim=(.065,.36),xticks=[0,.5,1,1.5,2],yticks=[.1,.2,.3],
          xlabel='Coupling $g_N$',ylabel=r'Mean roughness $\langle Q_2\rangle$')
    a.legend(loc='upper right',ncol=2,fontsize=8,handlelength=.9,columnspacing=.7)
    panel(a,'A','Roughness decreases with coupling')
    vals=[r['smoothing_vs_measured_g0'] for r in d['primary_4rank_table']]
    vals += [d['lattice_12x8_test']['smoothing_12x8_vs_measured_g0'],d['extended_n32_smoothing_crosscheck']]
    labels=['$N=8$  |  8 × 4','$N=16$  |  8 × 4','$N=24$  |  8 × 4','$N=32$  |  12 × 4',
            '$N=24$  |  12 × 8','$N=32$  |  extended']
    yy=[0,1,2,3,4.6,5.6]
    for i,(v,y) in enumerate(zip(vals,yy)):
        mu=v['smoothing_pct'];col=COLORS[i] if i<4 else INK
        b.errorbar(mu,y,xerr=[[mu-v['ci95_mc_lo']],[v['ci95_mc_hi']-mu]],fmt=['o','s','^','D','v','D'][i],
                   mfc=col if i<4 else 'white',mec=col,ecolor=col,ms=5,capsize=3,lw=1.15)
    b.axhline(3.8,color=GRID,lw=.8)
    b.set(yticks=yy,yticklabels=labels,xlim=(47.48,48.18),xticks=[47.5,47.7,47.9,48.1],
          ylim=(6.15,-.65),xlabel='Reduction at $g_N=2$ (%)')
    b.tick_params(axis='y',length=0);b.spines['left'].set_visible(False)
    b.text(.02,.29,'Follow-up checks',transform=b.transAxes,fontsize=8,color=GRAY)
    panel(b,'B','Similar effect across tested sizes')
    fig.text(.52,.035,'Bars: approximate 95% Monte Carlo intervals  •  Right panel uses an expanded percentage scale',ha='center',fontsize=8,color=GRAY)
    save(fig,'fig_followup_roughness')


def diagnostics():
    d=read('modal_hmc_followup_analysis.json')
    rr=[d[k] for k in ['extended_n32_batch1','extended_n32_batch2','extended_n48_batch1','extended_n48_batch2']]
    rr.append(d['lattice_12x8_test']['lat_12x8_g2'])
    labels=['$N=32$, global','$N=32$, +local','$N=48$, global','$N=48$, +local','$N=24$, 12 × 8']
    fig,axs=plt.subplots(1,3,figsize=(WIDTH,3.35),gridspec_kw={'width_ratios':[1,1,1.1]})
    fig.subplots_adjust(left=.20,right=.985,top=.80,bottom=.22,wspace=.36)
    for a,field in zip(axs[:2],['Q2','rho2_eig_mean']):
        xx=[r['diagnostics'][field]['rhat_rank'] for r in rr]
        a.axvline(1.01,color=GRAY,ls='--',lw=.8)
        a.plot(xx,range(5),'o',color=BLUE,ms=5)
        a.set(ylim=(4.6,-.6),yticks=range(5),yticklabels=labels if a==axs[0] else [],xlabel=r'Split-$\hat R$')
        a.tick_params(axis='y',length=0);a.spines['left'].set_visible(False)
        for y,x in enumerate(xx): a.annotate(f'{x:.4f}',(x,y),xytext=(5,0),textcoords='offset points',fontsize=7.8,va='center')
    # Label values explicitly so close-to-one differences remain inspectable.
    axs[0].set(xlim=(.999,1.0185),xticks=[1,1.01]);axs[0].set_xticklabels(['1.00','1.01'])
    axs[1].set(xlim=(.99,1.49),xticks=[1,1.2,1.4]);axs[1].set_xticklabels(['1.0','1.2','1.4'])
    for metric,dy,col,m,lab in [('ess_bulk',-.1,BLUE,'o','Bulk'),('ess_tail',.1,RUST,'s','Tail')]:
        xx=[r['diagnostics']['rho2_eig_mean'][metric] for r in rr]
        axs[2].scatter(xx,np.arange(5)+dy,s=24,color=col,marker=m,label=lab)
    axs[2].axvline(400,color=GRAY,ls='--',lw=.8)
    axs[2].set(xscale='log',xlim=(10,20000),xticks=[10,100,1000,10000],ylim=(4.6,-.6),
               yticks=range(5),yticklabels=[],xlabel='Effective sample size')
    axs[2].tick_params(axis='y',length=0);axs[2].spines['left'].set_visible(False)
    axs[2].legend(loc='upper center',bbox_to_anchor=(.5,-.25),ncol=2,handletextpad=.3,columnspacing=.7)
    for ax,l,title in zip(axs,'ABC',['Roughness','Local mean area','Local-area information']):panel(ax,l,title)
    fig.text(.20,.055,r'Dashed targets: $\hat R<1.01$ and ESS ≥ 400. The larger-lattice area estimate still fails.',fontsize=8,color=GRAY)
    save(fig,'fig6_modal_b200_hmc')


def loops():
    d=read('p2_uv1_2loop_gs_cemz_results.json')['part_a_2loop_finiteness']
    rows=d['fuzzy_only_2loop_divergence_scan'];bb=d['cross_focusing_all_loop_rescue']
    fig,(a,b)=plt.subplots(1,2,figsize=(WIDTH,3.25))
    fig.subplots_adjust(left=.095,right=.985,bottom=.20,top=.82,wspace=.40)
    xx=[r['L_par'] for r in rows]
    for key,col,m,lab in [('I_2loop_derivative_sunset',RUST,'s','Two-loop diagnostic'),('I_1loop_W2',BLUE,'o','One-loop diagnostic')]:
        a.loglog(xx,[r[key] for r in rows],marker=m,color=col,label=lab)
    a.set(xlabel=r'Radial/time cutoff $\Lambda_\parallel/\Lambda_\perp$',ylabel='Diagnostic magnitude',xticks=xx)
    a.xaxis.set_major_formatter(ScalarFormatter());a.minorticks_off();a.legend(loc='upper left',fontsize=8)
    panel(a,'A','Angular cutoff alone is insufficient')
    x=[r['B_cross_focusing'] for r in bb];y=[r['exact_finite_2loop_sunset_sum'] for r in bb]
    b.semilogy(x,y,'o-',color=BLUE)
    for r in bb:
        b.annotate(f"{r['total_allowed_modes_at_Nrad_inf']} modes",(r['B_cross_focusing'],r['exact_finite_2loop_sunset_sum']),
                   xytext=(0,8),textcoords='offset points',ha='center',fontsize=8)
    b.set(xlabel='Single-mode cost cutoff $B$',ylabel='Two-loop diagnostic sum',xticks=x,xlim=(3.3,10.7),ylim=(120,1.3e6))
    b.minorticks_off();panel(b,'B','The regulated sum grows with budget')
    fig.text(.54,.035,'Finite at fixed cutoff does not establish a finite cutoff-removal limit.',ha='center',fontsize=8,color=GRAY)
    save(fig,'review_03_two_loop_budget')


def hawking():
    from review_spectrum import transmission
    w=np.linspace(.06,1.60,155)
    curves=[]
    fig,(a,b)=plt.subplots(1,2,figsize=(WIDTH,3.25))
    fig.subplots_adjust(left=.09,right=.985,bottom=.20,top=.82,wspace=.36)
    for j,col,ls in [(2,BLUE,'-'),(3,RUST,'--'),(4,TEAL,'-.'),(5,PURPLE,':'),(6,GRAY,(0,(5,2)))]:
        t=transmission(j,w);curves.append((2*j+1)/np.pi*w*t/np.expm1(8*np.pi*w))
        a.plot(w,t,color=col,ls=ls,lw=1.4,label=f'$j={j}$')
    total=np.sum(curves,axis=0);alpha=simpson(total,x=w)
    assert abs(alpha/3.836658052716681e-6-1)<1e-7
    a.set(xlim=(0,1.6),ylim=(0,1.04),xticks=[0,.4,.8,1.2,1.6],yticks=[0,.5,1],
          xlabel=r'Frequency $M\omega$',ylabel='Transmission probability')
    a.legend(loc='lower right',ncol=2,fontsize=8)
    b.fill_between(w,total*1e6,color=LIGHT);b.plot(w,total*1e6,color=BLUE,lw=1.6)
    b.set(xlim=(0,1.6),ylim=(0,None),xticks=[0,.4,.8,1.2,1.6],xlabel=r'Frequency $M\omega$',ylabel=r'$10^6\,d\alpha/d(M\omega)$')
    b.text(.97,.95,'Integrated coefficient\n'+r'$\alpha=3.83666\times10^{-6}$'+'\nPage: '+r'$3.8\times10^{-6}$',
           transform=b.transAxes,ha='right',va='top',fontsize=8.5,linespacing=1.5)
    panel(a,'A','Curvature filters the outgoing waves')
    panel(b,'B','Thermally weighted graviton spectrum')
    fig.text(.54,.035,r'Known Schwarzschild benchmark  •  $G=c=\hbar=1$  •  Radiated power $P=\alpha/M^2$',ha='center',fontsize=8,color=GRAY)
    save(fig,'review_04_graviton_spectrum')


def corrected_pair_sum():
    """Art-only evaluation of the corrected projector; no new data file written."""
    out=[]
    for s in [100.,300.,1000.,3000.,10000.,30000.]:
        J=int(np.floor(np.sqrt(s/np.log(3))));limit=s/2;mb=0j;mi=0j
        for j in range(2,J+1):
            t2=(j-1)*(j+2)/2;disc=limit**2-4*t2
            if disc>=0:
                lo=max(1,math.ceil((limit-np.sqrt(disc))/2));hi=math.floor((limit+np.sqrt(disc))/2)
                k=np.arange(lo,hi+1,dtype=float);omega=k+t2/k
                assert np.all(2*omega<=s+1e-8)
                mb+=(2*j+1)*np.sum(np.exp(-2j*omega/s)/omega**2)/s**2
            kd=max(4000,int(1.5*s));k=np.arange(1,kd+1,dtype=float);omega=k+t2/k
            z=2j*(kd+.5)/s;tail=(np.exp(-z)-z*exp1(z))/(kd+.5)
            mi+=(2*j+1)*(np.sum(np.exp(-2j*omega/s)/omega**2)+tail)/s**2
        out.append([s,mb.real,abs(mb.imag),abs(mi-mb)])
    return np.array(out)


def strengthen_contrast(fig):
    """Keep legacy category encodings while darkening low-contrast line colors."""
    from matplotlib.colors import to_rgba
    mapping={'#E69F00':RUST,'#009E73':TEAL,'#D55E00':'#AD397B','#CC79A7':PURPLE,'#56B4E9':'#326E92'}
    def changed(c):
        try:r=np.asarray(to_rgba(c),dtype=float)
        except (ValueError,TypeError):return c
        for old,new in mapping.items():
            if np.allclose(r[:3],to_rgba(old)[:3],atol=1e-6):
                q=list(to_rgba(new));q[3]=r[3];return q
        return c
    for ax in fig.axes:
        for line in ax.lines:line.set_color(changed(line.get_color()))
        for col in ax.collections:
            for getter,setter in [(col.get_facecolors,col.set_facecolors),(col.get_edgecolors,col.set_edgecolors)]:
                vals=getter()
                if len(vals):setter([changed(c) for c in vals])
        for patch in ax.patches:
            patch.set_facecolor(changed(patch.get_facecolor()))
            patch.set_edgecolor(changed(patch.get_edgecolor()))


def supplementary():
    """Restyle existing benchmark panels, preserving data except stated fixes."""
    import orx_figstyle
    captured=[]
    old_save=orx_figstyle.save
    orx_figstyle.save=lambda fig,stem: captured.append((fig,Path(stem).name))
    scripts=['plot_path2_frontier.py','plot_v2_unification_and_clock.py',
             'plot_path2_breakthrough.py','plot_path2_uv_frontier.py','plot_original_problem_solved.py']
    for script in scripts: runpy.run_path(str(HERE/script),run_name='figure_source')
    orx_figstyle.save=old_save
    titles={
      'fig2_unified_matrix_gravity':['Selected matrix-product errors','Known black-hole potentials','Truncated boost test','Illustrative radius trajectory'],
      'fig5_v2_unification_and_clock':['Cubic-vertex convergence','Projected dynamics','Corrected two-particle cutoff','Quadratic-clock bounds'],
      'fig1_breakthrough_1loop':['Quadratic clock bounds','Angular-regulator dependence','Schwarzschild entropy benchmark','Entropy-coefficient residuals'],
      'fig3_uv_completion':['Transverse-only loop diagnostics','Assumed rank-weight model','Model crossover in AdS','Imposed pixel-size criterion'],
      'fig4_original_resolution':['Gaussian reference amplitude','A finite Gram-matrix example','Regulated correlation example','Raychaudhuri backreaction test']}
    legends={
     'fig2_unified_matrix_gravity':[['One-sided product','Symmetric wave, $j=4$','Symmetric wave, $j=2$','Symmetric Kerr term'],
       ['Scalar reference','Regge–Wheeler','Zerilli'],[r'$\eta=0.5$',r'$\eta=1.0$',r'$\eta=1.5$'],
       ['Imposed minimum','Classical focusing','Bounce model']],
     'fig5_v2_unification_and_clock':[['One-sided (2,2,2)','Symmetric (2,2,2)','Symmetric (4,2,2)'],
       ['Self-energy error','Transition-probability error'],[],['Second moment','Chernoff','Fourth moment']],
     'fig1_breakthrough_1loop':[['Local: floor','Coupled: floor','Local: variance','Coupled: variance'],
       ['Graviton reference','Maxwell reference','Sharp cutoff','Smooth: graviton','Smooth: Maxwell'],
       ['Mixed reference','Mixed estimate','Singlet reference','Singlet estimate','Corner-only term'],
       ['Final / flat estimate','Five-radius fit','Exact identity']],
     'fig3_uv_completion':[['Two-loop diagnostic','One-loop diagnostic','Mode cutoff $B=4$'],
       ['Flat-space weight','AdS weight'],['Classical reference','Coefficient $+106/45$','Coefficient $-3/2$'],
       ['Imposed radius threshold',r'$s_{\min}=\ln(1+d_0)$',r'$s_{\min}=\ln d_0$']],
     'fig4_original_resolution':[['Real part','Leading real asymptote','Imaginary magnitude','Leading log asymptote'],
       ['First 100 magnitudes','40 numerical null residuals'],['Reference $u^{-4}$','Regulated model','Imposed length scale'],
       ['Classical residual','Shear variance only','With the cross-term']]}
    # Standalone labels avoid source-paper equation or chapter numbering.
    ylabs={
      'fig2_unified_matrix_gravity':['Matrix-product error',r'Potential $M^2V$',r'Wavepacket infidelity',r'Radius $\rho/\ell_P$'],
      'fig5_v2_unification_and_clock':['Cubic-vertex error','Observable error','Gaussian amplitude magnitude','Clock-projector failure bound'],
      'fig1_breakthrough_1loop':['Dimensionless clock statistic',r'Logarithmic coefficient $c_{\log}$','Relative subleading entropy',''],
      'fig3_uv_completion':['Diagnostic magnitude',r'$10^{-4}\ln W_N$','Black-hole probability in the model',r'Pixel radius $\rho_{\min}/\ell_P$'],
      'fig4_original_resolution':['Gaussian amplitude magnitude','Gram eigenvalue magnitude','Area-focusing correlator','Raychaudhuri residual']}
    xlabs={
      'fig2_unified_matrix_gravity':['Matrix rank $N$',r'Tortoise coordinate $r_*/M$','Matrix rank $N$','Affine parameter $v$ (model units)'],
      'fig5_v2_unification_and_clock':['Matrix rank $N$','Matrix rank $N$','Dimensionless area $s$','Dimensionless area $s$'],
      'fig1_breakthrough_1loop':['Dimensionless area $s$',r'Cutoff ratio $\alpha$','Dimensionless area $s$','Signed benchmark difference (%)'],
      'fig3_uv_completion':[r'Cutoff ratio $\Lambda_\parallel/\Lambda_\perp$','Matrix rank $N$',r'Temperature $TL$ ($L=25\ell_P$)','Local corner-state dimension $d_0$'],
      'fig4_original_resolution':['Dimensionless area $s$','Eigenvalue index',r'Affine separation $u/\ell_P$',r'Expansion parameter $\hbar$ (model units)']}
    for fig,stem in captured:
        style()
        strengthen_contrast(fig)
        fig.set_layout_engine(None);fig.set_size_inches(WIDTH,6.8)
        for i,ax in enumerate(fig.axes):
            # Remove old panel letters and long incidental annotations.
            for t in list(ax.texts):t.remove()
            ax.set_position([[.10,.66,.37,.25],[.61,.66,.37,.25],[.10,.20,.37,.25],[.61,.20,.37,.25]][i])
            ax.set_xlabel(xlabs[stem][i],fontsize=9);ax.set_ylabel(ylabs[stem][i],fontsize=9)
            ax.tick_params(labelsize=8)
            for line in ax.lines:line.set_markersize(4)
            handles,labels=ax.get_legend_handles_labels()
            if ax.get_legend():ax.get_legend().remove()
            labels=legends[stem][i]
            if labels:
                ax.legend(handles[:len(labels)],labels,loc='upper left',bbox_to_anchor=(0,-.34),
                          fontsize=8,ncol=2 if len(labels)>2 else 1,handlelength=1.5,
                          columnspacing=.75,handletextpad=.5,borderaxespad=0)
            panel(ax,'ABCD'[i],titles[stem][i])
        if stem=='fig2_unified_matrix_gravity':
            ax=fig.axes[2]
            ax.axhspan(5e-16,1e-14,color=LIGHT,zorder=-1)
            ax.text(.98,.10,'Floating-point floor',transform=ax.transAxes,ha='right',fontsize=8,color=GRAY)
        if stem=='fig3_uv_completion':
            fig.axes[1].set_ylim(-2.4,8.2)
        if stem=='fig5_v2_unification_and_clock':
            ax=fig.axes[2];ax.clear();style();q=corrected_pair_sum()
            for y,col,m,lab in [(1,BLUE,'o','Real part'),(2,RUST,'s','Imaginary magnitude'),(3,TEAL,'^','Cutoff-tail difference')]:
                ax.loglog(q[:,0],q[:,y],marker=m,color=col,ms=4,label=lab)
            ax.set(xlabel='Dimensionless area $s$',ylabel='Gaussian amplitude magnitude',xticks=[100,1000,10000])
            ax.minorticks_off();ax.legend(loc='upper left',bbox_to_anchor=(0,-.34),fontsize=8,ncol=2)
            panel(ax,'C',titles[stem][2])
        if stem=='fig4_original_resolution':
            ax=fig.axes[1];ax.axhspan(1e-19,1e-14,color=LIGHT,zorder=-1)
            ax.text(.04,.36,'Floating-point residuals below',transform=ax.transAxes,ha='left',fontsize=8,color=GRAY)
        # One compact scope line replaces promotional claims in the artwork.
        fig.text(.10,.985,'ADDITIONAL BENCHMARKS',fontsize=8,color=GRAY,va='top')
        save(fig,stem)


def main():
    OUT.mkdir(exist_ok=True);style()
    for f in [geometry,state_count,roughness,diagnostics,loops,hawking]:f()
    supplementary()
    print('Rendered',len(STEMS),'vector figures and review PNGs; numerical input files unchanged.')

if __name__=='__main__':main()
