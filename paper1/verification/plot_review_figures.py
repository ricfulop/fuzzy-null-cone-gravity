"""Render four explanatory vector figures from the accompanying numerical data."""
from pathlib import Path
import csv
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.ticker import FuncFormatter
import numpy as np
from scipy.integrate import simpson

ROOT = Path(__file__).resolve().parents[1]
DATA, OUT = ROOT/'verification/review_data', ROOT/'figs'
DATA.mkdir(exist_ok=True)
OUT.mkdir(exist_ok=True)
INK, MUTED, GRID = '#183449', '#49616C', '#D7E2E3'
TEAL, LIGHT, ORANGE = '#006E71', '#C6E9E5', '#A64723'
INDIGO = '#365A94'
from orx_figstyle import use_style, save as publication_save, WIDE
use_style()
INK = MUTED = "black"

def sheet(title):
    return plt.figure(figsize=(WIDE,3.65),facecolor='white')

def save(fig, name):
    if name in ('01_causal_diamond','02_state_count'):
        positions=(.035,.526) if name.startswith('01') else (.10,.60)
        for x,letter in zip(positions,('A','B')):
            fig.text(x,.90,letter,fontsize=10,fontweight='bold')
    for text in fig.findobj(matplotlib.text.Text):
        text.set_fontsize(min(10,text.get_fontsize()))
    publication_save(fig, str(OUT/('review_'+name)))

def grid(ax):
    ax.grid(False)
    ax.tick_params(length=3,pad=3)

def label(ax,x,y,s,**kw):
    return ax.text(x,y,s,transform=ax.transAxes,**kw)

def geometry():
    fig=sheet('A gravitational clock built from light rays')
    ax=fig.add_axes([.035,.055,.46,.79]); ax.set(xlim=(-1.48,1.55),ylim=(-1.22,1.24),aspect='equal'); ax.axis('off')
    ax.fill([0,1,0,-1],[-1,0,1,0],color='#E8F4F2',zorder=0)
    ax.plot([0,-1,0,1,0],[-1,0,1,0,-1],color=TEAL,lw=2)
    # An ellipse represents the spherical cross-section of the 4D diamond.
    theta=np.linspace(0,np.pi,121)
    ax.plot(np.cos(theta),.15*np.sin(theta),'--',color=MUTED,lw=.9)
    ax.plot(np.cos(theta),-.15*np.sin(theta),color=TEAL,lw=1.25)
    ax.scatter([0,0],[-1,1],s=23,color=INK,zorder=4)
    ax.text(.10,1.00,'Later event',va='center'); ax.text(.10,-1,'Earlier event',va='center')
    ax.text(0,.015,'Boundary sphere',ha='center',va='center',fontsize=8,backgroundcolor='#E8F4F2')
    ax.annotate('Light ray',xy=(-.50,.50),xytext=(-1.41,.76),ha='left',va='center',arrowprops=dict(arrowstyle='-',color=MUTED,lw=.7))
    ax.annotate('',xy=(-1.28,.08),xytext=(-1.28,-.42),arrowprops=dict(arrowstyle='->',color=MUTED,lw=.85))
    ax.text(-1.28,-.48,'Time',ha='center',va='top',color=MUTED,fontsize=7.5)
    ray=fig.add_axes([.526,.49,.43,.29]); ray.set(xlim=(0,1),ylim=(-.52,.52));ray.axis('off')
    x=np.linspace(0,1,260)
    ray.axvspan(.23,.78,color=LIGHT,alpha=.70)
    for h in (-1,-.5,0,.5,1):
        y=h*(.44-.32*x**1.55)
        ray.plot(x,y,color=TEAL,lw=1.15)
    label(ray,.01,1.08,'Light rays focus',fontsize=9.5,fontweight='semibold')
    ray.text(.51,-.49,'Averaging interval',ha='center',va='center',fontsize=7.5)
    taper=fig.add_axes([.526,.14,.43,.22]);taper.set(xlim=(0,1),ylim=(-.04,1.18))
    u=(x-.23)/.55; inside=(u>0)&(u<1);g=np.zeros_like(u)
    g[inside]=np.exp(1-1/(4*u[inside]*(1-u[inside])))
    taper.fill_between(x,g,color=LIGHT);taper.plot(x,g,color=TEAL,lw=1.5)
    taper.set(xticks=[],yticks=[]);taper.spines['left'].set_visible(False)
    taper.set_xlabel('Position along the light ray',labelpad=5)
    taper.set_ylabel('Weight',labelpad=3)
    label(taper,.01,1.10,'Smooth clock weight',fontsize=9.5,fontweight='semibold')
    save(fig,'01_causal_diamond')

def count():
    fig=sheet('At a budget of four, the spin-2 sector has 31 states')
    ax=fig.add_axes([.10,.22,.45,.59]);ax.set(xlim=(.6,6.4),ylim=(1.55,5.5))
    ks, js=np.meshgrid(np.arange(1,7),np.arange(2,6))
    cost=ks+(js-1)*(js+2)/(2*ks)
    good=cost<=4
    assert sum((2*(2*js[good]+1)))==30 and (2*cost.min()>4)
    ax.scatter(ks[~good],js[~good],marker='x',color='#92A2A6',s=27,lw=1)
    ax.scatter(ks[good],js[good],s=55,color=TEAL,zorder=3)
    ax.set(xticks=range(1,7),yticks=range(2,6),xlabel='Radial mode  $k$',ylabel='Angular mode  $j$')
    ax.tick_params(length=0)
    label(ax,0,1.08,'Allowed single-graviton modes',fontsize=9.6,fontweight='semibold')
    ax.scatter([],[],s=28,color=TEAL,label='Allowed')
    ax.scatter([],[],s=22,marker='x',color='#829599',label='Over budget')
    ax.legend(loc='upper center',bbox_to_anchor=(.50,-.25),ncol=2,frameon=False,fontsize=7.4,handletextpad=.25,columnspacing=.9,borderaxespad=0)
    side=fig.add_axes([.60,.15,.37,.64]);side.axis('off')
    label(side,0,1.0,'Count the independent states',va='top',fontsize=9.6,fontweight='semibold')
    label(side,0,.81,'3 radial modes\n5 angular orientations\n2 polarizations',fontsize=10.5,va='top',linespacing=1.48)
    side.plot([0,.91],[.35,.35],color=GRID,transform=side.transAxes,lw=1)
    label(side,0,.29,'30 one-graviton states',fontsize=11,color=TEAL,fontweight='semibold',va='top')
    label(side,0,.14,'+ 1 vacuum state',fontsize=10.5,va='top')
    label(side,0,-.05,'Two gravitons cost at least 6.',color=MUTED,fontsize=8,va='top')
    save(fig,'02_state_count')

def loops():
    fig=sheet('The two-loop sum grows with the excitation budget')
    d=np.genfromtxt(DATA/'two_loop.csv',delimiter=',',names=True)
    ax=fig.add_axes([.155,.19,.73,.61]);grid(ax)
    x,y=d['budget'],d['sum']
    ax.set_yscale('log');ax.set(xlim=(3.5,10.95),ylim=(130,1.2e6),xticks=x,yticks=[1e3,1e4,1e5,1e6])
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v,pos:f'{v:,.0f}'))
    ax.minorticks_off();ax.set(xlabel='Allowed excitation budget  $B$',ylabel='Two-loop sum')
    ax.plot(x,y,lw=1.5,color=TEAL);ax.scatter(x,y,s=41,color=TEAL,edgecolor='white',lw=.8,zorder=3)
    for a,b in zip(x,y):
        ax.annotate(f'{b:,.0f}',(a,b),xytext=(0,8),ha='center',textcoords='offset points',fontsize=8.6,fontweight='semibold')
    save(fig,'03_two_loop_budget')

def hawking():
    fig=sheet('Gravitons emitted by a nonrotating black hole')
    d=np.genfromtxt(DATA/'spectrum.csv',delimiter=',',names=True)
    w,total=d['M_omega'],d['total']
    a=simpson(total,x=w)
    ax=fig.add_axes([.105,.20,.53,.61]);grid(ax)
    ax.fill_between(w,1e6*total,color=LIGHT)
    ax.plot(w,1e6*total,color=TEAL,lw=1.8)
    ax.set(xlim=(0,1.60),ylim=(0,None),xticks=[0,.4,.8,1.2,1.6],xlabel='Frequency  $M\\omega$',ylabel=r'$10^6\,d\alpha/d(M\omega)$')
    label(ax,0,1.08,'Power per unit frequency',fontsize=9.6,fontweight='semibold')
    bx=fig.add_axes([.67,.15,.30,.65]);bx.axis('off')
    label(bx,0,1,'Integrated power',va='top',fontsize=9.6,fontweight='semibold')
    label(bx,0,.79,rf'${a*1e6:.5f}\times10^{{-6}}$',fontsize=16,color=TEAL,va='top')
    label(bx,0,.61,r'Power $P=α/M^2$',fontsize=7.8,va='top')
    bx.plot([0,.96],[.49,.49],transform=bx.transAxes,color=GRID,lw=.9)
    label(bx,0,.41,r'$3.8\times10^{-6}$',fontsize=14,color=INK,va='top')
    label(bx,0,.26,"Page (1976), Table I\nPublished total, in the same units",fontsize=7.65,va='top',linespacing=1.35,color=MUTED)
    save(fig,'04_graviton_spectrum')

def prepare_data():
    source=ROOT/'verification'
    rows=json.loads((source/'p2_uv1_2loop_gs_cemz_results.json').read_text())['part_a_2loop_finiteness']['cross_focusing_all_loop_rescue']
    with (DATA/'two_loop.csv').open('w') as f:
        writer=csv.writer(f);writer.writerow(['budget','modes','sum'])
        writer.writerows((r['B_cross_focusing'],r['total_allowed_modes_at_Nrad_inf'],r['exact_finite_2loop_sunset_sum']) for r in rows)
    from review_spectrum import main as spectrum
    spectrum()

if __name__ == '__main__':
    prepare_data()
    for render in (geometry,count,loops,hawking):
        render()
