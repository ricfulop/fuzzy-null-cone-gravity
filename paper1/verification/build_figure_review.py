"""Extract current manuscript figures/captions into a standalone LaTeX review.

Run after compiling paper_reading.tex so section references can be resolved.
Build the emitted figures_review.tex twice with pdflatex. No numeric data changes.
"""
from pathlib import Path
import re
ROOT=Path(__file__).resolve().parents[1]
s=(ROOT/'paper.tex').read_text();aux=(ROOT/'paper_reading.aux').read_text()
refs=dict(re.findall(r'\\newlabel\{([^}]+)\}\{\{([^}]+)\}',aux))
macros=s[s.index('\\newcommand'):s.index('\\emergencystretch')]
head=r'''\documentclass[11pt,letterpaper]{article}
\usepackage[margin=0.7in]{geometry}
\usepackage{graphicx,amsmath,amssymb,bm,mathrsfs,caption}
\usepackage[colorlinks=true,linkcolor=blue]{hyperref}
\captionsetup{font=small,labelfont=bf}
\setlength{\parindent}{0pt}
\emergencystretch=1.5em
'''+macros+'\n\\begin{document}\n'
blocks=[]
for i,m in enumerate(re.finditer(r'\\begin\{figure\*\}.*?\\end\{figure\*\}',s,re.S)):
    x=m.group()
    x=re.sub(r'\\begin\{figure\*\}\[[^]]+\]',r'\\begin{figure}[!ht]',x).replace('\\end{figure*}','\\end{figure}')
    x=x.replace('width=\\textwidth','width=\\linewidth,height=6.6in,keepaspectratio')
    x=x.replace('\\cite{Page1976}','(Page, 1976, Table I)')
    x=re.sub(r'\\ref\{(sec:[^}]+)\}',lambda q:refs[q.group(1)],x)
    if i==0:x='\\textbf{Figure review: finite-state null-cone model}\\par\n\\small Main figures 1--6; additional benchmarks 7--11. Captions match the manuscript.\\normalsize\n'+x
    blocks.append(x+'\n\\clearpage\n')
(ROOT/'figures_review.tex').write_text(head+'\n'.join(blocks)+'\\end{document}\n')
print('Wrote',len(blocks),'captioned figures.')
