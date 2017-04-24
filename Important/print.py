import datetime
import os
import sys
import shutil
from collections import OrderedDict

header = r"""
\documentclass{article}
\usepackage{amsmath}
\usepackage{minted}
\usepackage{fontspec}
\usepackage{xeCJK}
\usepackage{titlesec}
\usepackage{hyperref}

\newcommand{\sectionbreak}{\clearpage}
\newcommand{\subsectionbreak}{\clearpage}
\setmonofont{Fantasque Sans Mono}
\setCJKmainfont{Microsoft YaHei}
\setCJKmonofont{Microsoft YaHei}
\setminted[cpp]{linenos, breaklines, tabsize=2}

\title{SAJIHPTS Code Snippets}
\date{\today}
\author{SAJIHPTS Team}

\begin{document}
  \maketitle
  \newpage

  \tableofcontents
  \newpage
"""

footer = r"""
\end{document}
"""

os.chdir(sys.path[0])

if os.path.exists("dist"):
	shutil.rmtree("dist")
os.makedirs("dist")
os.chdir("dist")

dic = OrderedDict({})
dicT = OrderedDict({})

def texencode(s):
	return s.replace('&', r'\&').replace('_', r'\_')

for root, dirs, files in os.walk('..', topdown=True):
	dirs[:] = [d for d in dirs if d[0] != '.']
	for f in files:
		if(f.endswith(".cpp")):
			dic.setdefault(root, []).append(f)
		if(f.endswith(".txt")):
			dicT.setdefault(root, []).append(f)

with open('out.tex', 'w', encoding="utf8") as fh:
	fh.write("% Auto Generated on " + datetime.datetime.today().isoformat(' '))
	fh.write(header)
	for d in dic:
		fh.write("""
		\\section{%s}
		""" % texencode(d.split('\\')[-1]))
		for f in dic[d]:
			fh.write("""
			\\subsection{%s}
			\\inputminted{cpp}{%s}
			""" % (texencode(f.split('.')[0]), '"' + os.path.join(d, f).replace('\\', '/')) + '"')
	for d in dicT:
		fh.write("""
		\\section{%s}
		""" % texencode(d.split('\\')[-1]))
		for f in dicT[d]:
			fh.write("""
			\\subsection{%s}
			\\inputminted{text}{%s}
			""" % (texencode(f.split('.')[0]), '"' + os.path.join(d, f).replace('\\', '/')) + '"')
	fh.write(footer)

os.system("xelatex -shell-escape out.tex")
os.system("xelatex -shell-escape out.tex")
