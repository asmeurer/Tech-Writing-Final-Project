all: latex

latex:
	pdflatex presentation.tex
	pdflatex presentation.tex
	open -a Skim presentation.pdf

nonstopmode:
	pdflatex -interaction=nonstopmode presentation.tex
	pdflatex -interaction=nonstopmode presentation.tex
	open -a Skim presentation.pdf

fast:
	pdflatex -interaction=nonstopmode presentation.tex
	open -a Skim presentation.pdf

clean:
	-rm -f *.blg *.lof *.log *.lot *.pdf *.toc *.bbl *.aux *glg *.glo *.gls *.ist *.acn *.acr *.alg
