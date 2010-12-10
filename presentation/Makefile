all: extras latex

latex:
	pdflatex presentation.tex
	pdflatex presentation.tex
	open -a Skim presentation.pdf

nonstopmode: extrasnonstopmode
	pdflatex -interaction=nonstopmode presentation.tex
	pdflatex -interaction=nonstopmode presentation.tex
	open -a Skim presentation.pdf

fast:
	pdflatex -interaction=nonstopmode presentation.tex
	open -a Skim presentation.pdf

clean:
	-rm -f *.blg *.lof *.log *.lot *.pdf *.toc *.bbl *.aux *glg *.glo *.gls *.ist *.acn *.acr *.alg

extras:
	pdflatex algorithmfull.tex
	pdflatex algorithmfull.tex
	pdflatex polyrischdenocancel2full.tex
	pdflatex polyrischdenocancel2full.tex

extrasnonstopmode:
	pdflatex -interaction=nonstopmode algorithmfull.tex
	pdflatex -interaction=nonstopmode algorithmfull.tex
	pdflatex -interaction=nonstopmode polyrischdenocancel2full.tex
	pdflatex -interaction=nonstopmode polyrischdenocancel2full.tex
