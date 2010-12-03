all:
	pdflatex test.tex
	makeglossaries test
	pdflatex test.tex
	makeglossaries test
	pdflatex test.tex
	pdflatex test.tex
	open -a Skim test.pdf

clean:
	-rm -f *.blg *.lof *.log *.lot *.pdf *.toc *.bbl *.aux *glg *.glo *.gls *.ist *.acn *.acr *.alg
