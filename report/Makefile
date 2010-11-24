all:
	pdflatex report.tex
	bibtex report
	pdflatex report.tex
	pdflatex report.tex
	makeglossaries report
	pdflatex report.tex
	open -a Skim report.pdf

letter:
	pdflatex letter.tex
	open -a Skim letter.pdf

nonstopmode:
	pdflatex -interaction nonstopmode report.tex
	bibtex report
	pdflatex -interaction nonstopmode report.tex
	pdflatex -interaction nonstopmode report.tex
	makeglossaries report
	pdflatex -interaction nonstopmode report.tex
	open -a Skim report.pdf

clean:
	-rm -f *.blg *.lof *.log *.lot *.pdf *.toc *.bbl *.aux *glg *.glo *.gls *.ist *.acn *.acr *.alg
