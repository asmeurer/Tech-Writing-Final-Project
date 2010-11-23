all:
	pdflatex report.tex
	bibtex report
	pdflatex report.tex
	pdflatex report.tex
	open -a Skim report.pdf

clean:
	rm *.blg *.lof *.log *.lot *.pdf *.toc *.bbl *.aux
