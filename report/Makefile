all:
	pdflatex report.tex
	bibtex report
	pdflatex report.tex
	pdflatex report.tex
	makeglossaries report
	pdflatex report.tex
	open -a Skim report.pdf
