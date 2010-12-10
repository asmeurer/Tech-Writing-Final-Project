all:samples python letter latex

latex:
	pdflatex report.tex
	bibtex report
	pdflatex report.tex
	makeglossaries report
	pdflatex report.tex
	makeglossaries report
	pdflatex report.tex
	pdflatex report.tex
	open -a Skim report.pdf

letter:
	pdflatex letter.tex
	pdflatex letter.tex

letteropen: letter
	open -a Skim letter.pdf

nonstopmode: python
	pdflatex -interaction nonstopmode report.tex
	bibtex report
	pdflatex -interaction nonstopmode report.tex
	makeglossaries report
	pdflatex -interaction nonstopmode report.tex
	makeglossaries report
	pdflatex -interaction nonstopmode report.tex
	pdflatex -interaction nonstopmode report.tex
	open -a Skim report.pdf

fast:
	pdflatex -interaction nonstopmode report.tex
	open -a Skim report.pdf

python:
	pygmentize -f latex polyrischdenocancel2.py >| polyrischdenocancel2.tex

clean:
	-rm -f *.blg *.lof *.log *.lot *.pdf *.toc *.bbl *.aux *glg *.glo *.gls *.ist *.acn *.acr *.alg polyrischdenocancel2.tex

samples:
	xelatex sample1.tex
	xelatex sample1.tex
	xelatex sample2.tex
	xelatex sample2.tex
	xelatex sample3.tex
	xelatex sample3.tex
	xelatex sample4.tex
	xelatex sample4.tex
	xelatex sample5.tex
	xelatex sample5.tex
	xelatex sample6.tex
	xelatex sample6.tex
	xelatex sample7.tex
	xelatex sample7.tex
	xelatex sample8.tex
	xelatex sample8.tex
	xelatex sample9.tex
	xelatex sample9.tex
	xelatex sample10.tex
	xelatex sample10.tex
	xelatex sample11.tex
	xelatex sample11.tex
	xelatex sample12.tex
	xelatex sample12.tex
	xelatex sample13.tex
	xelatex sample13.tex
	xelatex sample14.tex
	xelatex sample14.tex
	xelatex sample15.tex
	xelatex sample15.tex
	xelatex sample16.tex
	xelatex sample16.tex
	xelatex sample17.tex
	xelatex sample17.tex
