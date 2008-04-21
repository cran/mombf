mombfmanual: mombf.Rnw
	echo 'Sweave("mombf.Rnw")' | R --vanilla --quiet
	pdflatex mombf.tex
	bibtex mombf
	pdflatex mombf.tex
	pdflatex mombf.tex
	rm mombf.tex
	rm mombf.log
	rm mombf.aux
	rm mombf.out
	rm mombf.blg
	rm mombf.bbl
	rm mombf-fig*.pdf
	rm mombf-fig*.eps
