mombfmanual: mombf.Rnw
	echo 'Sweave("mombf.Rnw")' | R --vanilla --quiet
	R CMD pdflatex mombf.tex
	R CMD bibtex mombf
	R CMD pdflatex mombf.tex
	R CMD pdflatex mombf.tex
	rm mombf.tex
	rm mombf.log
	rm mombf.aux
	rm mombf.out
	rm mombf.blg
	rm mombf.bbl
	rm mombf-fig*.pdf
	rm mombf-fig*.eps
	rm Rplots.*

