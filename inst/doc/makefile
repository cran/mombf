mombfmanual: mombf_manual.Rnw
	echo 'Sweave("mombf_manual.Rnw")' | R --vanilla --quiet
	pdflatex mombf_manual.tex
	bibtex mombf_manual
	pdflatex mombf_manual.tex
	pdflatex mombf_manual.tex
	rm mombf_manual.tex
	rm mombf_manual.log
	rm mombf_manual.aux
	rm mombf_manual.out
	rm mombf_manual.blg
	rm mombf_manual.bbl
	rm mombf_manual-fig*.pdf
	rm mombf_manual-fig*.eps
	rm Rplots.ps
