all: centipede-tutorial.pdf centipede-tutorial.html

R_OPTS=--no-save --no-restore --no-init-file --no-site-file # vanilla, but with --environ

centipede-tutorial.html: centipede-tutorial.Rmd
	R ${R_OPTS} -e "rmarkdown::render('centipede-tutorial.Rmd', output_format = 'html_document')"

centipede-tutorial.pdf: centipede-tutorial.Rmd
	R ${R_OPTS} -e "rmarkdown::render('centipede-tutorial.Rmd', output_format = 'pdf_document')"

clean:
	rm -f centipede-tutorial.pdf centipede-tutorial.html
