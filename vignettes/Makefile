all:
	"$(R_HOME)/bin$(R_ARCH_BIN)/Rscript" -e "library(knitr); knit('fmcsR.Rmd')"
	pandoc --self-contained --output fmcsR.html -H header.html --from markdown --smart --table-of-contents --template html.template fmcsR.md;

clean:
	rm -f *.md
