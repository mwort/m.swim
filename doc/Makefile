
# A Makefile to publish and convert the documentation
include ../test/Makefile

URL = ~/www/m.swim
MD_SRC = home.md developing_debugging.md

help_pages:
	GRASS_BATCH_JOB=./build_help_pages.grass $(GRASS) $(GRASSTESTLOC)/PERMANENT


html: help_pages
	mkdir -p $@
	for f in $(MD_SRC); do pandoc -s --toc --template template.html --css template.css --toc-depth 2 -f gfm -t html $$f > $@/$${f%.*}.html ; done
	cp -r m.swim.*.html template.css img $@
	cd $@; \
	mv home.html index.html ;\
	find . -name "m.swim.*.html" -exec sed -i.backup -e "s:grass_logo.png:img/grass_logo.png:g" '{}' + ; \
	find . -name "*.html" -exec sed -i.backup -e "s/.md/.html/g" '{}' + ; \
	find . -name "*.html" -exec sed -i.backup -e "s/home//g" '{}' + ; \
	rm *.backup


publish: html
	rsync -av html/ $(URL)


clean:
	rm -rf html m.swim.*.html
