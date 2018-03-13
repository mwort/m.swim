
# A Makefile to publish and convert the documentation
include ../test/Makefile

GITHUBWIKI = git@github.com:mwort/m.swim.wiki.git
# GITLABWIKI = git@gitlab.pik-potsdam.de:wortmann/m.swim.wiki.git
# git remote add gitlab $(GITLABWIKI); \


help_pages: m.swim.%.html

m.swim.%.html:
	GRASS_BATCH_JOB=./build_help_pages.grass $(GRASS) $(GRASSTESTLOC)/PERMANENT


wiki: help_pages
	git init wiki
	for f in m.swim.*.html; do pandoc -f html -t markdown_github $$f > wiki/$${f%.*}.md ; done
	cd wiki; \
	cp -r ../*.md ../img . ; \
	mv README.md _Sidebar.md ; \
	find . -name "m.swim.*.md" -exec sed -i.backup -e "s:grass_logo.png:img/grass_logo.png:g" '{}' + ; \
	find . -name "*.md" -exec sed -i.backup -e "s/.md//g" '{}' + ; \
	find . -name "*.md" -exec sed -i.backup -e "s/home//g" '{}' + ; \
	rm *.backup


publish: wiki
	cd wiki; \
	git remote add github $(GITHUBWIKI); \
	git add *; git commit -m 'Wiki files.' ; \
	git push github master --force


clean:
	rm -rf wiki m.swim.*.html
