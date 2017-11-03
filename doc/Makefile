
# A Makefile to publish and convert the documentation

GITHUBWIKI = https://github.com/mwort/m.swim.wiki.git
# GITLABWIKI = git@gitlab.pik-potsdam.de:wortmann/m.swim.wiki.git
#	git remote set-url --add --push origin $(GITLABWIKI); \


wiki:
	git init wiki
	cd wiki; \
	git remote add origin $(GITHUBWIKI); \
	cp -r ../*.md ../img . ; \
	mv README.md _Sidebar.md ; \
	find . -name "*.md" -exec sed -i '' -e "s/.md//g" '{}' + ; \
	git add *; git commit -m 'Wiki files.' ; \
	git push origin master --force

clean:
	rm -rf wiki
