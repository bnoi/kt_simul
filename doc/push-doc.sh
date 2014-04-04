#!/usr/bin/env sh

cd ../
git fetch upstream gh-pages:gh-pages
git checkout gh-pages
git rm -fr .
git checkout master .
cd doc/
rm -fr source/generated/
sphinx-autogen source/*.rst
make html
mv build ../
cd ../
rm -fr $(find . -type f -maxdepth 1 | grep -v build |grep -v \.git)
rm -fr .gitignore
mv build/html/* .
rm -fr build/
touch .nojekyll
git add .
git ci -am "Generated gh-pages for `git log master -1 --pretty=short --abbrev-commit`" && git push upstream gh-pages && git checkout master
cd ../
