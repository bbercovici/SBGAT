
if cd doc/; then
	git rm -rf *
	git commit -m "purge before updates"
	doxygen ../Doxyfile
	git add *
	git commit -m "updates to documentation"
	git push origin gh-pages
else
	echo "Command failed"
fi
