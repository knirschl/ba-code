find . -type f -name ba.\* -exec rm {} \;
find . -type f -name rf_\* -exec rm {} \;
find . -name "*.generax_pick*" -exec rename -vn ".generax_pick" "" {} \;