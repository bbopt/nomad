#!/bin/sh
here="`dirname \"$0\"`"
echo "cd-ing to $here"
cd "$here" || exit 1

cd html/;


mv tree.html test; sed s/'<\/body>'/'<br \/> <p><a href="http:\/\/www.gerad.ca\/nomad" target="_parent"><img src="..\/nomad.jpg" width="100" \/><\/a><\/p><p><a href="http:\/\/www.gerad.ca\/nomad" target="_parent">http:\/\/www.gerad.ca\/nomad<\/a><\/p><\/body>'/g test > tree.html;

rm test;

mv doxygen.css test; sed s/'background: white'/'background-image:url("..\/NOMAD_background.jpg"); background-repeat:repeat-y;'/g test > doxygen.css;

rm test;


