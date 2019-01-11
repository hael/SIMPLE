echo "" > public/css/index.css
echo "" > public/js/index.js

for cssfolder in `find server -name "css"`; do
    cat $cssfolder/*.css >> public/css/index.css
done

for jsfolder in `find server -name "js"`; do
    cat $jsfolder/*.js >> public/js/index.js
done

for imgfolder in `find server -name "img"`; do
    cp $imgfolder/* public/img/
done

#cat external/LiteMol/js/LiteMol-plugin.min.js >> public/js/index.js
#cat external/LiteMol/css/LiteMol-plugin.min.css >> public/css/index.css
#cat external/Chartist/dist/chartist.min.js >> public/js/index.js
#cat external/Chartist/dist/chartist.min.css >> public/css/index.css
