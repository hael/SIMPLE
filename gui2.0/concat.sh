for cssfolder in `find src -name "css"`; do
    cat $cssfolder/*.css >> dist/public/css/index.css
done

for jsfolder in `find dist -name "js"`; do
    cat $jsfolder/*.js >> dist/public/js/index.js
done

for imgfolder in `find src -name "img"`; do
    cp $imgfolder/* dist/public/img
done

for viewsfolder in `find src -name "views"`; do
    cp $viewsfolder/*.pug dist/views/
done

cat external/LiteMol/js/LiteMol-plugin.min.js >> dist/public/js/index.js
cat external/LiteMol/css/LiteMol-plugin.min.css >> dist/public/css/index.css
cat external/Chartist/dist/chartist.min.js >> dist/public/js/index.js
cat external/Chartist/dist/chartist.min.css >> dist/public/css/index.css
cp -r src/themes/* dist/public/css/
