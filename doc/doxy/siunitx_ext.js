 
// MathJax.Hub.Config({
//   extensions: ["tex2jax.js","[siunitx]/siunitx.js"],
//   jax: ["input/TeX","output/HTML-CSS"],
//   tex2jax: {inlineMath: [["$","$"],["\\(","\\)"]]},
//   TeX: {extensions: ["AMSmath.js","AMSsymbols.js", "sinuitx.js"]}
// });
// MathJax.Ajax.config.path['siunitx']  = 'http://rawgit.com/burnpanck/MathJax-siunitx/master/';
  MathJax.Hub.Config({
   extensions: ["tex2jax.js","[Contrib]/siunitx/siunitx.js"],
   jax: ["input/TeX","output/HTML-CSS"],
   tex2jax: {inlineMath: [["$","$"],["\\(","\\)"]]},
   TeX: {extensions: ["AMSmath.js","AMSsymbols.js", "sinuitx.js"]}
 });
 MathJax.Ajax.config.path["Contrib"] = "https://cdn.mathjax.org/mathjax/contrib";