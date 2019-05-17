const cannyEdgeDetector = require("canny-edge-detector");
const image = require('image-js');
 
image.Image.load('test.jpeg').then((img) => {
  const grey = img.grey();
  const edge = cannyEdgeDetector(grey, {gaussianBlur: 10, lowThreshold: 5, highThreshold:100});
  return edge.save('edge.jpeg');
})
