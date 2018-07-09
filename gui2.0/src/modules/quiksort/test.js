const {Image} = require('image-js');
const {Roi} = require('image-js');
const {Stack} = require('image-js');

/*var particlediameter = 70
Image.load('test2.jpeg').then((image) => {
 // var blur = image.gaussianFilter({radius: 20, sigma: 5}).grey()
 var gaussianradius = 10
 var blur
 var mask
 var roiManager
 var rois
 var exit = true
 while(exit){  
	gaussianradius += 1
	blur = image.gaussianFilter({radius: gaussianradius, sigma: 3}).grey()
	mask = blur.mask({algorithm : "triangle"})
	roiManager = image.getRoiManager();
	roiManager.fromMask(mask);
	rois = roiManager.getRois({
	   positive: true,
	   negative: false,
	   minSurface: (particlediameter / 4) * (particlediameter / 4)
		});
		
	for (var ROI of rois){
		if(ROI.meanX > ((image.width / 2) - (particlediameter/2)) && ROI.meanX < ((image.width / 2) + (particlediameter/2)) && ROI.meanY >((image.width / 2) - (particlediameter/2)) && ROI.meanY < ((image.width / 2) + (particlediameter/2))){
			console.log(gaussianradius)
			exit = false
			break
		}
	}
 }
 
 
for (var ROI of rois){
	if(ROI.meanX > ((image.width / 2) - (particlediameter/2)) && ROI.meanX < ((image.width / 2) + (particlediameter/2)) && ROI.meanY >((image.width / 2) - (particlediameter/2)) && ROI.meanY < ((image.width / 2) + (particlediameter/2))){
		image.paintMasks(ROI.getMask(), {})	
	}
 image.save('img.jpg');
}
})*/

/*Image.load('test2.jpeg').then((image) => {
	image.grey().blurFilter().save("blackhat.jpg")
	var north = image
	var east = image.rotate(90)
	var south = image.rotate(180)
	var west = image.rotate(270)
	var average = Image.createFrom(image)
	
    for (let j = 0; j < image.data.length; j++) {
		average.data[j] = (north.data[j] + west.data[j] + south.data[j] + east.data[j]) / 4
    }
    average.save('average.jpg');
})*/


/*var particlediameter = 70
Image.load('test3.jpeg').then((image) => {
	var grey = image.grey()
	var blur = grey.gaussianFilter({radius: 7, sigma: 5})
	var mask = blur.mask({algorithm : "triangle"})
	//var open = blur.open({iterations : 1})
	var close = blur.close({iterations : 3, kernel : [[0,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0],
														[0,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0],
														[0,0,0,0,1,1,1,1,1,1,1,1,1,0,0,0,0],
														[0,0,0,1,1,1,1,1,1,1,1,1,1,1,0,0,0],
														[0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0],
														[0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0],
														[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
														[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
														[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
														[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
														[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
														[0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0],
														[0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0],
														[0,0,0,1,1,1,1,1,1,1,1,1,1,1,0,0,0],
														[0,0,0,0,1,1,1,1,1,1,1,1,1,0,0,0,0],
														[0,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0],
														[0,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0]
	]})
	var mask2 = close.mask({algorithm : "triangle"})
	mask.save('img.jpg');
	mask2.save('img2.jpg');

})*/

var minrad = 10
var maxrad = 40

Image.load('test.jpeg').then((image) => {
	var grey = image.grey()
	var blur = grey.gaussianFilter({radius: 5, sigma: 1})
	var radii = []
	var thresholdarray = []
	var highthreshold
	var lowthreshold
	
	for(var x = 0; x < blur.width; x++){
		for(var y = 0; y < blur.height; y++){
			var rad2 = Math.pow((x - (blur.width / 2 )), 2) + Math.pow((y - (blur.height / 2 )), 2)
			radii.push([x, y, Math.pow(rad2, 0.5)])
		}
	}
	
	for(var z = 0; z < radii.length; z++){
		if(radii[z][2] > maxrad){
			blur.setPixelXY(radii[z][0], radii[z][1], [0])
		}else if(radii[z][2] < minrad){
			blur.setPixelXY(radii[z][0], radii[z][1], [255])
		}else{
			thresholdarray.push(blur.getPixelXY(radii[z][0], radii[z][1])[0])
		}
	}
	
	thresholdarray.sort(function(a, b){return b-a})
	highthreshold = thresholdarray[minrad * minrad * 10]
	lowthreshold = thresholdarray[thresholdarray.length - (minrad * minrad * 10)]
	
	console.log(highthreshold, lowthreshold)
	
	for(var z = 0; z < radii.length; z++){
		if(radii[z][2] < maxrad && radii[z][2] > minrad){
			if(blur.getPixelXY(radii[z][0], radii[z][1]) > highthreshold){
				blur.setPixelXY(radii[z][0], radii[z][1], [255])
			}else if(blur.getPixelXY(radii[z][0], radii[z][1]) < lowthreshold){
				blur.setPixelXY(radii[z][0], radii[z][1], [0])
			}
		}
	}
	
	var close = blur.close({iterations : 3, kernel : [[0,0,0,1,1,1,0,0,0],
													  [0,0,1,1,1,1,1,0,0],
													  [0,1,1,1,1,1,1,1,0],
													  [1,1,1,1,1,1,1,1,1],
													  [1,1,1,1,1,1,1,1,1],
													  [1,1,1,1,1,1,1,1,1],
													  [0,1,1,1,1,1,1,1,0],
													  [0,0,1,1,1,1,1,0,0],
													  [0,0,0,1,1,1,0,0,0]
	]})
	close.save('blur.jpg');
	/*for(var pixel of grad){
		if(pixel[2] > 35){
			backgrad += pixel[3]
			backgradcount++
		}
	}
	backgrad /= backgradcount
	console.log(backgrad)
	
	backgrad = 0
	backgradcount = 0
	for(var pixel of grad){
		if(pixel[2] < 35){
			backgrad += pixel[3]
			backgradcount++
		}
	}
	backgrad /= backgradcount*/

	
	
})

