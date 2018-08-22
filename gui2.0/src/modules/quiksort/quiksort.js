const {Image} = require('image-js');
const {Stack} = require('image-js');

var minrad = 10
var maxrad = 60

Image.load('test1.jpeg').then((image) => {
	var blur = image.grey()
	for(var x = 0; x < 5; x++){
			for(var y = 0; y < 5; y++){
				var cropped = blur.crop({
					x:x,
					y:y,
					width: (image.width - 5),
					height: (image.height - 5)
				});	
				ring(cropped)
			}
	}
})

function ring (blur) {	
	var histogram = []
	var histogrampop = []
	for(var i = 0; i < Math.floor((blur.width / 2) / 2); i++){
		histogram.push(0)
		histogrampop.push(0)
	}
	for(var x = 0; x < blur.width; x++){
		for(var y = 0; y < blur.height; y++){
			var rad = Math.pow(Math.pow((x - (blur.width / 2 )), 2) + Math.pow((y - (blur.height / 2 )), 2), 0.5)
			histogram[Math.floor(rad / 2)] += Number(blur.getPixelXY(x,y))
			histogrampop[Math.floor(rad / 2)]++
		}
	}
	for(var i = 0; i < Math.floor((blur.width / 2) / 2); i++){
		histogram[i] /= histogrampop[i]	
	}
	var minmean = 0
	var minmeanpop = 0
	for(var i = 0; i <= Math.floor(minrad / 2); i++){
		minmean += histogram[i]
		minmeanpop++
	}
	minmean /= minmeanpop
	
	var outmean = 0
	var outmeanpop = 0
	for(var i = Math.floor(maxrad / 2);i < Math.floor((blur.width / 2) / 2); i++){
		outmean += histogram[i]
		outmeanpop++
	}
	outmean /= outmeanpop
	threshold = (minmean + outmean) / 2
	
	console.log(minmean, outmean, threshold)

	for(var i = Math.floor(minrad / 2) - 3; i < Math.floor((blur.width / 2) / 2); i++){
		if(histogram[i]  < threshold && histogram[i+1]  < threshold && histogram[i+2]  < threshold){
			console.log(i * 2)
			break
		}
	}	
}
