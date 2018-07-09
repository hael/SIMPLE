const {Image} = require('image-js');
const {Stack} = require('image-js');

var minrad = 10
var maxrad = 60

Image.load('test1.jpeg').then((image) => {
	var blur = image.grey()
	var average = new Uint32Array(blur.data.length / 4);
	
	for(var angle = 1; angle < 360; angle++){
		var rotated = blur.rotate(angle)
		var cropped = rotated.crop({x : ((rotated.width - blur.width) / 2), y : ((rotated.height - blur.height) / 2), width : (blur.width/2), height : (blur.height/2)})
		for(var pixel = 0; pixel < cropped.data.length; pixel++){
			average[pixel] += cropped.data[pixel]
		}
	}
	
	var cropped = blur.crop({width : (blur.width/2), height : (blur.height/2)})
	
	for(var pixel = 0; pixel < average.length; pixel++){
		cropped.data[pixel] = Math.floor(average[pixel]/360)
	}
	var rotated = cropped.rotate(180)
	var xmasscen = 0
	var ymasscen = 0
	var pixcount = 0
	for(var x = 0; x < rotated.width; x++){
		for(var y = 0; y < rotated.height; y++){
			var radius = Math.pow((Math.pow(x, 2) + Math.pow(y, 2)), 0.5)
			if(radius > (cropped.width - 10)){
				xmasscen += x * cropped.getPixelXY(x, y)
				ymasscen += y * cropped.getPixelXY(x, y)
				pixcount++
			}
		}
	}
	xmasscen /= pixcount
	ymasscen /= pixcount
	console.log(xmasscen, ymasscen)
	
	
	//var blur2 = blur.gaussianFilter({radius: 5, sigma: 1})
	//var blur2 = blur.resize({width:50, InterpolationAlgorithm : 'bilinear'})
	rotated.save('blur.jpg')

//	average.save("stl.jpg") ("huang" | "intermodes" | "isodata" | "li" | "maxentropy" | "mean" | "minerror" | "minimum" | "moments" | "otsu" | "percentile" | "renyientropy" | "shanbhag" | "triangle" | "yen") 
/*	
	queue.push([(blur.width / 2), (blur.height / 2 ), 0])
	
	while(queue.length > 0){
		if(pixels[queue[0][0]][queue[0][1]][0] != 0 && pixels[queue[0][0]][queue[0][1]][1] < maxrad){
			//blur.setPixelXY(queue[0][0], queue[0][1], [255])
			var pixelthreshold = pixels[queue[0][0]][queue[0][1]][0] - queue[0][2] -2
			pixels[queue[0][0]][queue[0][1]][0] = 0
			var repeat = true
			//North
			if(pixels[queue[0][0] - 1][queue[0][1]][0] > pixelthreshold){
				repeat = false
				queue.push([queue[0][0] - 1, queue[0][1], 0])
			}
			//NorthEast
			if(pixels[queue[0][0] - 1][queue[0][1] - 1][0] > pixelthreshold){
				repeat = false
				queue.push([queue[0][0] - 1, queue[0][1] - 1, 0])
			}
			//NorthWest
			if(pixels[queue[0][0] - 1][queue[0][1] + 1][0] > pixelthreshold){
				repeat = false
				queue.push([queue[0][0] - 1, queue[0][1] + 1, 0])
			}
			//East
			if(pixels[queue[0][0]][queue[0][1] - 1][0] > pixelthreshold){
				repeat = false
				queue.push([queue[0][0], queue[0][1] - 1, 0])
			}
			//West
			if(pixels[queue[0][0]][queue[0][1] + 1][0] > pixelthreshold){
				repeat = false
				queue.push([queue[0][0], queue[0][1] + 1, 0])
			}
			//South
			if(pixels[queue[0][0] + 1][queue[0][1]][0] > pixelthreshold){
				repeat = false
				queue.push([queue[0][0] + 1, queue[0][1], 0])
			}
			//SouthEast
			if(pixels[queue[0][0] + 1][queue[0][1] - 1][0] > pixelthreshold){
				repeat = false
				queue.push([queue[0][0] + 1, queue[0][1] - 1, 0])
			}
			//SouthWest
			if(pixels[queue[0][0] + 1][queue[0][1] + 1][0] > pixelthreshold){
				repeat = false
				queue.push([queue[0][0] + 1, queue[0][1] + 1, 0])
			}
			if(repeat){
				console.log("repeat")
				queue.push([queue[0][0], queue[0][1], queue[0][2]])
			}else{
				pixels[queue[0][0]][queue[0][1]][0] = 0
			}
		}
		queue.shift()
	}
	
	/*for(var x = 0; x < blur.width; x++){
		for(var y = 0; y < blur.height; y++){
			if(pixels[x][y][1] < minrad){
				blur.setPixelXY(x, y, [255])
			}
		}
	}

	blur2.save('blur.jpg');*/

})

