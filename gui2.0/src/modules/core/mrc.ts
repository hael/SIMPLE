import * as fs from 'fs'
import * as binaryfile from 'binary-file'
import * as sharp from 'sharp'

export default class MRC{

	constructor() {

	}
	
	readHeader(file) {
		var header = {}
		const datain = new binaryfile(file, 'r', true)
		return datain.open()
			.then(function(){
				return datain.read(1024, 0)
			})
			.then(function(data){
				header['data'] = data
				return datain.close()
			})
			.then(function(){
				header['nx'] = header['data'].readInt32LE(0)
				header['ny'] = header['data'].readInt32LE(4)
				header['nz'] = header['data'].readInt32LE(8)
				header['mode'] = header['data'].readInt32LE(12)
				header['nsymbt'] = header['data'].readInt32LE(92)
				return (header)
			})
			.catch(function (err) {
				console.log(`There was an error: ${err}`);
			})
	}
	
	readData(file, frame=0) {
		var mrcfile
		var datafile
		return this.readHeader(file)
			.then((data) => {
				mrcfile = data
				datafile = new binaryfile(file, 'r', true)
				return datafile.open()
			})
			.then(() => {
				return datafile.seek(1024 + mrcfile['nsymbt'] + Number(frame) * mrcfile['nx'] * mrcfile['ny'] * 4)
			})
			.then(() => {
				if(mrcfile.mode == 2){
					return datafile.read(mrcfile['nx'] * mrcfile['ny'] * 4, 1024 + mrcfile['nsymbt'] + Number(frame) * mrcfile['nx'] * mrcfile['ny'] * 4)
					
				}
			})
			.then((val) => {
				var elements = mrcfile.nx * mrcfile.ny
				mrcfile['data'] = []
				mrcfile['data'].length = elements
				for(var i = 0; i < elements; i++){
					mrcfile['data'][i] = val.readFloatLE(i * 4)
				}
				return datafile.close()
			})
			.then(() => {
				return(mrcfile)
			})
	}
		
	
	toJPEG(arg){
		return this.readData(arg['stackfile'], arg['frame'])
			.then((mrcfile) => {
				var mean = 0
				var variance = 0
				var tmp
				
				for(var i = 0; i < mrcfile.nx * mrcfile.ny; i++){
					tmp = mean
					mean = tmp + mrcfile['data'][i]
				}
				mean /= mrcfile.nx * mrcfile.ny

				for(var i = 0; i < mrcfile.nx * mrcfile.ny; i++){
					tmp = variance
					variance = tmp + Math.pow((mrcfile['data'][i] - mean), 2)
				}
				variance /= mrcfile.nx * mrcfile.ny
				var sd = Math.sqrt(variance)
				var pixels = []
				pixels.length = mrcfile.nx * mrcfile.ny
				var norm
				
				for(var i = 0; i < mrcfile.nx * mrcfile.ny; i++){
					norm = (mrcfile['data'][i] - mean) / sd
					if(norm > 10){
						pixels[i] = 255
					} else if (norm < -10){
						pixels[i] = 0
					} else {
						pixels[i] = Math.round((norm * 12.8) + 128)
					}
				}

				var pixbuf =  Buffer.from(pixels)

				return sharp(pixbuf, { raw : {
					width : mrcfile.nx,
					height : mrcfile.ny,
					channels : 1
					}
				}).resize(Number(arg['width']))
				  .jpeg()
				  .toBuffer()
			})
			.then(function (image) {
				return({image : image})
			})
	}
}
