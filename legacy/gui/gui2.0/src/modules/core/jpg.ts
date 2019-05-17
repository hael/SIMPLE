import * as sharp from 'sharp'

export default class JPG{

	constructor() {
	
	}
	
	toJPEG(arg){
		return sharp(arg['stackfile'])
			.resize(Number(arg['width']))
			.toBuffer()
			.then((image) => {
				return({image : image})
			})
	}
}
