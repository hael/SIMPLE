declare var d3
declare var StackBlur
declare var define
class ManualPickingWidget {

	private widgetpopup
	private radius
	private boxsize
	private coords
	private canvas
	private particlecanvas
	private particlectx
	private ctx
	private img
	private micrographName
	private allCoords
	private ratio
	private microlist
	private selectedIndex
	
	constructor(){
		this.widgetpopup = document.getElementById('widgetpopup')
	}
	
	refresh(){
		var request = {
			mod : "manualpick",
			fnc : "getManualPickWidget",
			arg : {}
		}
		request['arg']['micfile'] = (<HTMLInputElement>document.getElementById('inputpath')).value
		return postAjaxPromise(request)
			.then((response) => response.json())
			.then ((json) => {
				this.widgetpopup.innerHTML = json.html
				this.widgetpopup.className = "guinierwidgetpopup"
				this.coords = []
				this.allCoords = []
				this.canvas = document.getElementById("pickingcanvas")
				this.radius = document.getElementById("manualpickradius")
				this.boxsize = document.getElementById("manualpickboxsize")
				this.ctx = this.canvas.getContext("2d")
				this.particlecanvas = document.getElementById("particlecanvas")
				this.particlectx = this.particlecanvas.getContext("2d")
				this.img = new Image();
				this.img.addEventListener('load',() => {
					this.ctx.drawImage(this.img,0,0)
					this.plotCoordinates()
					}, false)
				this.microlist = <HTMLInputElement>document.getElementById("micrographlist")
				this.microlist.options[0].click()
			})
	}
	
	view(element){
		this.refresh()
	}
	
	hide(){
		this.widgetpopup.innerHTML = ""
		this.widgetpopup.className = ""
	}
	
	public plotCoordinates() {
		var r = this.radius.value
		var l = this.boxsize.value
		this.ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);
		this.ctx.drawImage(this.img, 0, 0);
		
		this.ctx.lineWidth = 4;
		this.ctx.lineJoin = "round";
		this.ctx.strokeStyle = "magenta";
		
		var plotR = document.getElementById('showbox')
		if(plotR['checked'] == true){
			for(var circle of this.coords){
				this.ctx.strokeRect(circle.xcoord-(l/2), circle.ycoord-(l/2), l, l);
			}
		}
		
		var plotC = document.getElementById('showcircle')
		if(plotC['checked'] == true){
			for(var circle of this.coords){
				this.ctx.beginPath();
				this.ctx.arc(circle.xcoord, circle.ycoord, r, 0, 2*Math.PI);
				this.ctx.stroke();
			}
		}
		
		
	}
	
	public clickHandler2(event){
	//	!function(e,t){"object"==typeof exports&&"undefined"!=typeof module?t(exports):"function"==typeof define&&define.amd?define(["exports"],t):t(e.StackBlur=e.StackBlur||{})}(this,function(e){"use strict";function t(e,t){if(!(isNaN(t)||t<1)){t|=0;var o,i,u,x,a,d,l,c,s,p,h,v=e.data,y=e.width,w=e.height,b=t+t+1,j=y-1,k=w-1,m=t+1,B=m*(m+1)/2,N=new r,S=N;for(u=1;u<b;u++)if(S=S.next=new r,u==m)var _=S;S.next=N;var g=null,M=null;l=d=0;var O=n[t],P=f[t];for(i=0;i<w;i++){for(p=c=0,s=m*(h=v[d]),c+=B*h,S=N,u=0;u<m;u++)S.r=h,S=S.next;for(u=1;u<m;u++)x=d+((j<u?j:u)<<2),c+=(S.r=h=v[x])*(m-u),p+=h,S=S.next;for(g=N,M=_,o=0;o<y;o++)v[d]=c*O>>P,c-=s,s-=g.r,x=l+((x=o+t+1)<j?x:j)<<2,p+=g.r=v[x],c+=p,g=g.next,s+=h=M.r,p-=h,M=M.next,d+=4;l+=y}for(o=0;o<y;o++){for(p=c=0,d=o<<2,s=m*(h=v[d]),c+=B*h,S=N,u=0;u<m;u++)S.r=h,S=S.next;for(a=y,u=1;u<=t;u++)d=a+o<<2,c+=(S.r=h=v[d])*(m-u),p+=h,S=S.next,u<k&&(a+=y);for(d=o,g=N,M=_,i=0;i<w;i++)x=d<<2,v[x]=c*O>>P,c-=s,s-=g.r,x=o+((x=i+m)<k?x:k)*y<<2,c+=p+=g.r=v[x],g=g.next,s+=h=M.r,p-=h,M=M.next,d+=y}return e}}function r(){this.r=0,this.next=null}var n=[512,512,456,512,328,456,335,512,405,328,271,456,388,335,292,512,454,405,364,328,298,271,496,456,420,388,360,335,312,292,273,512,482,454,428,405,383,364,345,328,312,298,284,271,259,496,475,456,437,420,404,388,374,360,347,335,323,312,302,292,282,273,265,512,497,482,468,454,441,428,417,405,394,383,373,364,354,345,337,328,320,312,305,298,291,284,278,271,265,259,507,496,485,475,465,456,446,437,428,420,412,404,396,388,381,374,367,360,354,347,341,335,329,323,318,312,307,302,297,292,287,282,278,273,269,265,261,512,505,497,489,482,475,468,461,454,447,441,435,428,422,417,411,405,399,394,389,383,378,373,368,364,359,354,350,345,341,337,332,328,324,320,316,312,309,305,301,298,294,291,287,284,281,278,274,271,268,265,262,259,257,507,501,496,491,485,480,475,470,465,460,456,451,446,442,437,433,428,424,420,416,412,408,404,400,396,392,388,385,381,377,374,370,367,363,360,357,354,350,347,344,341,338,335,332,329,326,323,320,318,315,312,310,307,304,302,299,297,294,292,289,287,285,282,280,278,275,273,271,269,267,265,263,261,259],f=[9,11,12,13,13,14,14,15,15,15,15,16,16,16,16,17,17,17,17,17,17,17,18,18,18,18,18,18,18,18,18,19,19,19,19,19,19,19,19,19,19,19,19,19,19,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24];e.R=t,Object.defineProperty(e,"__esModule",{value:!0})});
		var x = event.offsetX;
		var y = event.offsetY;
		var xclick = x
		var yclick = y
		var r = this.radius.value
		var l = this.boxsize.value
		
		x = this.ratio * x
		y = this.ratio * y
		
		var arrayLength = this.coords.length;
		var totalCoords = this.allCoords.length;

		if (r > 0) {
			for (var i = arrayLength - 1; i > -1; i--) {	
				var coordpair = this.coords[i];
				if ((Math.pow((x - coordpair.xcoord), 2) + Math.pow((y - coordpair.ycoord), 2)) < Math.pow(r, 2)) {
					this.coords.splice(i, 1);
				}
			}
		}
		if (r > 0) {
			for (var i = totalCoords - 1; i > -1; i--) {	
				coordpair = {"xcoord":this.allCoords[i].xcoord, "ycoord":this.allCoords[i].ycoord}
				if ((Math.pow((x - coordpair.xcoord), 2) + Math.pow((y - coordpair.ycoord), 2)) < Math.pow(r, 2) && this.allCoords[i].name == this.micrographName) {
					this.allCoords.splice(i, 1);
				}
			}
		}
	
		if (l > 0) {
			var newArrayLength = this.coords.length;
			for (var i = newArrayLength - 1; i >- 1; i--) {	
				coordpair = this.coords[i];
				if (x - (l/2) < coordpair.xcoord && x + (l/2) > coordpair.xcoord && y - (l/2) < coordpair.ycoord && y + (l/2) > coordpair.ycoord) {
					this.coords.splice(i, 1);
				}
			}
		}
		if (l > 0) {
			var newArrayLength = this.allCoords.length;
			for (var i = newArrayLength - 1; i >- 1; i--) {	
				coordpair = {"xcoord":this.allCoords[i].xcoord, "ycoord":this.allCoords[i].ycoord}
				if (x - (l/2) < coordpair.xcoord && x + (l/2) > coordpair.xcoord && y - (l/2) < coordpair.ycoord && y + (l/2) > coordpair.ycoord && this.allCoords[i].name == this.micrographName) {
					this.allCoords.splice(i, 1);
				}
			}
		}
		
		
		
		if (arrayLength == this.coords.length) {
			this.ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);
			this.ctx.drawImage(this.img,0,0);
		
			//iteration 1
			var newImg = this.ctx.getImageData( x - (l/2) , y - (l/2) , l , l);
			this.particlecanvas.width = l
			this.particlecanvas.height = l
			this.particlectx.clearRect(0, 0, l, l)
			this.particlectx.putImageData(newImg, 0, 0);
		
			var cont = document.getElementById('showcontours')
			if(cont['checked'] == true){
				var values = new Array( l * l )
				var contours = d3.contours().size([l, l])
				var projection = d3.geoIdentity().scale(1)
				var path = d3.geoPath(projection, this.particlectx)
				StackBlur.R(newImg, 7)
				for(var i = 0; i < (l * l); i++){
					values[i] = 1 - newImg.data[(i << 2)] / 255;
				}
				var threshold
				var calccontours
		
				for(threshold = 50; threshold < 200; threshold += 1){
					calccontours = contours.thresholds([(threshold) / 255])(values)	
					var area = path.area(calccontours[0])
					if(area < (l * l / 3)){
						console.log(area)
						console.log(threshold)
						break
					}
				}
		
				var fakecontour = contours.thresholds([(threshold) / 255])(values)
				fakecontour[0]['coordinates'] = []
				for(var i = 0; i < calccontours[0]['coordinates'].length; i++){
					fakecontour[0]['coordinates'] = []
					fakecontour[0]['coordinates'][0] = calccontours[0]['coordinates'][i]
					var area = path.area(fakecontour[0])
					var centroid = path.centroid(fakecontour[0])
					var radius = 80
					if(area < 1000){
						calccontours[0]['coordinates'][i] = []	
						}else if ( Math.pow((centroid[0]-(l/2)) , 2) + Math.pow((centroid[1]-(l/2)), 2) > Math.pow(radius, 2)) { 
						calccontours[0]['coordinates'][i] = []
					}
					
				}
				var centroid = path.centroid(calccontours[0])
				console.log("Coordinates: " + centroid[0] + " " + centroid[1])
				
				if (isNaN(centroid[0])){
					x = xclick
					y = yclick
					this.particlectx.clearRect(0, 0, l, l)
					document.getElementById("partFound").innerHTML ="No particle found on last click."
					this.plotCoordinates()
					return
				} else {
					document.getElementById("partFound").innerHTML =""
				}
				
				
			//iteration 2
				x += centroid[0] - l/2
				y += centroid[1] - l/2
				
				
				newImg = this.ctx.getImageData( x - (l/2) , y - (l/2) , l , l);
				this.particlecanvas.width = l
				this.particlecanvas.height = l
				this.particlectx.clearRect(0, 0, l, l)
				this.particlectx.putImageData(newImg, 0, 0);
				var values = new Array( l * l )
				var contours = d3.contours().size([l, l])
				var projection = d3.geoIdentity().scale(1)
				var path = d3.geoPath(projection, this.particlectx)
				StackBlur.R(newImg, 7)
				for(var i = 0; i < (l * l); i++){
					values[i] = 1 - newImg.data[(i << 2)] / 255;
				}
				var threshold
				var calccontours
		
				for(threshold = 50; threshold < 200; threshold += 1){
					calccontours = contours.thresholds([(threshold) / 255])(values)	
					var area = path.area(calccontours[0])
					if(area < (l * l / 3)){
						break
					}
				}
		
				var fakecontour = contours.thresholds([(threshold) / 255])(values)
				fakecontour[0]['coordinates'] = []
				for(var i = 0; i < calccontours[0]['coordinates'].length; i++){
					fakecontour[0]['coordinates'] = []
					fakecontour[0]['coordinates'][0] = calccontours[0]['coordinates'][i]
					var area = path.area(fakecontour[0])
					var centroid = path.centroid(fakecontour[0])
					var radius = 60
					if(area < 1000){
						calccontours[0]['coordinates'][i] = []	
					}else if ( Math.pow((centroid[0]-(l/2)) , 2) + Math.pow((centroid[1]-(l/2)), 2) > Math.pow(radius, 2)) { 
						calccontours[0]['coordinates'][i] = []
					}
				}
				
				this.particlectx.lineWidth = 2;
				this.particlectx.lineJoin = "round";
				this.particlectx.strokeStyle = "magenta";
				this.particlectx.beginPath();
				path(calccontours[0]);
				var centroid = path.centroid(calccontours[0])
				console.log("Coordinates: " + centroid[0])
				console.log(centroid[1])
			
			
				this.particlectx.stroke();
			} else {
				document.getElementById("partFound").innerHTML =""
			}
			
			this.coords.push({"xcoord":x, "ycoord":y})	
			this.allCoords.push({"name":this.micrographName, "xcoord":(x-(l/2)), "ycoord":(y-(l/2)), "boxsize":l})
			
			this.plotCoordinates()
		
		} else {
			this.particlectx.clearRect(0, 0, l, l)
			this.plotCoordinates()
		}
		
		document.getElementById("partNo").innerHTML ="Number of particles: " + this.coords.length + "/" + this.allCoords.length
		
	}
	
	
	
	public changeMicrograph(path, xdim, ydim, micrographname){
		
		//canvas dimensions defined here
		this.canvas.width = xdim
		this.canvas.height = ydim
		this.micrographName = micrographname
		this.img.src = "/image?stackfile=" + path + "&frame=0&width=" + xdim
		
		//find previous boxes
		this.coords = []
		for(var asdf of this.allCoords){
			if(asdf.name == micrographname){
				this.coords.push({ "xcoord":(asdf.xcoord + (asdf.boxsize/2)), "ycoord":(asdf.ycoord + (asdf.boxsize/2)) })
			}
		}
		document.getElementById("partNo").innerHTML ="Number of particles: " + this.coords.length + "/" + this.allCoords.length
		this.ratio = xdim / this.canvas.scrollWidth
		var l = this.boxsize.value
		this.particlectx.clearRect(0, 0, l, l)
		
	}


	public prevMicrograph(event){
		this.microlist = <HTMLInputElement>document.getElementById("micrographlist")
		this.selectedIndex = this.microlist.selectedIndex
		
		if( this.selectedIndex == 0 ){
		} else {
			this.microlist.selectedIndex = this.selectedIndex - 1
			this.microlist.options[this.selectedIndex - 1].click()
		}
		
	}

	public nextMicrograph(event){
		this.microlist = <HTMLInputElement>document.getElementById("micrographlist")
		this.selectedIndex = this.microlist.selectedIndex
		
		if( this.selectedIndex == this.microlist.length-1 ){
			console.log("NOPE")
		} else {
			this.microlist.selectedIndex = this.selectedIndex + 1
			this.microlist.options[this.selectedIndex + 1].click()
		}
	
	}

	public donePicking(event){
		console.log(this.allCoords);
		(<HTMLInputElement>document.getElementById('keypick')).value = JSON.stringify(this.allCoords)
		this.hide()
	}
	
}

var manualpickingwidget
window.addEventListener('load', () => {manualpickingwidget = new ManualPickingWidget()})
!function(e,t){"object"==typeof exports&&"undefined"!=typeof module?t(exports):"function"==typeof define&&define.amd?define(["exports"],t):t(e.StackBlur=e.StackBlur||{})}(this,function(e){"use strict";function t(e,t){if(!(isNaN(t)||t<1)){t|=0;var o,i,u,x,a,d,l,c,s,p,h,v=e.data,y=e.width,w=e.height,b=t+t+1,j=y-1,k=w-1,m=t+1,B=m*(m+1)/2,N=new r,S=N;for(u=1;u<b;u++)if(S=S.next=new r,u==m)var _=S;S.next=N;var g=null,M=null;l=d=0;var O=n[t],P=f[t];for(i=0;i<w;i++){for(p=c=0,s=m*(h=v[d]),c+=B*h,S=N,u=0;u<m;u++)S.r=h,S=S.next;for(u=1;u<m;u++)x=d+((j<u?j:u)<<2),c+=(S.r=h=v[x])*(m-u),p+=h,S=S.next;for(g=N,M=_,o=0;o<y;o++)v[d]=c*O>>P,c-=s,s-=g.r,x=l+((x=o+t+1)<j?x:j)<<2,p+=g.r=v[x],c+=p,g=g.next,s+=h=M.r,p-=h,M=M.next,d+=4;l+=y}for(o=0;o<y;o++){for(p=c=0,d=o<<2,s=m*(h=v[d]),c+=B*h,S=N,u=0;u<m;u++)S.r=h,S=S.next;for(a=y,u=1;u<=t;u++)d=a+o<<2,c+=(S.r=h=v[d])*(m-u),p+=h,S=S.next,u<k&&(a+=y);for(d=o,g=N,M=_,i=0;i<w;i++)x=d<<2,v[x]=c*O>>P,c-=s,s-=g.r,x=o+((x=i+m)<k?x:k)*y<<2,c+=p+=g.r=v[x],g=g.next,s+=h=M.r,p-=h,M=M.next,d+=y}return e}}function r(){this.r=0,this.next=null}var n=[512,512,456,512,328,456,335,512,405,328,271,456,388,335,292,512,454,405,364,328,298,271,496,456,420,388,360,335,312,292,273,512,482,454,428,405,383,364,345,328,312,298,284,271,259,496,475,456,437,420,404,388,374,360,347,335,323,312,302,292,282,273,265,512,497,482,468,454,441,428,417,405,394,383,373,364,354,345,337,328,320,312,305,298,291,284,278,271,265,259,507,496,485,475,465,456,446,437,428,420,412,404,396,388,381,374,367,360,354,347,341,335,329,323,318,312,307,302,297,292,287,282,278,273,269,265,261,512,505,497,489,482,475,468,461,454,447,441,435,428,422,417,411,405,399,394,389,383,378,373,368,364,359,354,350,345,341,337,332,328,324,320,316,312,309,305,301,298,294,291,287,284,281,278,274,271,268,265,262,259,257,507,501,496,491,485,480,475,470,465,460,456,451,446,442,437,433,428,424,420,416,412,408,404,400,396,392,388,385,381,377,374,370,367,363,360,357,354,350,347,344,341,338,335,332,329,326,323,320,318,315,312,310,307,304,302,299,297,294,292,289,287,285,282,280,278,275,273,271,269,267,265,263,261,259],f=[9,11,12,13,13,14,14,15,15,15,15,16,16,16,16,17,17,17,17,17,17,17,18,18,18,18,18,18,18,18,18,19,19,19,19,19,19,19,19,19,19,19,19,19,19,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24];e.R=t,Object.defineProperty(e,"__esModule",{value:!0})});
