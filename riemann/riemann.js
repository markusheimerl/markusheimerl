/*
* Copyright 2021 Markus Heimerl, OTH Regensburg
* Licensed under CC BY-NC 4.0
*
* ANY USE OF THIS SOFTWARE MUST COMPLY WITH THE
* CREATIVE COMMONS ATTRIBUTION-NONCOMMERCIAL 4.0 INTERNATIONAL LICENSE
* AGREEMENTS
*/

// --- COMPLEX NUMBERS ---
class Complex {
	constructor(real, imaginary) {
		this.real = real;
		this.imaginary = imaginary;
	}
}

function trans(a) {
	return (typeof (a) == "number" || typeof (a) == "float") ? new Complex(a, 0) : a;
}

function sub_complex(a, b) {
	a = trans(a);
	b = trans(b);
	return new Complex(a.real - b.real, a.imaginary - b.imaginary);
}

function add_complex(a, b) {
	a = trans(a);
	b = trans(b);
	return new Complex(a.real + b.real, a.imaginary + b.imaginary);
}

function mult_complex(a, b) {
	a = trans(a);
	b = trans(b);
	return new Complex(a.real * b.real + a.imaginary * b.imaginary * (-1), a.real * b.imaginary + a.imaginary * b.real);
}

function div_complex(a, b) {
	a = trans(a);
	b = trans(b);

	var denom = (b.real * b.real + b.imaginary * b.imaginary);
	return new Complex((a.real * b.real + a.imaginary * b.imaginary) / denom, (a.imaginary * b.real - a.real * b.imaginary) / denom);
}

function pow_complex(b, e) {
	var real = Math.pow(b, e.real) * Math.cos(e.imaginary * Math.log(b));
	var imaginary = Math.pow(b, e.real) * Math.sin(e.imaginary * Math.log(b));
	return new Complex(real, imaginary);
}
// --- ---

function clamp(val, min, max) {
	return val > max ? max : val < min ? min : val;
}

// --- CartesianCoordinateSystem ---
class CartesianCoordinateSystem {
	/**
	 * [CartesianCoordinateSystem Creates a scrollable, zoomable, blank standard cartesian coordinate system]
	 * @param  {string} canvasname [The ID of the canvas on which to draw]
	 * @param  {[number]} scalex [Length of one unit in x direction]
	 * @param  {[number]} scaley [Length of one unit in y direction]
	 * @param  {[string]} xkey [label for x axis]
	 * @param  {[string]} ykey [label for y axis]
	 * @param  {[number]} xdense [density for value labels along x axis, 0 to turn off]
	 * @param  {[number]} ydense [density for value labels along y axis, 0 to turn off]
	 * @param  {[number]} translatex [move viewport this many pixels into x direction]
	 * @param  {[number]} translatey [move viewport this many pixels into y direction]
	 */
	constructor(canvasname, scalex, scaley, xdense, ydense, translatex, translatey, cartfunc, guifunc, sliderinpfunc, slidposxstar, slidposystar, sliderlen, totslidnum) {

		this.jquerycanvas = $("#" + canvasname);
		this.canvas = this.jquerycanvas[0];
		this.canvas.width = this.jquerycanvas.parent().width() - this.jquerycanvas.parent().width() * 0.0;
		this.ctx = this.canvas.getContext('2d');
		this.width = this.canvas.width;
		this.height = this.canvas.height;
		this.scalex = scalex;
		this.scaley = scaley;
		this.xkey = "real";
		this.ykey = "imaginary";
		this.xdense = xdense;
		this.ydense = ydense;
		this.translatex = translatex;
		this.translatey = translatey;
		this.mousedown = false;
		this.mdowncoords = [0, 0];
		this.drawCartFunction = cartfunc;
		this.drawGui = guifunc;
		this.handleSliderInput = sliderinpfunc;
		this.running = false;
		this.n_max = 150;
		this.data = {};
		this.current_n = 1;
		this.sliderposxStart = slidposxstar;
		this.sliderposyStart = slidposystar;
		this.sliderlength = sliderlen;
		this.totalSliderNum = totslidnum;
		this.slideryspace = 30;
		this.pausebtnx = 5;
		this.pausebtny = this.height - 55;
		this.pausebtnwidth = 50;
		this.pausebtnheight = 50;
		this.navbtnwidth = this.pausebtnwidth;
		this.navbtnheight = this.navbtnwidth;
		this.navx = this.width - this.navbtnwidth - 5;
		this.navy = this.height - 105;

		$("#" + canvasname)[0].addEventListener("mousedown", function (event) {
			this.mousedown = true;
			this.mdowncoords[0] = event.clientX - this.canvas.getBoundingClientRect().left;
			this.mdowncoords[1] = event.clientY - this.canvas.getBoundingClientRect().top;
		}.bind(this));

		$("#" + canvasname)[0].addEventListener("mouseup", function (event) {
			this.mousedown = false;
		}.bind(this));

		$("#" + canvasname)[0].addEventListener("click", function (event) {
			if (this.checkslider()) {
				this.handleSliderInput(event, this);
			} else if (this.checkPauseButton()) {
				this.handlePause();
			} else if (this.checkZoomSlider()) {
				this.handleZoom(event);
			} else {
				this.checkAndHandleNavigation();
			}
		}.bind(this));

		$("#" + canvasname)[0].addEventListener("mousemove", function (event) {
			if (this.mousedown == true) {
				if (this.checkslider()) {
					this.handleSliderInput(event, this);
				} else if (this.checkZoomSlider()) {
					this.handleZoom(event);
				}
			}
		}.bind(this));
	}

	checkslider() {
		return this.mdowncoords[0] > this.sliderposxStart && this.mdowncoords[0] < this.sliderposxStart + this.sliderlength &&
			this.mdowncoords[1] > this.sliderposyStart - this.slideryspace / 2 && this.mdowncoords[1] < this.sliderposyStart + this.slideryspace * this.totalSliderNum;
	}

	checkPauseButton() {
		return this.mdowncoords[0] > this.pausebtnx && this.mdowncoords[0] < this.pausebtnx + this.pausebtnwidth &&
			this.mdowncoords[1] > this.pausebtny && this.mdowncoords[1] < this.pausebtny + this.pausebtnheight;
	}

	checkZoomSlider() {
		return this.mdowncoords[0] > this.pausebtnx + this.pausebtnwidth + 10 && this.mdowncoords[0] < this.pausebtnx + this.pausebtnwidth + 10 + 140 &&
			this.mdowncoords[1] > this.pausebtny + this.pausebtnheight / 2 - 5 && this.mdowncoords[1] < this.pausebtny + this.pausebtnheight / 2 + 5;
	}

	checkAndHandleNavigation() {
		if (this.mdowncoords[0] > this.navx && this.mdowncoords[0] < this.navx + this.navbtnwidth && this.mdowncoords[1] > this.navy && this.mdowncoords[1] < this.navy + this.navbtnheight) {
			this.translatex += 20;
		} else if (this.mdowncoords[0] > this.navx - this.navbtnwidth && this.mdowncoords[0] < this.navx &&
			this.mdowncoords[1] > this.navy + this.navbtnheight && this.mdowncoords[1] < this.navy + this.navbtnheight * 2) {
			this.translatey += 20;
		} else if (this.mdowncoords[0] > this.navx - this.navbtnwidth * 2 && this.mdowncoords[0] < this.navx - this.navbtnwidth * 2 + this.navbtnwidth &&
			this.mdowncoords[1] > this.navy && this.mdowncoords[1] < this.navy + this.navbtnheight) {
			this.translatex -= 20;
		} else if (this.mdowncoords[0] > this.navx - this.pausebtnwidth && this.mdowncoords[0] < this.navx &&
			this.mdowncoords[1] > this.navy - this.navbtnheight && this.mdowncoords[1] < this.navy) {
			this.translatey -= 20; // up
		}
	}

	handlePause() {
		this.running ? this.stop() : this.continu();
	}

	handleZoom(event) {
		var zoom = handleSlider(event, this.pausebtnx + this.pausebtnwidth + 10, 140, 3, 500, 1, 0, this.canvas.getBoundingClientRect().left);
		this.scalex = zoom;
		this.scaley = zoom;
	}

	start() { this.current_n = 1; this.running = true; this.drawLoop(); }

	continu() { this.running = true; this.drawLoop(); }

	stop() { this.running = false; }

	drawFirstFrame() {
		this.start();
		this.stop();
	}

	coordToScreen(x, y) {
		return [this.width / 2 + this.scalex * x - this.translatex, this.height / 2 + this.scaley * y * (-1) - this.translatey];
	}

	screenToCoord(x, y) {
		var leftestX = ((this.width / 5) / 2) * -1;
		var uppestY = ((this.height / 5) / 2);

		var xCoord = (x) / 5 + leftestX;
		var yCoord = (y) / 5 + uppestY;

		return [xCoord, yCoord];
	}

	updateOrientation() {
		this.canvas.width = this.canvas.width = this.jquerycanvas.parent().width();
		this.width = this.canvas.width = this.jquerycanvas.parent().width();
		this.navx = this.width - this.navbtnwidth - 5;
		this.navy = this.height - 105;
		this.drawFirstFrame();
	}

	drawNavigationButton(x, y, w, h, xoff, yoff, x1, x2, x3, y1, y2, y3) {
		this.ctx.fillStyle = 'rgba(0, 0, 0, 0.3)';
		this.ctx.fillRect(x, y, w, h);
		this.ctx.fill();
		this.ctx.strokeRect(x, y, w, h);
		this.ctx.stroke();
		this.ctx.fillStyle = 'rgba(0, 0, 0, 1)';
		this.ctx.beginPath();
		this.ctx.moveTo(this.navx + xoff + x1 * this.navbtnwidth, this.navy + yoff + y1 * this.navbtnheight);
		this.ctx.lineTo(this.navx + xoff + x2 * this.navbtnwidth, this.navy + yoff + y2 * this.navbtnheight);
		this.ctx.lineTo(this.navx + xoff + x3 * this.navbtnwidth, this.navy + yoff + y3 * this.navbtnheight);
		this.ctx.closePath();
		this.ctx.fill();
	}

	drawLoop() {

		this.translatex *= -1;
		this.translatey *= -1;

		// white background
		this.ctx.fillStyle = 'rgba(255, 255, 255, 1)';
		this.ctx.fillRect(0, 0, this.width, this.height);

		this.ctx.translate(this.translatex, this.translatey);

		// grid
		var nwlines = this.width / this.scalex;
		var nhlines = this.height / this.scaley;

		var maxTrans = Math.max(this.translatex, this.translatey);
		var minTrans = Math.min(this.translatex, this.translatey);

		var txcoord = this.width / 2;
		var tycoord = this.height / 2;

		if (Math.min(this.width / 2 + this.translatex, 25) == this.width / 2 + this.translatex) {
			txcoord = 0 - this.translatex + 25;
		} else if (Math.max(this.width / 2 + this.translatex, this.width - 25) == this.width / 2 + this.translatex) {
			txcoord = this.width - this.translatex - 25;
		}

		if (Math.min(this.height / 2 + this.translatey, 25) == this.height / 2 + this.translatey) {
			tycoord = 0 - this.translatey + 25;
		} else if (Math.max(this.height / 2 + this.translatey, this.height - 25) == this.height / 2 + this.translatey) {
			tycoord = this.height - this.translatey - 25;
		}

		this.ctx.lineWidth = 1;
		this.ctx.strokeStyle = 'rgba(0, 0, 0, 0.3)';
		this.ctx.fillStyle = 'rgba(0, 0, 0, 1)';
		this.ctx.font = '10px times';

		// xgrid
		for (var i = (-1) * Math.floor(nwlines / 2) - Math.floor(maxTrans / this.scalex) - 1; i < nwlines / 2 - Math.floor(minTrans / this.scalex); i++) {
			this.ctx.beginPath();
			this.ctx.moveTo(this.width / 2 + this.scalex * i, 0 - maxTrans);
			this.ctx.lineTo(this.width / 2 + this.scalex * i, this.height - minTrans);
			if (i % this.xdense == 0) this.ctx.fillText(i, this.width / 2 + this.scalex * i + 1, tycoord + 9);
			this.ctx.stroke();
		}

		// ygrid
		for (var i = (-1) * Math.floor(nhlines / 2) - Math.floor(maxTrans / this.scaley) - 1; i < nhlines / 2 - Math.floor(minTrans / this.scaley); i++) {
			this.ctx.beginPath();
			this.ctx.moveTo(0 - maxTrans, this.height / 2 + this.scaley * i);
			this.ctx.lineTo(this.width - minTrans, this.height / 2 + this.scaley * i);
			if (i % this.ydense == 0) this.ctx.fillText(i * (-1), txcoord + 1, this.height / 2 + this.scaley * i);
			this.ctx.stroke();
		}

		// x coordinate
		this.ctx.beginPath();
		this.ctx.lineWidth = 2;
		this.ctx.strokeStyle = 'rgba(0, 0, 0, 1)';
		this.ctx.moveTo(5 - this.translatex, tycoord);
		this.ctx.lineTo(this.width - 10 - this.translatex, tycoord);
		this.ctx.stroke();

		// x coordinate arrow
		this.ctx.beginPath();
		this.ctx.moveTo(this.width - 10 - this.translatex, tycoord + 5);
		this.ctx.lineTo(this.width - 10 - this.translatex, tycoord - 5);
		this.ctx.lineTo(this.width - 5 - this.translatex, tycoord);
		this.ctx.fill();

		// x key
		this.ctx.font = '15px times';
		this.ctx.fillText(this.xkey, this.width - this.xkey.length * 7 - this.translatex, tycoord + 15);

		// y coordinate
		this.ctx.beginPath();
		this.ctx.moveTo(txcoord, 10 - this.translatey);
		this.ctx.lineTo(txcoord, this.height - 5 - this.translatey);
		this.ctx.stroke();

		// y coordinate arrow
		this.ctx.beginPath();
		this.ctx.moveTo(txcoord + 5, 10 - this.translatey);
		this.ctx.lineTo(txcoord - 5, 10 - this.translatey);
		this.ctx.lineTo(txcoord, 5 - this.translatey);
		this.ctx.fill();

		// y key
		this.ctx.beginPath();
		this.ctx.fillText(this.ykey, txcoord + 10, 15 - this.translatey);

		this.ctx.translate(-this.translatex, -this.translatey);

		this.translatex *= -1;
		this.translatey *= -1;

		// draw cartasian function
		this.drawCartFunction(this);

		// draw gui
		this.drawGui(this);

		// draw zoom slider
		drawSlider(this.scalex, this.pausebtnx + this.pausebtnwidth + 10, 140, this.pausebtny + this.pausebtnheight / 2, 3, 500, 1, "Zoom", this.ctx);

		// draw pause/play
		this.ctx.fillStyle = 'rgba(0, 0, 0, 0.3)';
		this.ctx.fillRect(this.pausebtnx, this.pausebtny, this.pausebtnwidth, this.pausebtnheight);
		this.ctx.strokeRect(this.pausebtnx, this.pausebtny, this.pausebtnwidth, this.pausebtnheight);
		this.ctx.fill();
		if (this.running) {
			this.ctx.lineWidth = 5;
			this.ctx.beginPath();
			this.ctx.moveTo(this.pausebtnx + 0.4 * this.pausebtnwidth, this.pausebtny + 0.1 * this.pausebtnheight);
			this.ctx.lineTo(this.pausebtnx + 0.4 * this.pausebtnwidth, this.pausebtny + 0.9 * this.pausebtnheight);
			this.ctx.stroke();
			this.ctx.beginPath();
			this.ctx.moveTo(this.pausebtnx + 0.6 * this.pausebtnwidth, this.pausebtny + 0.1 * this.pausebtnheight);
			this.ctx.lineTo(this.pausebtnx + 0.6 * this.pausebtnwidth, this.pausebtny + 0.9 * this.pausebtnheight);
			this.ctx.stroke();
		} else {
			this.ctx.fillStyle = 'rgba(0, 0, 0, 1)';
			this.ctx.beginPath();
			this.ctx.moveTo(this.pausebtnx + 0.3 * this.pausebtnwidth, this.pausebtny + 0.1 * this.pausebtnheight);
			this.ctx.lineTo(this.pausebtnx + 0.3 * this.pausebtnwidth, this.pausebtny + 0.9 * this.pausebtnheight);
			this.ctx.lineTo(this.pausebtnx + 0.8 * this.pausebtnwidth, this.pausebtny + 0.5 * this.pausebtnheight);
			this.ctx.closePath();
			this.ctx.fill();
		}

		// draw navigation
		this.ctx.lineWidth = 2;
		// right
		this.drawNavigationButton(this.navx, this.navy, this.navbtnwidth, this.navbtnheight, 0, 0, 0.3, 0.3, 0.8, 0.1, 0.9, 0.5);
		// down
		this.drawNavigationButton(this.navx - this.navbtnwidth, this.navy + this.navbtnheight, this.navbtnwidth,
			this.navbtnheight, -this.pausebtnwidth, this.pausebtnheight, 0.1, 0.9, 0.5, 0.3, 0.3, 0.8);
		// left
		this.drawNavigationButton(this.navx - this.navbtnwidth * 2, this.navy, this.navbtnwidth,
			this.navbtnheight, -this.pausebtnwidth * 2, 0, 0.7, 0.7, 0.2, 0.1, 0.9, 0.5);
		// up
		this.drawNavigationButton(this.navx - this.navbtnwidth, this.navy - this.navbtnheight, this.navbtnwidth,
			this.navbtnheight, -this.pausebtnwidth, -this.pausebtnheight, 0.1, 0.9, 0.5, 0.7, 0.7, 0.2);

		// check for nav hold down
		if (this.mousedown) this.checkAndHandleNavigation();

		// draw current_n
		this.ctx.fillStyle = "rgba(0, 0, 0, 1)";
		this.ctx.font = '20px times';
		this.ctx.fillText("n = " + this.current_n, 50, 50);

		if (this.current_n < this.n_max) this.current_n++;
		else this.current_n = 0;

		if (this.running) { setTimeout(function () { this.drawLoop(this.current_n, this.data); }.bind(this), 32); }
	}
}


// https://prime-numbers.info/list/primes
const primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997];

// --- ZETA HARMONIC ---
function drawZetaHarmonic(coordsys) {
	var ctx = coordsys.ctx;

	ctx.save();
	var travprim = 0;
	for (var x = 0; x < 1000; x++) {
		ctx.fillStyle = "rgba(0, 0, 0, 1)";
		// is x power of prime? then incr travprim

		for (var i = 0; i < primes.length; i++) {
			for (var k = 1; k < 5; k++) {
				if (x == Math.pow(primes[i], k)) {
					travprim += Math.log(primes[i]);
				}
			}
		}
		var origCoords = coordsys.coordToScreen(x, travprim);
		ctx.fillRect(origCoords[0], origCoords[1], 2, 2);
		ctx.stroke();
	}
	ctx.restore();


	ctx.save();
	for (var x = 0; x < 200; x++) {
		ctx.fillStyle = "rgba(255, 0, 0, 1)";

		var screencoords = coordsys.coordToScreen(x, chebyshev_psi_zero(x, coordsys.current_n));
		ctx.fillRect(screencoords[0], screencoords[1], 3, 3);
		ctx.stroke();
	}
	ctx.restore();
}

// http://empslocal.ex.ac.uk/people/staff/mrwatkin/zeta/encoding2.htm
var zeta_zeros = [14.134725142, 21.022039639, 25.010857580, 30.424876126, 32.935061588, 37.586178159, 40.918719012, 43.327073281, 48.005150881, 49.773832478, 52.970321478, 56.446247697, 59.347044003, 60.831778525, 65.112544048, 67.079810529, 69.546401711, 72.067157674, 75.704690699, 77.144840069, 79.337375020, 82.910380854, 84.735492981, 87.425274613, 88.809111208, 92.491899271, 94.651344041, 95.870634228, 98.831194218, 101.317851006, 103.725538040, 105.446623052, 107.168611184, 111.029535543, 111.874659177, 114.320220915, 116.226680321, 118.790782866, 121.370125002, 122.946829294, 124.256818554, 127.516683880, 129.578704200, 131.087688531, 133.497737203, 134.756509753, 138.116042055, 139.736208952, 141.123707404, 143.111845808, 146.000982487, 147.422765343, 150.053520421, 150.925257612, 153.024693811, 156.112909294, 157.597591818, 158.849988171, 161.188964138, 163.030709687, 165.537069188, 167.184439978, 169.094515416, 169.911976479, 173.411536520, 174.754191523, 176.441434298, 178.377407776, 179.916484020, 182.207078484, 184.874467848, 185.598783678, 187.228922584, 189.416158656, 192.026656361, 193.079726604, 195.265396680, 196.876481841, 198.015309676, 201.264751944, 202.493594514, 204.189671803, 205.394697202, 207.906258888, 209.576509717, 211.690862595, 213.347919360, 214.547044783, 216.169538508, 219.067596349, 220.714918839, 221.430705555, 224.007000255, 224.983324670, 227.421444280, 229.337413306, 231.250188700, 231.987235253, 233.693404179, 236.524229666, 237.769820481, 239.555477573, 241.049157796, 242.823271934, 244.070898497, 247.136990075, 248.101990060, 249.573689645, 251.014947795, 253.069986748, 255.306256455, 256.380713694, 258.610439492, 259.874406990, 260.805084505, 263.573893905, 265.557851839, 266.614973782, 267.921915083, 269.970449024, 271.494055642, 273.459609188, 275.587492649, 276.452049503, 278.250743530, 279.229250928, 282.465114765, 283.211185733, 284.835963981, 286.667445363, 287.911920501, 289.579854929, 291.846291329, 293.558434139, 294.965369619, 295.573254879, 297.979277062, 299.840326054, 301.649325462, 302.696749590, 304.864371341, 305.728912602, 307.219496128, 310.109463147, 311.165141530, 312.427801181, 313.985285731, 315.475616089, 317.734805942, 318.853104256, 321.160134309, 322.144558672, 323.466969558, 324.862866052, 327.443901262, 329.033071680, 329.953239728, 331.474467583, 333.645378525];

function chebyshev_psi_zero(x, limit) {

	var secndterm = Math.log(2 * Math.PI);
	var fourthterm = 0.5 * Math.log(1 - (1 / Math.pow(x, 2)));
	var thirdterm = 0;

	var sum = new Complex(0, 0);
	for (var rho = 0; rho < limit; rho++) {
		var rho_val = new Complex(0.5, zeta_zeros[rho]);
		var rho_val_conjugate = new Complex(0.5, (-1) * zeta_zeros[rho]);
		sum = add_complex(sum, div_complex(pow_complex(x, rho_val), rho_val));
		sum = add_complex(sum, div_complex(pow_complex(x, rho_val_conjugate), rho_val_conjugate));
	}

	return x - secndterm - sum.real - fourthterm;
}
// --- ---

// --- NONTRIVIAL_ZETAZEROS ---
function drawNontrivialZetazeros(coordsys) {
	var ctx = coordsys.ctx;

	ctx.save();
	for (var y = coordsys.data.imgmin; y < coordsys.data.imgmax; y += 0.05) {
		var x = 0.5;
		ctx.fillStyle = "rgba(0, 0, 0, 1)";
		var origCoords = coordsys.coordToScreen(x, y);
		ctx.fillRect(origCoords[0], origCoords[1], 2, 2);
		ctx.stroke();
		var z = zeta_extended(new Complex(x, y), coordsys.current_n);
		var destCoords = coordsys.coordToScreen(z.real, z.imaginary);
		if (z.real > -0.1 && z.real < 0.1 && z.imaginary > -0.1 && z.imaginary < 0.1) {
			ctx.fillStyle = "rgba(0, 0, 255, 1)";
		} else {
			ctx.fillStyle = "rgba(255, 0, 0, 1)";
		}
		ctx.fillRect(destCoords[0], destCoords[1], 2, 2);
		ctx.stroke();
	}
	ctx.restore();
}

function handleNontrivialZetaSlider(event, coordsys) {
	var mousey = event.clientY - this.canvas.getBoundingClientRect().top;
	var canvasBoundingLeft = this.canvas.getBoundingClientRect().left;

	if (mousey <= this.sliderposyStart + this.slideryspace / 2) {
		coordsys.data.imgmin = handleSlider(event, this.sliderposxStart, this.sliderlength, -50, 50, 0.05, 2, canvasBoundingLeft);
	} else if (mousey <= this.sliderposyStart + this.slideryspace + this.slideryspace / 2) {
		coordsys.data.imgmax = handleSlider(event, this.sliderposxStart, this.sliderlength, -50, 50, 0.05, 2, canvasBoundingLeft);
	} else if (mousey <= this.sliderposyStart + this.slideryspace * 2 + this.slideryspace / 2) {
		coordsys.n_max = handleSlider(event, this.sliderposxStart, this.sliderlength, 1, 2000, 1, 0, canvasBoundingLeft);
	}
}

function drawNontrivialZetaGui(coordsys) {
	drawSlider(coordsys.data.imgmin, this.sliderposxStart, this.sliderlength, this.sliderposyStart, -50, 50, 0.05, "img_min", coordsys.ctx);
	drawSlider(coordsys.data.imgmax, this.sliderposxStart, this.sliderlength, this.sliderposyStart + this.slideryspace, -50, 50, 0.05, "img_max", coordsys.ctx);
	drawSlider(coordsys.n_max, this.sliderposxStart, this.sliderlength, this.sliderposyStart + (this.slideryspace * 2), 0, 2000, 1, "n_max", coordsys.ctx);
}
// --- ---


// --- ANALYTIC CONTINUATION ---
function drawAnalyticContinuation(coordsys) {
	var ctx = coordsys.ctx;

	ctx.save();
	for (var y = coordsys.data.imgmin; y < coordsys.data.imgmax; y += 0.05) {
		for (var x = coordsys.data.realmin; x < coordsys.data.realmax; x += 0.05) {
			ctx.fillStyle = "rgba(0, 0, 0, 1)";
			var origCoords = coordsys.coordToScreen(x, y);
			ctx.fillRect(origCoords[0], origCoords[1], 2, 2);
			ctx.stroke();
			var cinp = new Complex(x, y);
			var z = zeta_extended(cinp, coordsys.current_n);
			var destCoords = coordsys.coordToScreen(z.real, z.imaginary);
			ctx.fillStyle = "rgba(255, 0, 0, 1)";
			ctx.fillRect(destCoords[0], destCoords[1], 2, 2);

			ctx.stroke();
		}
	}
	ctx.restore();
}

function handleAnalyticConSlider(event, coordsys) {
	var mousey = event.clientY - this.canvas.getBoundingClientRect().top;
	var canvasBoundingLeft = this.canvas.getBoundingClientRect().left;

	if (mousey <= this.sliderposyStart + this.slideryspace / 2) {
		coordsys.data.realmin = handleSlider(event, this.sliderposxStart, this.sliderlength, -5, 5, 0.05, 2, canvasBoundingLeft);
	} else if (mousey <= this.sliderposyStart + this.slideryspace + this.slideryspace / 2) {
		coordsys.data.realmax = handleSlider(event, this.sliderposxStart, this.sliderlength, -5, 5, 0.05, 2, canvasBoundingLeft);
	} else if (mousey <= this.sliderposyStart + this.slideryspace * 2 + this.slideryspace / 2) {
		coordsys.data.imgmin = handleSlider(event, this.sliderposxStart, this.sliderlength, -5, 5, 0.05, 2, canvasBoundingLeft);
	} else if (mousey <= this.sliderposyStart + this.slideryspace * 3 + this.slideryspace / 2) {
		coordsys.data.imgmax = handleSlider(event, this.sliderposxStart, this.sliderlength, -5, 5, 0.05, 2, canvasBoundingLeft);
	} else if (mousey <= this.sliderposyStart + this.slideryspace * 4 + this.slideryspace / 2) {
		coordsys.n_max = handleSlider(event, this.sliderposxStart, this.sliderlength, 1, 2000, 1, 0, canvasBoundingLeft);
	}
}

function drawAnalyticConGui(coordsys) {
	drawSlider(coordsys.data.realmin, this.sliderposxStart, this.sliderlength, this.sliderposyStart, -5, 5, 0.05, "real_min", coordsys.ctx);
	drawSlider(coordsys.data.realmax, this.sliderposxStart, this.sliderlength, this.sliderposyStart + (this.slideryspace * 1), -5, 5, 0.05, "real_max", coordsys.ctx);
	drawSlider(coordsys.data.imgmin, this.sliderposxStart, this.sliderlength, this.sliderposyStart + (this.slideryspace * 2), -5, 5, 0.05, "img_min", coordsys.ctx);
	drawSlider(coordsys.data.imgmax, this.sliderposxStart, this.sliderlength, this.sliderposyStart + (this.slideryspace * 3), -5, 5, 0.05, "img_max", coordsys.ctx);
	drawSlider(coordsys.n_max, this.sliderposxStart, this.sliderlength, this.sliderposyStart + (this.slideryspace * 4), 0, 2000, 1, "n_max", coordsys.ctx);
}

// https://en.wikipedia.org/wiki/Riemann_zeta_function#Representations
// converges for Re(z) > 0
function zeta_extended(s, limit) {
	var leftprod = div_complex(1, sub_complex(s, 1));
	var sum = new Complex(0, 0);
	for (var n = 1; n < limit; n++) {
		var powc = pow_complex(n + 1, s);
		var left_sub = div_complex(n, powc);
		var right_sub = div_complex(sub_complex(n, s), pow_complex(n, s));
		sum = add_complex(sum, sub_complex(left_sub, right_sub));
	}
	return mult_complex(leftprod, sum);
}
// --- ---

// --- TRANSFORMATION ---
function drawTransformation(coordsys) {
	var ctx = coordsys.ctx;

	ctx.save();
	for (var y = coordsys.data.imgmin; y < coordsys.data.imgmax; y += 0.05) {
		for (var x = coordsys.data.realmin; x < coordsys.data.realmax; x += 0.05) {

			ctx.fillStyle = "rgba(0, 0, 0, 1)";
			var origCoords = coordsys.coordToScreen(x, y);
			ctx.fillRect(origCoords[0], origCoords[1], 2, 2);
			ctx.stroke();

			var z = zeta(new Complex(x, y), coordsys.current_n);

			var destCoords = coordsys.coordToScreen(z.real, z.imaginary);
			ctx.fillStyle = "rgba(255, 0, 0, 1)";
			ctx.fillRect(destCoords[0], destCoords[1], 2, 2);
			ctx.stroke();
		}
	}
	ctx.restore();
}

function handleTransformationSlider(event, coordsys) {
	var mousey = event.clientY - this.canvas.getBoundingClientRect().top;
	var canvasBoundingLeft = this.canvas.getBoundingClientRect().left;

	if (mousey <= this.sliderposyStart + this.slideryspace / 2) {
		coordsys.data.realmin = handleSlider(event, this.sliderposxStart, this.sliderlength, -5, 5, 0.05, 2, canvasBoundingLeft);
	} else if (mousey <= this.sliderposyStart + this.slideryspace + this.slideryspace / 2) {
		coordsys.data.realmax = handleSlider(event, this.sliderposxStart, this.sliderlength, -5, 5, 0.05, 2, canvasBoundingLeft);
	} else if (mousey <= this.sliderposyStart + this.slideryspace * 2 + this.slideryspace / 2) {
		coordsys.data.imgmin = handleSlider(event, this.sliderposxStart, this.sliderlength, -5, 5, 0.05, 2, canvasBoundingLeft);
	} else if (mousey <= this.sliderposyStart + this.slideryspace * 3 + this.slideryspace / 2) {
		coordsys.data.imgmax = handleSlider(event, this.sliderposxStart, this.sliderlength, -5, 5, 0.05, 2, canvasBoundingLeft);
	} else if (mousey <= this.sliderposyStart + this.slideryspace * 4 + this.slideryspace / 2) {
		coordsys.n_max = handleSlider(event, this.sliderposxStart, this.sliderlength, 1, 2000, 1, 0, canvasBoundingLeft);
	}
}

function drawTransformationGui(coordsys) {
	drawSlider(coordsys.data.realmin, this.sliderposxStart, this.sliderlength, this.sliderposyStart, -5, 5, 0.05, "real_min", coordsys.ctx);
	drawSlider(coordsys.data.realmax, this.sliderposxStart, this.sliderlength, this.sliderposyStart + (this.slideryspace * 1), -5, 5, 0.05, "real_max", coordsys.ctx);
	drawSlider(coordsys.data.imgmin, this.sliderposxStart, this.sliderlength, this.sliderposyStart + (this.slideryspace * 2), -5, 5, 0.05, "img_min", coordsys.ctx);
	drawSlider(coordsys.data.imgmax, this.sliderposxStart, this.sliderlength, this.sliderposyStart + (this.slideryspace * 3), -5, 5, 0.05, "img_max", coordsys.ctx);
	drawSlider(coordsys.n_max, this.sliderposxStart, this.sliderlength, this.sliderposyStart + (this.slideryspace * 4), 0, 2000, 1, "n_max", coordsys.ctx);
}
// --- ---


// --- ZETA CONVERGENCE ---
function drawZeta(coordsys) {
	var ctx = coordsys.ctx;
	ctx.save();
	ctx.fillStyle = 'rgba(' + coordsys.data.real * 255 / 5 + ', 0, ' + coordsys.data.img * 255 / 5 + ', ' + 255 + ')';
	for (var limit = 1; limit < coordsys.current_n; limit++) {
		var z = zeta(new Complex(coordsys.data.real, coordsys.data.img), limit);
		var z_screencoords = coordsys.coordToScreen(z.real, z.imaginary);
		ctx.fillRect(z_screencoords[0], z_screencoords[1], 2, 2);
		ctx.stroke();
	}
	ctx.restore();
}

function handleZetaSlider(event, coordsys) {
	var mousey = event.clientY - this.canvas.getBoundingClientRect().top;
	var canvasBoundingLeft = this.canvas.getBoundingClientRect().left;

	if (mousey <= this.sliderposyStart + this.slideryspace / 2) {
		coordsys.data.real = handleSlider(event, this.sliderposxStart, this.sliderlength, -5, 5, 0.05, 2, canvasBoundingLeft);
	} else if (mousey <= this.sliderposyStart + this.slideryspace + this.slideryspace / 2) {
		coordsys.data.img = handleSlider(event, this.sliderposxStart, this.sliderlength, -5, 5, 0.05, 2, canvasBoundingLeft);
	} else if (mousey <= this.sliderposyStart + this.slideryspace * 2 + this.slideryspace / 2) {
		coordsys.n_max = handleSlider(event, this.sliderposxStart, this.sliderlength, 1, 2000, 1, 0, canvasBoundingLeft);
	}
}

function drawZetaGui(coordsys) {
	drawSlider(coordsys.data.real, this.sliderposxStart, this.sliderlength, this.sliderposyStart, -5, 5, 0.05, "a", coordsys.ctx);
	drawSlider(coordsys.data.img, this.sliderposxStart, this.sliderlength, this.sliderposyStart + (this.slideryspace * 1), -5, 5, 0.05, "b", coordsys.ctx);
	drawSlider(coordsys.n_max, this.sliderposxStart, this.sliderlength, this.sliderposyStart + (this.slideryspace * 2), 0, 2000, 1, "n_max", coordsys.ctx);
}

function zeta(z, limit) {
	var zres = new Complex(0, 0);
	for (var x = 1; x <= limit; x++) {
		var ii = z.imaginary * Math.log(1 / x);
		var pp = Math.pow(1 / x, z.real);
		zres.real += pp * Math.cos(ii);
		zres.imaginary += pp * Math.sin(ii);
	}
	return zres;
}

// --- CREATE MEDIA RECORDER ---
function createMediaRecorder(canvasname, videoplayername, buttonname) {
	var mediaRecorder = new MediaRecorder($("#" + canvasname)[0].captureStream(60));
	var chunks = [];
	mediaRecorder.ondataavailable = function (e) { chunks.push(e.data); };
	mediaRecorder.onstop = function (e) {
		var blob = new Blob(chunks, { 'type': 'video/mp4' });
		chunks = [];
		var videoURL = URL.createObjectURL(blob);
		$("#" + videoplayername)[0].src = videoURL;
		setTimeout(function () {
			$("#" + canvasname).css("display", "none");
			$("#" + videoplayername).css("display", "block");
			$("#" + buttonname).css("display", "block");
		}
			, 2000);
	};
	return mediaRecorder;
}
// --- ---
/*
var transMediaRec = createMediaRecorder("transformation", "transVideoPlayer", "transButton");
transMediaRec.start();
setTimeout(function (){transMediaRec.stop();}, 5000);
addSimReturn("transformation", "transVideoPlayer", "transButton", transformation_coordsys, drawTransformation);
*/

// --- ADD RETURN TO SIMULATION BUTTON LISTENER ---
function addSimReturn(canvasname, videoplayername, buttonname, coordsys, loopfunc) {
	$("#" + buttonname).click(function () {
		coordsys.drawLoop(0, loopfunc, -1);
		$("#" + videoplayername).css("display", "none");
		$("#" + buttonname).css("display", "none");
		$("#" + canvasname).css("display", "block");
	});
}
// --- ---

// --- DRAW SLIDER ---
function drawSlider(value, sliderposxStart, sliderlength, sliderposy, bot, top, step, key, ctx) {
	var sliderposxEnd = sliderposxStart + sliderlength;

	ctx.save();
	ctx.beginPath();
	ctx.moveTo(sliderposxStart, sliderposy);
	ctx.lineTo(sliderposxEnd, sliderposy);
	ctx.stroke();
	ctx.beginPath();

	var totalrange = Math.abs(Math.max(bot, top) - Math.min(bot, top));
	var totalunits = totalrange / step;
	var segment = value / step;

	var sliderlength = sliderposxEnd - sliderposxStart;
	var sliderstep = sliderlength / totalunits;

	var xx = segment * sliderstep;

	var offset = (-bot / step) * sliderstep;

	ctx.arc(sliderposxStart + xx + offset, sliderposy, 4, 0, 2 * Math.PI);
	ctx.fill();
	ctx.font = '20px times';
	ctx.fillText(key + " = " + value, sliderposxStart + xx + offset, sliderposy - 10);
	ctx.restore();
}
// --- ---

// --- HANDLE SLIDER ---
function handleSlider(event, x, length, bot, top, step, decimalprecision, boundingClientRectLeft) {

	var sliderposxEnd = x + length;

	var totalrange = Math.abs(Math.max(bot, top) - Math.min(bot, top));
	var totalunits = totalrange / step;
	var sliderlength = sliderposxEnd - x;
	var sliderstep = length / totalunits;
	var offset = (-bot / step) * sliderstep;

	var xx = event.clientX - boundingClientRectLeft - x - offset;
	var segment = xx / sliderstep;

	return clamp(parseFloat((segment * step).toFixed(decimalprecision)), bot, top);
}
// --- ---


window.addEventListener("orientationchange", function () {
	setTimeout(function () { zeta_convergence_coordsys.updateOrientation(); }, 500);
	setTimeout(function () { transformation_coordsys.updateOrientation(); }, 500);
	setTimeout(function () { analytic_cont.updateOrientation(); }, 500);
	setTimeout(function () { nontrivial_zetazeros.updateOrientation(); }, 500);
	setTimeout(function () { zeta_harmonic.updateOrientation(); }, 500);
});


var zeta_convergence_coordsys = new CartesianCoordinateSystem("zeta_convergence", 300, 300, 1, 1, 0, 0, drawZeta, drawZetaGui, handleZetaSlider, 50, 100, 140, 3);
zeta_convergence_coordsys.data.real = 1.3;
zeta_convergence_coordsys.data.img = 2;
zeta_convergence_coordsys.drawFirstFrame();


var transformation_coordsys = new CartesianCoordinateSystem("transformation", 100, 100, 2, 2, 100, 0, drawTransformation, drawTransformationGui, handleTransformationSlider, 50, 100, 140, 5);
transformation_coordsys.data.realmin = 1.05;
transformation_coordsys.data.realmax = 3;
transformation_coordsys.data.imgmin = -3;
transformation_coordsys.data.imgmax = 3;
transformation_coordsys.drawFirstFrame();


var analytic_cont = new CartesianCoordinateSystem("analytic_continuation", 100, 100, 2, 2, 100, 0, drawAnalyticContinuation, drawAnalyticConGui, handleAnalyticConSlider, 50, 100, 140, 5);
analytic_cont.data.realmin = 0.05;
analytic_cont.data.realmax = 1;
analytic_cont.data.imgmin = -3;
analytic_cont.data.imgmax = 3;
analytic_cont.drawFirstFrame();


var nontrivial_zetazeros = new CartesianCoordinateSystem("nontrivzetazeros", 100, 100, 2, 2, 100, 0, drawNontrivialZetazeros, drawNontrivialZetaGui, handleNontrivialZetaSlider, 50, 100, 140, 3);
nontrivial_zetazeros.data.imgmin = -50;
nontrivial_zetazeros.data.imgmax = 50;
nontrivial_zetazeros.drawFirstFrame();


var zeta_harmonic = new CartesianCoordinateSystem("zetaharmonic", 5, 5, 5, 5, 210, -130, drawZetaHarmonic, function () { }, function () { }, -1, -1, -1, -1);
zeta_harmonic.xkey = "x";
zeta_harmonic.ykey = "y";
zeta_harmonic.drawFirstFrame();

