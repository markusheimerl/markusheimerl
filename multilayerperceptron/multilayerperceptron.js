/*
* Copyright 2021 Markus Heimerl, OTH Regensburg
* Licensed under CC BY-NC 4.0
*
* ANY USE OF THIS SOFTWARE MUST COMPLY WITH THE
* CREATIVE COMMONS ATTRIBUTION-NONCOMMERCIAL 4.0 INTERNATIONAL LICENSE
* AGREEMENTS
*/
var malloc = Module.cwrap('malloc', 'number', ['number']);
var predict = Module.cwrap('predict', 'number', ['string', 'number']);
var getImages = Module.cwrap('getImages', 'number', ['string', 'number']);
var evaluateNetwork = Module.cwrap('evaluateNetwork', 'number', ['string', 'number']);
var printImage = Module.cwrap('printImage', 'number', ['number']);
var getAllNeuronsAndWeights = Module.cwrap('getAllNeuronsAndWeights', 'number', ['string']);
var printSingleNeuron = Module.cwrap('printSingleNeuron', 'number', ['number']);
var getWeightValueAsString = Module.cwrap('getWeightValueAsString', 'string', ['number', 'number']);
var getBiasValueAsString = Module.cwrap('getBiasValueAsString', 'string', ['number']);
var getNeuronValueAsString = Module.cwrap('getNeuronValueAsString', 'string', ['number']);
var predictHardcoded = Module.cwrap('predictHardcoded', 'number', ['number', 'number']);
var stepFeedNetwork = Module.cwrap('stepFeedNetwork', 'number', ['number', 'number', 'number']);
var stepFeedForward = Module.cwrap('stepFeedForward', 'number', ['number', 'number', 'number', 'number']);
var fastForwardFeedForwardOfThisLayer = Module.cwrap('fastForwardFeedForwardOfThisLayer', 'number', ['number', 'number', 'number', 'number']);
var stepApplySigmoid = Module.cwrap('stepApplySigmoid', 'number', ['number', 'number']);

/*
if( /Android|webOS|iPhone|iPad|Mac|Macintosh|iPod|BlackBerry|IEMobile|Opera Mini/i.test(navigator.userAgent) ) {
	$("#multilayerperceptroncanvas")[0].style.display = "none";
	$("#multilayerperceptronfeedforwardtopinfo")[0].style.display = "none";
	$("#multilayerperceptronfeedforward")[0].style.display = "none";
}
*/

var canvas = $("#multilayerperceptroncanvas")[0];
var ctx = canvas.getContext("2d");
var neuronsobject = {};
var globalinputdatapointer = -1;
var infoimage = -1;
var ffphasestring = "None";
var running = false;

function readIntFromHeap(u8addr) {
	var addr = 0;
	addr = addr | HEAPU8[u8addr + 3];
	addr = addr << 8;
	addr = addr | HEAPU8[u8addr + 2];
	addr = addr << 8;
	addr = addr | HEAPU8[u8addr + 1];
	addr = addr << 8;
	addr = addr | HEAPU8[u8addr + 0];
	return addr;
}

function printNeuron(arr, neuronindex) {
	printSingleNeuron(readIntFromHeap(arr + 4 * neuronindex));
}

function getWeight(arr, neuronindex, weightindex, biasandweightsvaluebuffer) {
	return getWeightValueAsString(readIntFromHeap(arr + 4 * neuronindex), weightindex, biasandweightsvaluebuffer);
}

function getNeuron(arr, neuronindex, biasandweightsvaluebuffer) {
	return getNeuronValueAsString(readIntFromHeap(arr + 4 * neuronindex), biasandweightsvaluebuffer);
}

function getBias(arr, neuronindex, biasandweightsvaluebuffer) {
	return getBiasValueAsString(readIntFromHeap(arr + 4 * neuronindex), biasandweightsvaluebuffer);
}

function resizeTo(canvas, xscale, yscale, x_min, y_min) {
	var context = canvas.getContext("2d");
	var imageData = context.getImageData(0, 0, canvas.width, canvas.height);
	var newCanvas = $("<canvas>").attr("width", imageData.width/* *50 */).attr("height", imageData.height/* *50 had to be removed bc of firefox bug*/)[0];
	/*var newCanvas = document.createElement('canvas');
	newCanvas.width  = imageData.width;
	newCanvas.height = imageData.height;
	var newContext = newCanvas.getContext("2d");
	newContext.putImageData(imageData, 0, 0);*/

	newCanvas.getContext("2d").putImageData(imageData, 0, 0);

	context.clearRect(0, 0, canvas.width, canvas.height);
	context.scale(xscale, yscale);
	context.translate(-x_min, -y_min);
	context.drawImage(newCanvas, 0, 0);
}

function convertArrayToImageData(ctx, arr, x, y) {
	var imagedata = ctx.createImageData(x, y);
	var index = 3;
	for (var i = 0; i < 280 * 280; i++) {
		imagedata.data[index] = arr[i];
		index += 4;
	}
	return imagedata;
}

function convertImageDataToArray(imagedata) {
	var index = 3;
	var arr = [];
	for (var i = 0; i < imagedata.data.length / 4; i++) {
		arr[i] = imagedata.data[index];
		index += 4;
	}
	return arr;
}

function getMousePosition(e) {
	var mouseX = e.offsetX * canvas.width / canvas.clientWidth | 0;
	var mouseY = e.offsetY * canvas.height / canvas.clientHeight | 0;
	return { x: mouseX, y: mouseY };
}

var mousedown = false;
var touchevent = false;
$("#multilayerperceptroncanvas").mousedown(function (e) { mousedown = true; });
$("#multilayerperceptroncanvas").mouseup(function (e) { mousedown = false; });
$("#multilayerperceptroncanvas").mousemove(function (e) {
	if (mousedown) {
		if (touchevent) {
			ctx.lineTo(getMousePosition(e).x, getMousePosition(e).y);
			ctx.lineWidth = 10 * 2;
			ctx.stroke();
			ctx.beginPath();
			ctx.arc(getMousePosition(e).x, getMousePosition(e).y, 10, 0, Math.PI * 2);
			ctx.fill();
			ctx.beginPath();
			ctx.moveTo(getMousePosition(e).x, getMousePosition(e).y);
			touchevent = false;
		} else {
			ctx.beginPath();
			ctx.arc(e.offsetX, e.offsetY, 10, 0, 2 * Math.PI, false);
			ctx.fill();
		}
	}
});

// --- TOUCHEVENTS FOR DRAWING CANVAS ---
$("#multilayerperceptroncanvas")[0].addEventListener("touchstart", function (e) {
	var touch = e.touches[0];
	touchevent = true;
	var mouseEvent = new MouseEvent("mousedown", {
		clientX: touch.clientX,
		clientY: touch.clientY
	});
	$("#multilayerperceptroncanvas")[0].dispatchEvent(mouseEvent);
}, false);

$("#multilayerperceptroncanvas")[0].addEventListener("touchend", function (e) {
	var touch = e.touches[0];
	var mouseEvent = new MouseEvent("mouseup", {
		clientX: touch.clientX,
		clientY: touch.clientY
	});
	$("#multilayerperceptroncanvas")[0].dispatchEvent(mouseEvent);
}, false);
// --- ---



$("#recognizestart").click(function (e) {
	$(this).prop("disabled", true);
	var imagedata = ctx.getImageData(0, 0, canvas.width, canvas.height);
	var minified = new Uint8ClampedArray(imagedata.data.length / 4);
	var imgdataindex = 3;
	for (var i = 0; i < imagedata.data.length / 4; i++) {
		minified[i] = imagedata.data[imgdataindex];
		imgdataindex += 4;
	}

	// --- scale boudning box to 200x200 ---
	var x_min = 9999999;
	for (var i = 0; i < 280 * 280; i++) {
		if (minified[i] > 0 && ((i % 280) < x_min)) {
			x_min = i % 280;
		}
	}
	var x_max = -1;
	for (var i = 0; i < 280 * 280; i++) {
		if (minified[i] > 0 && ((i % 280) > x_max)) {
			x_max = i % 280;
		}
	}

	var y_min = 0;
	for (var i = 0; i < 280 * 280; i++) {
		if (minified[i] > 0) {
			y_min = parseInt(i / 280);
			break;
		}
	}
	var y_max = 0;
	for (var i = 280 * 280; i >= 0; i--) {
		if (minified[i] > 0) {
			y_max = parseInt(i / 280);
			break;
		}
	}

	var xscale = 200 / (x_max - x_min);
	var yscale = 200 / (y_max - y_min);
	var scale = 0;
	if (xscale > yscale) scale = yscale;
	else scale = xscale;
	resizeTo($("#multilayerperceptroncanvas")[0], scale, scale, x_min, y_min);
	var imgdat = $("#multilayerperceptroncanvas")[0].getContext("2d").getImageData(0, 0, canvas.width, canvas.height);
	minified = convertImageDataToArray(imgdat);
	// --- ---

	// --- center around center of mass ---
	// thank god: https://stackoverflow.com/a/37448909/9298528
	indexofmass = 0;
	x_sum = 0;
	y_sum = 0;
	datapoints = 0;
	for (var i = 0; i < 280 * 280; i++) {
		if (minified[i] > 0) {
			x_sum += i % 280;
			y_sum += parseInt(i / 280);
			datapoints++;
		}
	}
	x_sum = x_sum / datapoints;
	y_sum = y_sum / datapoints;
	// --- ---

	// --- move center of mass to 14,14 | in other words: center digit ---
	var xoff = 140 - x_sum;
	if (xoff < 0) xoff = 280 + xoff;
	for (var j = 0; j < xoff; j++) {
		for (var i = 280 * 280; i >= 0; i--) {
			minified[i] = minified[i - 1];
		}
	}

	var yoff = 140 - y_sum;
	var dir = 1;
	if (yoff < 0) dir = -1;

	if (dir == 1) {
		for (var j = 0; j < Math.abs(yoff); j++) {
			for (var i = 280 * 280 - 1; i >= 0; i--) {
				minified[i] = minified[i - 280];
			}
		}
	} else if (dir == -1) {
		for (var j = 0; j < Math.abs(yoff); j++) {
			for (var i = 0; i < 280 * 280; i++) {
				minified[i] = minified[i + 280];
			}
		}
	}
	// --- ---

	// --- minify ---
	var minifiedimagedata = ctx.createImageData(28, 28);
	for (var i = 3; i < 28 * 28 * 4; i += 4) {
		var total = 0;
		for (var j = 0; j < 100; j++) {
			var x = j % 10 + (((i - 3) / 4) * 10) % 280;
			var y = parseInt(j / 10) * 280 + (parseInt(((i - 3) / 4) / 28) * 10) * 280;
			total = total + minified[x + y];
		}

		total = total / 100;
		minifiedimagedata.data[i] = total;
	}
	infoimage = minifiedimagedata;
	// --- ---

	//$("#multilayerperceptroncanvas2")[0].getContext("2d").putImageData(minifiedimagedata, 0, 0);

	// --- upscale minfied version ---
	var upscaledimagedata = ctx.createImageData(280, 280);
	var upscaled = [];
	// blow up every pixel 10x10
	for (var i = 3; i < 28 * 28 * 4; i += 4) {
		// get pixel color for the next 100 colorings
		var col = minifiedimagedata.data[i];
		// color 10x10
		for (var j = 0; j < 10 * 10; j++) {
			var x = j % 10 + (((i - 3) / 4) * 10) % 280;
			var y = parseInt(j / 10) * 280 + (parseInt(((i - 3) / 4) / 28) * 10) * 280;
			upscaled[x + y] = col;
		}
	}

	var upindex = 3;
	for (var i = 0; i < 280 * 280; i++) {
		upscaledimagedata.data[upindex] = upscaled[i];
		upindex += 4;
	}
	$("#multilayerperceptroncanvas")[0].getContext("2d").putImageData(upscaledimagedata, 0, 0);
	// --- ---

	//--- FEED USER INPUT ---
	var datapointer = malloc(28 * 28);
	var minimgdataindex = 3;
	for (var i = 0; i < 28 * 28; i++) {
		HEAPU8[datapointer + i] = minifiedimagedata.data[minimgdataindex];
		minimgdataindex += 4;
	}
	//predict("assets/100_epoch_64_batchsize_relu_alpha0point5_784_200_10_95_13%.network", datapointer);
	predict("assets/100_epoch_64_batchsize_relu_alpha0point5_784_200_10_0_to_1_normalization_95_09%.network", datapointer);
	//predictHardcoded(neuronsobject.allneuronsarray, datapointer); // THIS IS THE THING THAT WORKS
	//printNeuron(neuronsobject.allneuronsarray, 784+199);
	globalinputdatapointer = datapointer;
	ffphasestring = "Feed Input Data";
	running = true;
	initNeurons();
});

/*
TODO:
- little arrow that says "network predicts you wrote a "3"
- Speed slider
- draw neuron thats pointed to by weightindexthatsred blue
- being able to repeat recognition without reloading the page
- make equation PNGs
- remove topoffset variable
- https://developer.mozilla.org/de/docs/Web/JavaScript/Reference/Statements/try...catch
DONE:
- USE NETWORK THATS TRAINED ON 0 TO 1 NORMALISATION NAMED 100_epoch_64_batchsize_relu_alpha0point5_784_200_10_0_to_1_normalization.network - DONE
- fast forward this layer button - DONE
- delete function.html - DONE
*/

// --- FEED FORWARD VISUALISATION ---
var loadingtext = "Loading...";
initNeurons();
function initNeurons() {
	try {
		neuronsobject.allneuronsarray = getAllNeuronsAndWeights("assets/100_epoch_64_batchsize_relu_alpha0point5_784_200_10_0_to_1_normalization_95_09%.network");
		feedforwardvisualisation();
	}
	catch (e) {
		loadingtext += ".";
		var canvasfftopinfo = $("#multilayerperceptronfeedforwardtopinfo")[0];
		var canvasmlpff = $("#multilayerperceptronfeedforward")[0];

		canvasfftopinfo.getContext("2d").clearRect(0, 0, canvasfftopinfo.width, canvasfftopinfo.height);
		canvasfftopinfo.getContext("2d").save();
		canvasfftopinfo.getContext("2d").font = "20px Arial";
		canvasfftopinfo.getContext("2d").fillText(loadingtext, 0, 30);
		canvasfftopinfo.getContext("2d").restore();

		canvasmlpff.getContext("2d").clearRect(0, 0, canvasmlpff.width, canvasmlpff.height);
		canvasmlpff.getContext("2d").save();
		canvasmlpff.getContext("2d").font = "20px Arial";
		canvasmlpff.getContext("2d").fillText(loadingtext, 0, 30);
		canvasmlpff.getContext("2d").restore();

		console.log(e);
		setTimeout(initNeurons, 500);
	}
}

function feedforwardvisualisation() {

	// so the loop should just keep on drawing. And draw a certain neuron red. We specify with this variable here which one that is.
	var neuronthatsredindex = -1;
	// also which of its weights is red (only neuron is red when this is -1)
	var weightthatsredindex = -1;

	var predictedneuron = -1;

	var neurondiameter = 40;
	var animationsteppingspeed = 100;

	var navbuttonwidth = 30;
	var navbuttonheight = 25;

	var fastforwardthislayerflag = false;

	// TOP CANVAS
	var canvastop = $("#multilayerperceptronfeedforwardtopinfo")[0];
	var ctxtop = canvastop.getContext("2d");
	var equationimage = new Image();
	var legendimage = new Image();
	legendimage.src = "legend.PNG";
	var fastforwardx = 920;
	var fastforwardy = 100;

	$("#multilayerperceptronfeedforwardtopinfo").click(function (e) {
		// fast forward button
		if (e.offsetX > fastforwardx && e.offsetX < fastforwardx + navbuttonheight * 2 + 8 && e.offsetY > fastforwardy - navbuttonwidth && e.offsetY < fastforwardy) {
			fastforwardthislayerflag = true;
		}
	});

	topinfoLoop();
	function topinfoLoop() {
		ctxtop.clearRect(0, 0, canvastop.width, canvastop.height);

		// phase
		ctxtop.lineWidth = 1;
		ctxtop.fillStyle = 'rgba(0, 0, 0, 1)';
		ctxtop.font = "20px Arial";
		ctxtop.save();

		ctxtop.fillText("Phase: ", 20, 40);
		ctxtop.fillText(ffphasestring, 20, 60);

		// info

		if (infoimage != -1 && (ffphasestring == "Feed Input Data" || ffphasestring == "Feed Forward")) {
			ctxtop.fillText("Info: ", 220, 40);
			ctxtop.putImageData(infoimage, 220, 50);
			ctxtop.strokeRect(220 - 1, 50 - 1, 31, 31);
			ctxtop.fillStyle = 'rgba(255, 0, 0, 1)';
			ctxtop.fillRect(220 + neuronthatsredindex % 28, 50 + parseInt(neuronthatsredindex / 28), 2, 2);
		}
		ctxtop.fillStyle = 'rgba(0, 0, 0, 1)';

		// equation
		if (ffphasestring == "Feed Input Data") {
			ctxtop.fillText("Equation:", 420, 40);
			equationimage.src = "equation_feedinputdataphase.png";
		}
		ctxtop.drawImage(equationimage, 520, 25);

		// speed
		ctxtop.fillText("Speed:", 20, 100);
		ctxtop.fillText(animationsteppingspeed, 20, 120);
		ctxtop.restore();

		// legend
		ctxtop.fillText("Legend: ", 420, 90);
		ctxtop.drawImage(legendimage, 420 + 90, 90 - 30);

		// fast forward this layer
		ctxtop.fillText("Fast Forward this layer: ", 700, 90);
		ctxtop.fillStyle = 'rgba(0, 0, 0, 0.3)';
		ctxtop.lineWidth = 2;
		ctxtop.beginPath();
		ctxtop.moveTo(fastforwardx, fastforwardy);
		ctxtop.lineTo(fastforwardx, fastforwardy - navbuttonwidth);
		ctxtop.lineTo(fastforwardx + navbuttonheight, fastforwardy - navbuttonwidth / 2);
		ctxtop.closePath();
		ctxtop.stroke();
		ctxtop.fill();
		ctxtop.beginPath();
		ctxtop.moveTo(fastforwardx + navbuttonheight + 3, fastforwardy);
		ctxtop.lineTo(fastforwardx + navbuttonheight + 3, fastforwardy - navbuttonwidth);
		ctxtop.lineTo(fastforwardx + navbuttonheight + navbuttonheight + 3, fastforwardy - navbuttonwidth / 2);
		ctxtop.closePath();
		ctxtop.stroke();
		ctxtop.fill();
		ctxtop.lineWidth = 3;
		ctxtop.beginPath();
		ctxtop.moveTo(fastforwardx + navbuttonheight + navbuttonheight + 7, fastforwardy);
		ctxtop.lineTo(fastforwardx + navbuttonheight + navbuttonheight + 7, fastforwardy - navbuttonwidth);
		ctxtop.stroke();

		if (running) requestAnimationFrame(topinfoLoop);
	}



	// BOTTOM CANVAS
	var topoffset = 0;
	var neuronoffset = 20;
	var bottomoffset = 50;
	var canvas = $("#multilayerperceptronfeedforward")[0];
	var ctx = canvas.getContext("2d");
	ctx.strokeStyle = "#000000";
	var xlayers = [100, 500, 900];
	var xnumberofneuronsperlayer = [784, 200, 10];
	var neuronbetweenspace = 5;
	var offset = [0, 90, 0];
	var biasandweightsvaluebuffer = malloc(10);

	$("#multilayerperceptronfeedforward").click(function (e) {
		for (var i = 0; i < 3; i++) {
			// top nav button
			if (e.offsetX > xlayers[i] - navbuttonwidth / 2 && e.offsetX < xlayers[i] + navbuttonwidth / 2 && e.offsetY > topoffset + neuronoffset - navbuttonheight / 2 && e.offsetY < topoffset + neuronoffset + navbuttonheight - navbuttonheight / 2) {
				offset[i] -= 9;
			}
			// bottom nav button
			if (e.offsetX > xlayers[i] - navbuttonwidth / 2 && e.offsetX < xlayers[i] + navbuttonwidth / 2 && e.offsetY < canvas.height + neuronoffset - bottomoffset - navbuttonheight / 2 + navbuttonheight && e.offsetY > canvas.height + neuronoffset - bottomoffset - navbuttonheight / 2) {
				offset[i] += 9;
			}
		}
	});

	// we kinda have to get all the weights
	// neuronsobject.allneuronsarray = getAllNeuronsAndWeights("assets/100_epoch_64_batchsize_relu_alpha0point5_784_200_10_95_13%.network");


	var sigmoidsteppingnetworkindex = 0;
	sigmoidStepping();
	function sigmoidStepping() {
		if (ffphasestring == "Apply sigmoid") {
			if (sigmoidsteppingnetworkindex == 10) {

				// find out what was predicted
				var prediction = -1;
				var max = 0.0;
				for (var i = 0; i < 10; i++) {
					var outputvalue = getNeuron(neuronsobject.allneuronsarray, 784 + 200 + i, biasandweightsvaluebuffer);
					if (outputvalue > max) { max = outputvalue; prediction = i; }
				}
				if (prediction == 9) offset[2] = 9;
				predictedneuron = prediction;

				neuronthatsredindex = -1; sigmoidsteppingnetworkindex = 0; ffphasestring = "Done"; setTimeout(sigmoidStepping, animationsteppingspeed);
				return;
			}
			neuronthatsredindex = 784 + 200 + sigmoidsteppingnetworkindex;
			stepApplySigmoid(neuronsobject.allneuronsarray, sigmoidsteppingnetworkindex);
			sigmoidsteppingnetworkindex++;
		}
		if (running) setTimeout(sigmoidStepping, animationsteppingspeed);
	}

	var weightindexlocal = -1;
	var feedforwardlayerindex = 0;
	feedForwardStepping();
	function feedForwardStepping() {
		if (ffphasestring == "Feed Forward") {

			if (neuronthatsredindex == -1) { neuronthatsredindex = 0; }

			if (feedforwardlayerindex == 0 && weightindexlocal == 199 && neuronthatsredindex < 783) { weightindexlocal = 0; neuronthatsredindex++; }
			else if (feedforwardlayerindex == 0 && weightindexlocal == 199 && neuronthatsredindex == 783) { neuronthatsredindex++; weightindexlocal = 0; feedforwardlayerindex = 1; offset[0] = 0; offset[1] = 0; }
			else if (feedforwardlayerindex == 1 && weightindexlocal == 9 && neuronthatsredindex < 784 + 199) { weightindexlocal = 0; neuronthatsredindex++; }
			else if (feedforwardlayerindex == 1 && weightindexlocal == 9 && neuronthatsredindex == 784 + 199) { weightindexlocal = -1; neuronthatsredindex = -1; feedforwardlayerindex = 0; offset[1] = 0; offset[2] = 0; ffphasestring = "Apply sigmoid"; }
			else { weightindexlocal++; }

			if (neuronthatsredindex != -1) {
				if (feedforwardlayerindex == 0) {
					if (Math.abs(neuronthatsredindex - offset[0]) >= 9) offset[0] = neuronthatsredindex;
					if (Math.abs(weightindexlocal - offset[1]) >= 9) offset[1] = weightindexlocal;
				} else if (feedforwardlayerindex == 1) {
					if (Math.abs(neuronthatsredindex - 784 - offset[1]) >= 9) offset[1] = neuronthatsredindex - 784;
					if (Math.abs(weightindexlocal - offset[2]) >= 9) offset[2] = weightindexlocal;
				}

				stepFeedForward(neuronsobject.allneuronsarray, neuronthatsredindex, weightindexlocal, feedforwardlayerindex);
				weightthatsredindex = weightindexlocal % 9;

				if (fastforwardthislayerflag) {
					fastForwardFeedForwardOfThisLayer(neuronsobject.allneuronsarray, neuronthatsredindex, weightindexlocal, feedforwardlayerindex);
					if (feedforwardlayerindex == 0) {
						neuronthatsredindex = 784;
						weightindexlocal = 0;
						offset[0] = 0;
						feedforwardlayerindex++;
					} else if (feedforwardlayerindex == 1) {
						weightindexlocal = -1; neuronthatsredindex = -1; feedforwardlayerindex = 0; offset[0] = 0; offset[1] = 0; offset[2] = 0; ffphasestring = "Apply sigmoid";
					}
					fastforwardthislayerflag = false;
				}
			}
		}
		if (running) setTimeout(feedForwardStepping, animationsteppingspeed);
	}


	var stepfeednetworkindex = 0;
	feedStepping();
	function feedStepping() {
		if (globalinputdatapointer != -1) {
			neuronthatsredindex = stepfeednetworkindex;
			if (Math.abs(neuronthatsredindex - offset[0]) == 9) offset[0] = neuronthatsredindex;
			stepFeedNetwork(neuronsobject.allneuronsarray, stepfeednetworkindex, globalinputdatapointer);
			stepfeednetworkindex++;
			if (fastforwardthislayerflag) {
				// determine how many neurons are left
				// run for loop for rest and set indexes to last neuron so animation can go to next phase like it just ran through
				for (var i = stepfeednetworkindex; i < 784; i++) {
					stepFeedNetwork(neuronsobject.allneuronsarray, stepfeednetworkindex, globalinputdatapointer);
					stepfeednetworkindex++;
				}
				animationsteppingspeed = 1;
				fastforwardthislayerflag = false;
			}
		}
		if (stepfeednetworkindex < 784) {
			if (running) setTimeout(feedStepping, animationsteppingspeed);
		} else {
			// reset
			ffphasestring = "Feed Forward";
			equationimage.src = "";
			//infoimage = -1;
			offset[0] = 0;
			neuronthatsredindex = -1;
			globalinputdatapointer = -1;
			stepfeednetworkindex = 0;
			if (running) setTimeout(feedStepping, animationsteppingspeed);
		}
	}


	feedforwarddrawloop();
	function feedforwarddrawloop() {
		// neurons and weights
		ctx.fillStyle = 'rgba(0, 0, 0, 1)';
		ctx.lineWidth = 1;
		ctx.clearRect(0, 0, canvas.width, canvas.height);

		for (var j = 0; j < 3; j++) {
			for (var i = Math.max(offset[j], 0); i < parseInt((canvas.height - topoffset - neuronoffset - bottomoffset) / (neurondiameter + neuronbetweenspace)) + offset[j]; i++) {
				if (i >= xnumberofneuronsperlayer[j]) break;

				var previouslayeroffset = xnumberofneuronsperlayer[j - 2];
				if (previouslayeroffset === undefined) previouslayeroffset = 0;

				var currentlayeroffset = xnumberofneuronsperlayer[j - 1];
				if (currentlayeroffset === undefined) currentlayeroffset = 0;

				// neuron
				ctx.beginPath();
				// if we are in the input layer tho, there needs to be a small box left to the neuron depicting its grayscale color
				if ((currentlayeroffset + previouslayeroffset) == 0) {
					ctx.strokeStyle = 'rgba(0, 0, 0, 1)';
					ctx.lineWidth = 1;
					ctx.strokeRect(xlayers[j] - neurondiameter * 1.5, neuronoffset + topoffset + (neurondiameter + neuronbetweenspace) + (i - offset[j]) * (neurondiameter + neuronbetweenspace) - neurondiameter / 4, neurondiameter / 2, neurondiameter / 2);
					var colgrayinprect = parseFloat(getNeuron(neuronsobject.allneuronsarray, previouslayeroffset + currentlayeroffset + i, biasandweightsvaluebuffer).slice(0, 5));
					ctx.fillStyle = 'rgba(0, 0, 0, ' + colgrayinprect + ')';
					ctx.fillRect(xlayers[j] - neurondiameter * 1.5, neuronoffset + topoffset + (neurondiameter + neuronbetweenspace) + (i - offset[j]) * (neurondiameter + neuronbetweenspace) - neurondiameter / 4, neurondiameter / 2, neurondiameter / 2);
				}

				ctx.fillStyle = 'rgba(226, 61, 195, 1)';
				ctx.fillText(getBias(neuronsobject.allneuronsarray, previouslayeroffset + currentlayeroffset + i, biasandweightsvaluebuffer).slice(0, 5), xlayers[j] - neurondiameter / 4, neuronoffset + topoffset + (neurondiameter + neuronbetweenspace) + (i - offset[j]) * (neurondiameter + neuronbetweenspace) + 2);
				ctx.fillStyle = 'rgba(237, 98, 98, 1)';
				ctx.fillText(getNeuron(neuronsobject.allneuronsarray, previouslayeroffset + currentlayeroffset + i, biasandweightsvaluebuffer).slice(0, 5), xlayers[j] - neurondiameter / 4, neuronoffset + topoffset + (neurondiameter + neuronbetweenspace) + (i - offset[j]) * (neurondiameter + neuronbetweenspace) + 12);
				ctx.fillStyle = 'rgba(0, 129, 42, 1)';
				ctx.fillText(i, xlayers[j] - neurondiameter / 4, neuronoffset + topoffset + (neurondiameter + neuronbetweenspace) + (i - offset[j]) * (neurondiameter + neuronbetweenspace) - 8);
				if (neuronthatsredindex == previouslayeroffset + currentlayeroffset + i) { ctx.lineWidth = 2; ctx.strokeStyle = 'rgba(255, 0, 0, 1)'; }
				else { ctx.lineWidth = 1; ctx.strokeStyle = 'rgba(0, 0, 0, 1)'; }
				if (784 + 200 + predictedneuron == previouslayeroffset + currentlayeroffset + i) { ctx.lineWidth = 2; ctx.strokeStyle = 'rgba(0, 255, 0, 1)'; }
				ctx.arc(xlayers[j], neuronoffset + topoffset + (neurondiameter + neuronbetweenspace) + (i - offset[j]) * (neurondiameter + neuronbetweenspace), neurondiameter / 2, 0, 2 * Math.PI);
				ctx.stroke();

				ctx.fillStyle = 'rgba(109, 173, 224, 1)';
				if (j < 2) {
					for (var k = 0 - (i - offset[j]); k < ((canvas.height - neuronoffset - bottomoffset) / (neurondiameter + neuronbetweenspace)) - 1 - (i - offset[j]); k++) {
						if ((i - offset[j] + k) > xnumberofneuronsperlayer[j + 1] - offset[j + 1] - 1) break;
						if (offset[j + 1] < 0 && (i - offset[j] + k) < Math.abs(offset[j + 1])) continue;

						if (neuronthatsredindex == previouslayeroffset + currentlayeroffset + i && weightthatsredindex == k + (i - offset[j])) {
							ctx.strokeStyle = 'rgba(255, 0, 0, 1)';
							ctx.lineWidth = 2;
						} else {
							ctx.strokeStyle = 'rgba(0, 0, 0, 1)';
							ctx.lineWidth = 1;
						}

						ctx.beginPath();
						// this weight goes from neuron i of layer j to neuron (i-offset[j]+k) of layer j+1
						//console.log("from neuron ", i ,"of layer ", j, "to neuron ", offset[j+1]+(i-offset[j]+k), "of layer ", j+1);
						ctx.fillText(getWeight(neuronsobject.allneuronsarray, currentlayeroffset + i, offset[j + 1] + (i - offset[j] + k), biasandweightsvaluebuffer), (xlayers[j] + neurondiameter / 2 + xlayers[j + 1] - neurondiameter / 2) / 2, (neuronoffset + topoffset + (neurondiameter + neuronbetweenspace) + (i - offset[j]) * (neurondiameter + neuronbetweenspace) + neuronoffset + topoffset + (neurondiameter + neuronbetweenspace) + (i - offset[j] + k) * (neurondiameter + neuronbetweenspace)) / 2);
						ctx.moveTo(xlayers[j] + neurondiameter / 2, neuronoffset + topoffset + (neurondiameter + neuronbetweenspace) + (i - offset[j]) * (neurondiameter + neuronbetweenspace));
						ctx.lineTo(xlayers[j + 1] - neurondiameter / 2, neuronoffset + topoffset + (neurondiameter + neuronbetweenspace) + (i - offset[j] + k) * (neurondiameter + neuronbetweenspace));
						ctx.stroke();
					}
				}
			}
		}

		ctx.lineWidth = 2;

		ctx.fillStyle = 'rgba(0, 0, 0, 0.3)';
		ctx.strokeStyle = 'rgba(0, 0, 0, 1)';
		// top ones
		for (var i = 0; i < 3; i++) {
			ctx.beginPath();
			ctx.moveTo(xlayers[i], topoffset + neuronoffset - navbuttonheight / 2);
			ctx.lineTo(xlayers[i] - navbuttonwidth / 2, topoffset + neuronoffset + navbuttonheight - navbuttonheight / 2);
			ctx.lineTo(xlayers[i] + navbuttonwidth / 2, topoffset + neuronoffset + navbuttonheight - navbuttonheight / 2);
			ctx.closePath();
			ctx.stroke();
			ctx.fill();
		}
		// bottom ones
		for (var i = 0; i < 3; i++) {
			ctx.beginPath();
			ctx.moveTo(xlayers[i], canvas.height + neuronoffset - bottomoffset - navbuttonheight / 2 + navbuttonheight);
			ctx.lineTo(xlayers[i] - navbuttonwidth / 2, canvas.height + neuronoffset - bottomoffset - navbuttonheight / 2);
			ctx.lineTo(xlayers[i] + navbuttonwidth / 2, canvas.height + neuronoffset - bottomoffset - navbuttonheight / 2);
			ctx.closePath();
			ctx.stroke();
			ctx.fill();
		}
		if (running) requestAnimationFrame(feedforwarddrawloop);
	}

}

// KATEX HOW TO DYNAMICALLY RENDER https://stackoverflow.com/questions/58040112/how-to-use-katex-auto-renderer-in-dynamically-changing-html
// https://stackoverflow.com/questions/12652769/rendering-html-elements-to-canvas

/* this is for... running the network in training mode. We are going to need this for backpropagation visualisation
var imgindex = 0;
var size
var testimages;

setTimeout(function(){
size = malloc(4);
testimages = getImages("assets/t10k-images.idx3-ubyte", size);
}, 1000);

//setTimeout(draw, 1500);
function draw(){
	// --- FEED TRAINING DATA ---
	
	var addr = 0;
	addr = addr | HEAPU8[testimages+3+4*imgindex];
	addr = addr << 8;
	addr = addr | HEAPU8[testimages+2+4*imgindex];
	addr = addr << 8;
	addr = addr | HEAPU8[testimages+1+4*imgindex];
	addr = addr << 8;
	addr = addr | HEAPU8[testimages+0+4*imgindex];
	console.log(imgindex, addr);
	
	// draw
	var testimagedata = ctx.createImageData(280, 280);
	var testimg = [];
	// blow up every pixel 10x10
	for(var i = 0; i < 28*28; i++){
		var col = HEAPU8[addr+i];
		for(var j = 0 ; j < 10*10; j++){
			var x = j%10 + (i*10)%280;
			var y = parseInt(j/10)*280 + (parseInt(i/28)*10)*280;
			testimg[x+y] = col;
		}
	}
	var testtindex = 3;
	for(var i = 0; i < 280*280; i++){
		testimagedata.data[testtindex] = testimg[i];
		testtindex += 4;
	}
	$("#multilayerperceptroncanvas4")[0].getContext("2d").putImageData(testimagedata, 0, 0);

	predict("assets/100_epoch_64_batchsize_relu_alpha0point5_784_200_10_95_13%.network", addr);
	imgindex++;
	//setTimeout(draw, 1000);
}
*/

/*
generally interesting code, but not to my use sadly -> rendering code SVG mixed with some vanilla HTML straight to canvas
const data = `
<svg xmlns='http://www.w3.org/2000/svg' width='200' height='200'>
  <foreignObject width='100%' height='100%'>
	<div xmlns='http://www.w3.org/1999/xhtml' style='font-size:40px'>
	  <em>I</em> like <span style='color:white; text-shadow:0 0 2px blue;'>CANVAS</span>
	  <p>asdasd</p>
	</div>
  </foreignObject>
</svg>
`;
const img = new Image();
const svg = new Blob([data], {type: "image/svg+xml;charset=utf-8"});
const url = URL.createObjectURL(svg);
img.onload = function() {
	console.log("ASDASD");
   $("#multilayerperceptroncanvas")[0].getContext("2d").drawImage(img, 0, 0);
  URL.revokeObjectURL(url);
};
img.src = url;
*/

