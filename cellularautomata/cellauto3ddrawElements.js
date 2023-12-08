/*
* Copyright 2021 Markus Heimerl, OTH Regensburg
* Licensed under CC BY-NC 4.0
*
* ANY USE OF THIS SOFTWARE MUST COMPLY WITH THE
* CREATIVE COMMONS ATTRIBUTION-NONCOMMERCIAL 4.0 INTERNATIONAL LICENSE
* AGREEMENTS
*/

/*
This draws the cube using only 8 verticies. Is not faster than the version where cubes are drawn with 36 verticies. No idea why.
If you know email me please at: markus (AT) markusheimerl (DOT) com
*/

cellularautomata3d();
function cellularautomata3d() {

	var z_crunchfactor = 1;
	var cellularworldsize = 70;
	var rulefactor = 2;

	function createCellGridWireframe(cellularworldsize, randomize) {
		var cellgrid = [];
		for (var x = 0; x < cellularworldsize; x++) {
			cellgrid[x] = [];
			for (var y = 0; y < cellularworldsize; y++) {
				cellgrid[x][y] = [];
				for (var z = 0; z < cellularworldsize / z_crunchfactor; z++) {
					if (randomize && z < cellularworldsize / 3 && y < cellularworldsize / 3 && x < cellularworldsize / 3) cellgrid[x][y][z] = Math.random() >= 0.5 ? 1 : 0;
					else cellgrid[x][y][z] = 0;

					// mark world edges
					if (y == 0 && z == 0) cellgrid[x][y][z] = 1;
					if (x == 0 && z == 0) cellgrid[x][y][z] = 1;
					if (x == 0 && y == 0) cellgrid[x][y][z] = 1;

					if (y == cellularworldsize - 1 && z == cellularworldsize / z_crunchfactor - 1) cellgrid[x][y][z] = 1;
					if (x == cellularworldsize - 1 && z == cellularworldsize / z_crunchfactor - 1) cellgrid[x][y][z] = 1;
					if (x == cellularworldsize - 1 && y == cellularworldsize - 1) cellgrid[x][y][z] = 1;

					if (y == cellularworldsize - 1 && z == 0) cellgrid[x][y][z] = 1;
					if (x == cellularworldsize - 1 && z == 0) cellgrid[x][y][z] = 1;
					if (x == cellularworldsize - 1 && y == 0) cellgrid[x][y][z] = 1;

					if (y == 0 && z == cellularworldsize / z_crunchfactor - 1) cellgrid[x][y][z] = 1;
					if (x == 0 && z == cellularworldsize / z_crunchfactor - 1) cellgrid[x][y][z] = 1;
					if (x == 0 && y == cellularworldsize - 1) cellgrid[x][y][z] = 1;
				}
			}
		}
		return cellgrid;
	}

	var cellgrid = createCellGridWireframe(cellularworldsize, true);

	cellularAutomataLogic();
	function cellularAutomataLogic() {
		if (running) {

			var cellgridnextframe = createCellGridWireframe(cellularworldsize, false);

			for (var x = 2; x < cellularworldsize - 2; x++) {
				for (var y = 2; y < cellularworldsize - 2; y++) {
					for (var z = 2; z < cellularworldsize - 2; z++) {
						// get cell status
						var status = cellgrid[x][y][z];

						// get neighbor count
						var zm1 = z - 1 < 2 ? cellularworldsize - 3 : z - 1;
						var zp1 = z + 1 == cellularworldsize - 2 ? 2 : z + 1;
						var ym1 = y - 1 < 2 ? cellularworldsize - 3 : y - 1;
						var yp1 = y + 1 == cellularworldsize - 2 ? 2 : y + 1;
						var xm1 = x - 1 < 2 ? cellularworldsize - 3 : x - 1;
						var xp1 = x + 1 == cellularworldsize - 2 ? 2 : x + 1;

						// below
						var neighborcount = cellgrid[x][y][zm1] + cellgrid[x][yp1][zm1] + cellgrid[x][ym1][zm1] + cellgrid[xp1][y][zm1]
							+ cellgrid[xm1][y][zm1] + cellgrid[xp1][yp1][zm1] + cellgrid[xp1][ym1][zm1] + cellgrid[xm1][ym1][zm1] + cellgrid[xm1][yp1][zm1];

						// level
						neighborcount += cellgrid[x][yp1][z] + cellgrid[x][ym1][z] + cellgrid[xp1][y][z]
							+ cellgrid[xm1][y][z] + cellgrid[xp1][yp1][z] + cellgrid[xp1][ym1][z] + cellgrid[xm1][ym1][z] + cellgrid[xm1][yp1][z];

						// above
						neighborcount += cellgrid[x][y][zp1] + cellgrid[x][yp1][zp1] + cellgrid[x][ym1][zp1] + cellgrid[xp1][y][zp1]
							+ cellgrid[xm1][y][zp1] + cellgrid[xp1][yp1][zp1] + cellgrid[xp1][ym1][zp1] + cellgrid[xm1][ym1][zp1] + cellgrid[xm1][yp1][zp1];

						// if cell alive and neiborcount 2,3 or 8, it stays alive
						if (status == 1 && (/*neighborcount == 2*rulefactor ||*/ neighborcount == 26 || neighborcount == 10 * rulefactor)) {
							cellgridnextframe[x][y][z] = 1;
						}
						// if cell is dead and neighborcount is 3, it is born
						else if (status == 0 && (/*neighborcount == 5 ||*/ neighborcount == 4)) {
							cellgridnextframe[x][y][z] = 1;
						}
						// in all other cases the cell remains dead or dies (already initialized as 0)
					}
				}
			}

			cellgrid = cellgridnextframe;
		}
		setTimeout(cellularAutomataLogic, 100);
	}



	// --- DRAWING CODE ---
	const vertexshadersource = `
		precision lowp float;

		attribute vec4 vertexposition;
		attribute vec4 color;

		uniform mat4 modelmatrix;
		uniform mat4 projectionmatrix;
		uniform mat4 viewmatrix;

		varying vec4 o_color;

		void main(){
			o_color = color;
			gl_Position = projectionmatrix * viewmatrix * modelmatrix * vertexposition;
		}
	`;

	const fragmentshadersource = `
		precision lowp float;
		varying vec4 o_color;
		void main(){
			gl_FragColor = o_color;
		}
	`;

	// --- GET CANVAS CONTEXT AND SETUP KEY LISTENERS ---
	var running = true;

	var fscanvas = $("#webglstudycanvas")[0];
	fscanvas.width = Math.pow(2, Math.floor(getBaseLog(2, $(window).width() * 0.88)));
	window.addEventListener("orientationchange", function () {
		setTimeout(function () {
			var newwidth = Math.pow(2, Math.floor(getBaseLog(2, $(window).width() * 0.88)));
			fscanvas.width = newwidth;
		}, 200);
	});
	var gl = fscanvas.getContext("webgl");
	const r3webgl = { ...parentr3webgl };
	r3webgl.gl = gl;

	$("#webglstudytoggle").click(toggle);
	function toggle() { running ? (running = false, $("#webglstudytoggle").html("Start")) : (running = true, $("#webglstudytoggle").html("Stop")); }
	// https://developer.mozilla.org/en-US/docs/Web/API/Pointer_Lock_API

	// --- ADD EVENT LISTENERS ---
	fscanvas.onclick = function () {
		fscanvas.requestPointerLock();
	}

	var viewxz = 0;
	var viewy = 0;

	document.addEventListener("pointerlockchange", function () {
		if (document.pointerLockElement === fscanvas) {
			document.addEventListener("mousemove", updatePosition, false);
		} else {
			document.removeEventListener("mousemove", updatePosition, false);
		}
	}, false);

	function updatePosition(e) {
		viewxz -= e.movementX * 0.1;
		viewy = Math.min(90, Math.max(-90, viewy - e.movementY * 0.1));
	}

	var keys = {};
	$("#webglstudycanvas").keydown(function (event) {
		keys[event.key] = true;
	});
	$("#webglstudycanvas").keyup(function (event) {
		keys[event.key] = false;
	});
	// --- ---

	// --- MAKE SHADERS AND PROGRAM ---
	const program = r3webgl.createShaderProgram(vertexshadersource, fragmentshadersource);
	gl.useProgram(program);

	// --- GET ALL ATTRIBUTE AND UNIFORM LOCATIONS
	const attribLocations = r3webgl.getAttribLocations(program, ["vertexposition", "color"]);
	const uniformLocations = r3webgl.getUniformLocations(program, ["modelmatrix", "viewmatrix", "projectionmatrix"]);

	// --- INIT 3D ---
	r3webgl.init3D();

	// --- MAKE BUFFER WITH POSITION DATA ---
	var positions = [
		-1.0, -1.0, -1.0,
		1.0, -1.0, -1.0,
		1.0, 1.0, -1.0,
		-1.0, 1.0, -1.0,
		-1.0, -1.0, 1.0,
		1.0, -1.0, 1.0,
		1.0, 1.0, 1.0,
		-1.0, 1.0, 1.0
	];
	var vertexbuffer = r3webgl.createBuffer(gl.ARRAY_BUFFER, positions);

	// --- MAKE VERTEX INDICES ---
	var indices = [
		5, 4, 0, 1, 5, 0, 6, 5, 1, 2, 6, 1,
		7, 6, 2, 3, 7, 2, 4, 7, 3, 0, 4, 3,
		6, 7, 4, 5, 6, 4, 1, 0, 3, 2, 1, 3
	];

	var indicesbuffer = gl.createBuffer();
	gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, indicesbuffer);
	gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(indices), gl.STATIC_DRAW);

	// --- GET COLOR DATA ---
	var colors = [
		0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
		1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0, 1.0,
		1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0,
		0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0
	];

	const colorsbuffer = gl.createBuffer();
	gl.bindBuffer(gl.ARRAY_BUFFER, colorsbuffer);
	gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(colors), gl.STATIC_DRAW);


	// --- SETUP CAMERA ---
	var camerapos = [30.0, 30.0, -100.0];

	requestAnimationFrame(drawScene);
	toggle();

	var then = 0;
	function drawScene(now) {
		if (running) {
			gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
			gl.clearColor(0.0, 0.0, 0.0, 1.0);

			// convert requestanimationframe timestamp to seconds
			now *= 0.001;
			// subtract the previous time from the current time
			var deltaTime = now - then;
			// remember the current time for the next frame
			then = now;

			// --- SETUP PROJECTION MATRIX --- (MAKE EVERYTHING 3D)
			var projectionmatrix = m4.createPerspectiveMatrix(m4.degreeToRadians(46.0), gl.canvas.clientWidth / gl.canvas.clientHeight, 1, 200000);
			gl.uniformMatrix4fv(uniformLocations.projectionmatrix, false, projectionmatrix);


			// --- SETUP LOOKAT MATRIX ---
			var lookatmatrix = m4.createIdentityMatrix();
			lookatmatrix = m4.mult(m4.createTranslationMatrix(camerapos[0] + Math.sin(m4.degreeToRadians(viewxz)), camerapos[1] + Math.sin(m4.degreeToRadians(viewy)), camerapos[2] + Math.cos(m4.degreeToRadians(viewxz))), lookatmatrix);
			var lookatposition = [lookatmatrix[12], lookatmatrix[13], lookatmatrix[14]];


			// --- FIRST PERSON CAMERA ---
			var movementspeed = 0.8;
			var factorws = (keys["w"] ? 1 : keys["s"] ? -1 : 0);
			lookatposition[0] += Math.sin(m4.degreeToRadians(viewxz)) * movementspeed * factorws;
			lookatposition[1] += Math.sin(m4.degreeToRadians(viewy)) * movementspeed * factorws;
			lookatposition[2] += Math.cos(m4.degreeToRadians(viewxz)) * movementspeed * factorws;
			camerapos[0] += Math.sin(m4.degreeToRadians(viewxz)) * movementspeed * factorws;
			camerapos[1] += Math.sin(m4.degreeToRadians(viewy)) * movementspeed * factorws;
			camerapos[2] += Math.cos(m4.degreeToRadians(viewxz)) * movementspeed * factorws;

			var factorad = (keys["d"] ? 1 : keys["a"] ? -1 : 0);
			var movcamvector = m4.cross([Math.sin(m4.degreeToRadians(viewxz)), Math.sin(m4.degreeToRadians(viewy)), Math.cos(m4.degreeToRadians(viewxz))], [0, 1, 0]);
			lookatposition[0] += movcamvector[0] * movementspeed * factorad;
			lookatposition[2] += movcamvector[2] * movementspeed * factorad;
			camerapos[0] += movcamvector[0] * movementspeed * factorad;
			camerapos[2] += movcamvector[2] * movementspeed * factorad;

			var factoreq = (keys["e"] ? movementspeed : keys["q"] ? -movementspeed : 0);
			lookatposition[1] += factoreq;
			camerapos[1] += factoreq;


			// --- SETUP VIEWMATRIX --- (MOVE THE WORLD INVERSE OF THE CAMERAMOVEMENT)
			var cameramatrix = m4.lookAt(camerapos, lookatposition, [0, 1, 0]);
			var viewmatrix = m4.inverse(cameramatrix);
			var viewmatrixlocation = gl.getUniformLocation(program, "viewmatrix");
			gl.uniformMatrix4fv(uniformLocations.viewmatrix, false, viewmatrix);


			// --- CONNECT BUFFERS TO ATTRIBUTES --- (only has to be done once since the only object vertex data we ever need is that of a cube)
			r3webgl.connectBufferToAttribute(gl.ARRAY_BUFFER, vertexbuffer, attribLocations.vertexposition, 3, true);
			r3webgl.connectBufferToAttribute(gl.ARRAY_BUFFER, colorsbuffer, attribLocations.color, 4, true);


			// -- DRAW ---
			console.time("drawloop");
			for (var x = 0; x < cellularworldsize; x++) {
				for (var y = 0; y < cellularworldsize; y++) {
					for (var z = 0; z < cellularworldsize; z++) {
						if (cellgrid[x][y][z] == 1) {
							gl.uniformMatrix4fv(uniformLocations.modelmatrix, false, r3webgl.createModelMatrix(x, y, z, 0, 0, 0, 0.5, 0.5, 0.5));
							gl.drawElements(gl.TRIANGLES, 36, gl.UNSIGNED_SHORT, 0);
						}
					}
				}
			}
			console.timeEnd("drawloop");

		}
		requestAnimationFrame(drawScene);
	}
}
