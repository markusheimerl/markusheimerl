/*
* Copyright 2021 Markus Heimerl, OTH Regensburg
* Licensed under CC BY-NC 4.0
*
* ANY USE OF THIS SOFTWARE MUST COMPLY WITH THE
* CREATIVE COMMONS ATTRIBUTION-NONCOMMERCIAL 4.0 INTERNATIONAL LICENSE
* AGREEMENTS
*/
cellularautomatonMulti();
function cellularautomatonMulti() {

	const vertexshadersource = `
			attribute vec2 a_vertexdata;
			void main(){
				gl_Position = vec4(a_vertexdata.xy, 0.0, 1.0);
			}
		`;
	const drawfragmentshadersource = `
			precision lowp float;
			uniform sampler2D texture;
			uniform vec2 dimensions;
			void main(){
				gl_FragColor = texture2D(texture, gl_FragCoord.xy/dimensions);
			}
		`;
	const logicfragmentshadersource = `
			precision lowp float;
			uniform sampler2D texture;
			uniform vec2 dimensions;

			int getValue(int x, int y){
				return (texture2D(texture, (gl_FragCoord.xy + vec2(x,y))/dimensions).g) > 0.0 ? 1 : 0;
			}

			void main(){
				int neighb = 0;

				for(int i = ` + $("#cellRuleInputMultiRadius").val() + `*-1; i < ` + $("#cellRuleInputMultiRadius").val() + `+1; i++){
					for(int j = ` + $("#cellRuleInputMultiRadius").val() + `*-1; j < ` + $("#cellRuleInputMultiRadius").val() + `+1; j++){
						if(!(i == 0 && j == 0)){
							neighb += getValue(i, j);
						}
					}
				}

				int current = getValue(0,0);
				if(` + parseRulestring($("#cellRuleInputMulti").val())[0] + `){gl_FragColor = vec4(0.7, 0.9, 0.1, 1.0);}
				else if(` + parseRulestring($("#cellRuleInputMulti").val())[1] + `){gl_FragColor = vec4(current, current, current, 1.0);}
				else{gl_FragColor = vec4(0.0, 0.0, 0.0, 1.0);}
			}
		`;


	var running = true;
	var scaleofonepixel = 1;

	$("#cellularautomatoncanvasmultineighbor")[0].width = Math.pow(2, Math.floor(getBaseLog(2, $(window).width() * 0.88)));
	var gl = $("#cellularautomatoncanvasmultineighbor")[0].getContext("webgl");

	// local copy of library https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Operators/Spread_syntax
	const r3webgl = { ...parentr3webgl };

	r3webgl.gl = gl;

	var programs = {
		drawprogram: r3webgl.createShaderProgram(vertexshadersource, drawfragmentshadersource),
		logicprogram: r3webgl.createShaderProgram(vertexshadersource, logicfragmentshadersource)
	};

	var drawuniformlocations = r3webgl.getUniformLocations(programs.drawprogram, ["dimensions"]);
	var logicuniformlocations = r3webgl.getUniformLocations(programs.logicprogram, ["dimensions"]);

	var textures = {
		front: createTexture(true),
		back: createTexture(false)
	}

	window.addEventListener("orientationchange", function () {
		toggleRunning();
		setTimeout(function () {
			var newwidth = Math.pow(2, Math.floor(getBaseLog(2, $(window).width() * 0.88)));
			$("#cellularautomatoncanvasmultineighbor")[0].width = newwidth;
			textures.front = createTexture(true);
			textures.back = createTexture(false);
			toggleRunning();
		}, 200);
	});

	function createTexture(noise) {
		var tex = gl.createTexture();
		gl.bindTexture(gl.TEXTURE_2D, tex);
		var data = [];
		if (noise) {
			for (var i = 0; i < gl.canvas.width * gl.canvas.height; i++) { var col = Math.round(Math.random() * Math.random() * 0.8) * 255; data.push(col, col, col); }
		} else {
			for (var i = 0; i < gl.canvas.width * gl.canvas.height; i++) { data.push(0, 0, 0); }
		}
		gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGB, gl.canvas.width / scaleofonepixel, gl.canvas.height / scaleofonepixel, 0, gl.RGB, gl.UNSIGNED_BYTE, new Uint8Array(data));
		gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
		gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
		return tex;
	}

	$("#cellToggleButtonApply").click(function () {
		const logicfragmentshadersource = `
			precision lowp float;
			uniform sampler2D texture;
			uniform vec2 dimensions;

			int getValue(int x, int y){
				return (texture2D(texture, (gl_FragCoord.xy + vec2(x,y))/dimensions).g) > 0.0 ? 1 : 0;
			}

			void main(){
				int neighb = 0;

				for(int i = ` + $("#cellRuleInputMultiRadius").val() + `*-1; i < ` + $("#cellRuleInputMultiRadius").val() + `+1; i++){
					for(int j = ` + $("#cellRuleInputMultiRadius").val() + `*-1; j < ` + $("#cellRuleInputMultiRadius").val() + `+1; j++){
						if(!(i == 0 && j == 0)){
							neighb += getValue(i, j);
						}
					}
				}
				int current = getValue(0,0);
				if(` + parseRulestring($("#cellRuleInputMulti").val())[0] + `){gl_FragColor = vec4(0.7, 0.9, 0.1, 1.0);}
				else if(` + parseRulestring($("#cellRuleInputMulti").val())[1] + `){gl_FragColor = vec4(current, current, current, 1.0);}
				else{gl_FragColor = vec4(0.0, 0.0, 0.0, 1.0);}
			}
		`;
		console.log(logicfragmentshadersource);
		programs.logicprogram = r3webgl.createShaderProgram(vertexshadersource, logicfragmentshadersource);
		drawuniformlocations = r3webgl.getUniformLocations(programs.drawprogram, ["dimensions"]);
		logicuniformlocations = r3webgl.getUniformLocations(programs.logicprogram, ["dimensions"]);
	});

	$("#cellToggleButtonMulti").click(toggleRunning);
	function toggleRunning() {
		if (running) {
			running = false;
			$("#cellToggleButtonMulti").html("Start");
		} else {
			running = true;
			$("#cellToggleButtonMulti").html("Stop");
		}
	}

	$("#cellRandButtonMulti").click(function () {
		textures.front = createTexture(true);
	});


	start();
	setTimeout(toggleRunning, 200);
	function start() {

		const quadvert = [1, 1, -1, 1, 1, -1, -1, -1,];
		var quadbuffer = r3webgl.createBuffer(gl.ARRAY_BUFFER, quadvert);
		r3webgl.connectBufferToAttribute(gl.ARRAY_BUFFER, quadbuffer, 0, 2, true);

		var fb = gl.createFramebuffer();


		mainloop();
		function mainloop() {
			if (running) {
				draw();
				play();
			}
			//requestAnimationFrame(mainloop);
			setTimeout(mainloop, 40);
		}

		function draw() {
			gl.bindFramebuffer(gl.FRAMEBUFFER, null);
			gl.bindTexture(gl.TEXTURE_2D, textures.front);
			gl.viewport(0, 0, gl.canvas.width, gl.canvas.height);
			gl.useProgram(programs.drawprogram);
			gl.uniform2f(drawuniformlocations.dimensions, gl.canvas.width, gl.canvas.height);
			gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);
		}

		function play() {
			gl.bindFramebuffer(gl.FRAMEBUFFER, fb);
			gl.viewport(0, 0, gl.canvas.width / scaleofonepixel, gl.canvas.height / scaleofonepixel);
			gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, textures.back, 0);
			gl.useProgram(programs.logicprogram);
			gl.uniform2f(logicuniformlocations.dimensions, gl.canvas.width / scaleofonepixel, gl.canvas.height / scaleofonepixel);
			gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);

			var swap = textures.back;
			textures.back = textures.front;
			textures.front = swap;
		}
	}
}
