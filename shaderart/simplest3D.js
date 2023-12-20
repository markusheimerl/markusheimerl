/*
* Copyright 2021 Markus Heimerl, OTH Regensburg
* Licensed under CC BY-NC 4.0
*
* ANY USE OF THIS SOFTWARE MUST COMPLY WITH THE
* CREATIVE COMMONS ATTRIBUTION-NONCOMMERCIAL 4.0 INTERNATIONAL LICENSE
* AGREEMENTS
*/

function getBaseLog(x, y) {
  return Math.log(y) / Math.log(x);
}

simplest3DFragmentShader();
function simplest3DFragmentShader(){

		const vertexshadersource = `
			attribute vec2 a_vertexdata;
			void main(){
				gl_Position = vec4(a_vertexdata.xy, 0.0, 1.0);
			}
		`;
		
		const fragmentshadersource = `
    precision highp float;

    uniform vec2 iResolution;
    uniform vec2 iMouse;
    uniform float iTime;
	
	struct sphere{
	  vec3 pos;
	  float radius;
	};

	float hardtransition(float x, float threashhold){
		if(x < threashhold) return 0.;
		else return 1.;
	}

	float distanceToPoint(vec3 CI, vec3 CP){
		return length(cross(CI,CP)) / length(CI);
	}

	void mainImage(out vec4 fragColor, in vec2 fragCoord){

		vec2 uv = fragCoord - iResolution/2.; // 0,0,0 is in center, left end of screen is -iResolution.x,0,0, right end on screen is 0,0,iResolution.x etc.
		uv.x /= iResolution.x; // 0,0,0 is in center, left end of screen is -1,0,0, right end on screen is 0,0,1, top of the screen is 0,1,0, bottom of screen is 0,-1,0
		uv.y /= iResolution.y;
		uv.x *= iResolution.x/iResolution.y; // we wouldnt need this line if canvas was a square. Fixes ratio issues. In my cases most of time the canvas is 1024x512. So x is 2 times y. This is what this calculation comes up with: iResolution.x/iResolution.y = 2.0

		vec3 cam = vec3(0.,0.,-2.);
		vec3 i = vec3(uv.x, uv.y, 0.);
		sphere s;
		s.pos = vec3(2.*cos(iTime),0.,5.+2.*sin(iTime));
		s.radius = 0.2;
		
		vec3 CI = i-cam;
		vec3 CP = s.pos-cam;
		
		float distancetopoint = distanceToPoint(CI, CP);
		distancetopoint = smoothstep(s.radius, s.radius+.01, distanceToPoint(CI, CP));

		fragColor = vec4(vec3(1.-distancetopoint), 1.0);
	}

    void main() {
      mainImage(gl_FragColor, gl_FragCoord.xy);
    }
  `;

	var running = true;

	$("#simplest3Dcanvas")[0].width = Math.pow(2, Math.floor(getBaseLog(2, $(window).width() * 0.88)));
	var gl = $("#simplest3Dcanvas")[0].getContext("webgl");
	
	// local copy of library https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Operators/Spread_syntax
	const r3webgl = {...parentr3webgl};
	r3webgl.gl = gl;
	var program = r3webgl.createShaderProgram(vertexshadersource, fragmentshadersource);
	var uniformlocations = r3webgl.getUniformLocations(program, ["iResolution", "iMouse", "iTime"]);

	window.addEventListener("orientationchange", function() {
		toggleRunning();
		setTimeout(function(){
			var newwidth = Math.pow(2, Math.floor(getBaseLog(2, $(window).width() * 0.88)));
			$("#simplest3Dcanvas")[0].width = newwidth;
			toggleRunning();
		}, 200);
	});

	$("#simplest3DToggleButton").click(toggleRunning);
	function toggleRunning(){
		if(running){
			running = false;
			$("#simplest3DToggleButton").html("Start");
		}else{
			running = true;
			$("#simplest3DToggleButton").html("Stop");
		}
	}

	start();
	//setTimeout(toggleRunning, 200);
	function start(){

		const quadvert = [1,  1, -1,  1, 1, -1, -1, -1,];
		var quadbuffer = r3webgl.createBuffer(gl.ARRAY_BUFFER, quadvert);
		r3webgl.connectBufferToAttribute(gl.ARRAY_BUFFER, quadbuffer, 0, 2, true);

		var mouseX = 0;
		var mouseY = 0;
		var then = 0;
		var time = 0;
		
		render(0);
		function render(now){
			if(running){
				now *= 0.001;
				const elapsedTime = Math.min(now - then, 0.1);
				time += elapsedTime;
				then = now;

				gl.viewport(0, 0, gl.canvas.width, gl.canvas.height);
				gl.useProgram(program);
				gl.uniform2f(uniformlocations.iResolution, gl.canvas.width, gl.canvas.height);
				gl.uniform2f(uniformlocations.iMouse, mouseX, mouseY);
				gl.uniform1f(uniformlocations.iTime, time);
				
				gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);
			}
			requestAnimationFrame(render);
		}
	}
}
