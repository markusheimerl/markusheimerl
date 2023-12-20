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

		#define MAX_STEPS 100
		#define MAX_DIST 1000.0
		#define SURF_DIST 0.01
		#define PI 3.14
		
		mat2 Rot(float a) {
			float s = sin(a);
			float c = cos(a);
			return mat2(c, -s, s, c);
		}
		
		float dBox(vec3 p, vec3 s) {
			p = abs(p)-s;
			return length(max(p, 0.))+min(max(p.x, max(p.y, p.z)), 0.);
		}

		float getDist(vec3 p){
			float planeDist = p.y;

			// rotating box
			vec3 bp = p;
			bp -= vec3(0,.75,8);		// translation
			bp.xz *= Rot(iTime);		// rotation
			float rotate = dBox(bp, vec3(.75));
		
			float d = min(rotate, planeDist);
			return d;
		}

		float rayMarch(vec3 ro, vec3 rd){
			float dO = 0.;
			
			for(int i = 0; i < MAX_STEPS; i++){
				vec3 p = ro + rd * dO;
				float dS = getDist(p);
				dO += dS;
				if(dO > MAX_DIST || dS < SURF_DIST) break;
			}
			
			return dO;
		}
		
		vec3 getNormal(vec3 p){
			float d = getDist(p);
			vec2 e = vec2(.01, 0);
			
			vec3 n = d - vec3(getDist(p-e.xyy), getDist(p-e.yxy), getDist(p-e.yyx));
			
			return normalize(n);
		}

		float getLight(vec3 p){
			vec3 lightPos = vec3(5, 15, 6);
			lightPos.xz += vec2(sin(iTime), cos(iTime))*2.;
			vec3 l = normalize(lightPos-p);
			vec3 n = getNormal(p);
			
			float dif = clamp(dot(n, l), 0., 1.);
			
			float d = rayMarch(p+n*SURF_DIST*2., l);
			if(d < length(lightPos - p)) dif *= 0.1;
			
			return dif;
		}

		void mainImage(out vec4 fragColor, in vec2 fragCoord){

			// --- RESHAPE COORDINATE SYSTEM TO 0,0,0 FOR THE CENTER OF CANVAS ---
			vec2 uv = fragCoord - iResolution/2.;
			uv.x /= iResolution.x;
			uv.y /= iResolution.y;
			uv.x *= iResolution.x/iResolution.y;
			
			vec3 col = vec3(0);
			
			vec3 ro = vec3(0,2,0);
			vec3 rd = normalize(vec3(uv.x, uv.y, 1));

			float d = rayMarch(ro, rd);
			
			vec3 p = ro + rd * d;
			
			float dif = getLight(p);
			col = vec3(dif);
			
			fragColor = vec4(col, 1.0);
		}

		void main() {
		  mainImage(gl_FragColor, gl_FragCoord.xy);
		}
	`;

	var running = true;

	$("#shapesandcolorscanvas")[0].width = Math.pow(2, Math.floor(getBaseLog(2, $(window).width() * 0.88)));
	var gl = $("#shapesandcolorscanvas")[0].getContext("webgl");
	
	// local copy of library https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Operators/Spread_syntax
	const r3webgl = {...parentr3webgl};
	r3webgl.gl = gl;
	var program = r3webgl.createShaderProgram(vertexshadersource, fragmentshadersource);
	var uniformlocations = r3webgl.getUniformLocations(program, ["iResolution", "iMouse", "iTime"]);

	window.addEventListener("orientationchange", function() {
		toggleRunning();
		setTimeout(function(){
			var newwidth = Math.pow(2, Math.floor(getBaseLog(2, $(window).width() * 0.88)));
			$("#shapesandcolorscanvas")[0].width = newwidth;
			toggleRunning();
		}, 200);
	});

	$("#shapesandcolorsButton").click(toggleRunning);
	function toggleRunning(){
		if(running){
			running = false;
			$("#shapesandcolorsButton").html("Start");
		}else{
			running = true;
			$("#shapesandcolorsButton").html("Stop");
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
		gl.viewport(0, 0, gl.canvas.width, gl.canvas.height);
		gl.useProgram(program);
		
		render(0);
		function render(now){
			if(running){
				now *= 0.001;
				const elapsedTime = Math.min(now - then, 0.1);
				time += elapsedTime;
				then = now;

				gl.uniform2f(uniformlocations.iResolution, gl.canvas.width, gl.canvas.height);
				gl.uniform2f(uniformlocations.iMouse, mouseX, mouseY);
				gl.uniform1f(uniformlocations.iTime, time);
				
				gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);
			}
			requestAnimationFrame(render);
		}
	}
}
