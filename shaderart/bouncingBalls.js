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
		#define MAX_DIST 100.0
		#define SURF_DIST 0.01
		#define PI 3.14


		float getDist(vec3 p){
			vec4 s = vec4(-1.3, 1.*sin((iTime+1.3*PI)*2.)+2., 4, 1);
			vec4 s2 = vec4(1.7, 1.*sin(iTime*2.)+2., 6, 1);
			
			float sphereDist = length(p-s.xyz) - s.w;
			float sphereDist2 = length(p-s2.xyz) - s2.w;
			float planeDist = p.y;
			
			float d = min(sphereDist2,min(sphereDist, planeDist));
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
			vec3 lightPos = vec3(0, 15, 6);
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
			
			vec3 ro = vec3(0,1,0);
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

	$("#ballsbouncingcanvas")[0].width = Math.pow(2, Math.floor(getBaseLog(2, $(window).width() * 0.88)));
	var gl = $("#ballsbouncingcanvas")[0].getContext("webgl");
	
	// local copy of library https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Operators/Spread_syntax
	const r3webgl = {...parentr3webgl};
	r3webgl.gl = gl;
	var program = r3webgl.createShaderProgram(vertexshadersource, fragmentshadersource);
	var uniformlocations = r3webgl.getUniformLocations(program, ["iResolution", "iMouse", "iTime"]);

	window.addEventListener("orientationchange", function() {
		toggleRunning();
		setTimeout(function(){
			var newwidth = Math.pow(2, Math.floor(getBaseLog(2, $(window).width() * 0.88)));
			$("#ballsbouncingcanvas")[0].width = newwidth;
			toggleRunning();
		}, 200);
	});

	$("#ballsbouncingButton").click(toggleRunning);
	function toggleRunning(){
		if(running){
			running = false;
			$("#ballsbouncingButton").html("Start");
		}else{
			running = true;
			$("#ballsbouncingButton").html("Stop");
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
