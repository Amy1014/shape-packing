uniform sampler2D colormap ;
uniform float lp ;
uniform vec3  center ;
uniform float scale ;
uniform vec3 U ;
uniform vec3 V ;
uniform vec3 W ;
varying vec3 xyz ;

float lp_length(vec3 v, float lp) {
   // Note: pow(x,y) returns undefined results if x < 0
   return (lp > 12.0) ? max(max(abs(v.x), abs(v.y)),abs(v.z)) : 
       pow(pow(abs(v.x),lp) + pow(abs(v.y),lp) + pow(abs(v.z),lp), 1.0 / lp) ;
}

void main() {
   vec3 d = xyz - center ;
   vec3 v = vec3(
   	dot(d,U), 
	dot(d,V),
	dot(d,W)	
   ) ;
   float l = lp_length(v, 2.0 * lp) * scale ;
   gl_FragColor = texture2D(colormap, vec2(l,0.0)) ;
}
