uniform sampler2D colormap ;
uniform float lp ;
uniform vec2  center ;
uniform float scale ;
uniform vec2 U ;
uniform vec2 V ;
varying vec2 xy ;

float lp_length(vec2 v, float lp) {
   // Note: pow(x,y) returns undefined results if x < 0
   return (lp > 12.0) ? max(abs(v.x), abs(v.y)) : pow(pow(abs(v.x),lp) + pow(abs(v.y),lp), 1.0 / lp) ;
}

void main() {
   vec2 d = xy - center ;
   float lU = length(U) ;
   float lV = length(V) ;
   vec2 v = vec2(lV*dot(d,normalize(U)), lU*dot(d,normalize(V))) ;
   float l = lp_length(v, 2.0 * lp) * scale ;
   gl_FragColor = texture2D(colormap, vec2(l,0.0)) ;
}
