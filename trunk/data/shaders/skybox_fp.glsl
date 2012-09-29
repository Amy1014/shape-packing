
uniform samplerCube light_map ;
varying vec3 cubecoord ;

vec3 rgbe_to_rgb(vec4 rgbe) {
    float scale = exp2(rgbe.a * 255.0 - 127.0);
    return scale * rgbe.rgb ;
}

vec3 lookup_lighting(samplerCube sampler, vec3 dir) {
    return textureCube(sampler, dir).rgb ;
}

void main() {
   gl_FragColor.rgb = lookup_lighting(light_map, cubecoord) ;
   gl_FragColor.a = 1.0 ;
}
