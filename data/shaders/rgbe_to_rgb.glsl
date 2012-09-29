
vec3 rgbe_to_rgb(vec4 rgbe) {
    float scale = exp2(rgbe.a * 255.0 - 127.0);
    return scale * rgbe.rgb ;
}
