#version 110

uniform samplerCube lighting_map ;
uniform samplerCube lighting_map_diff ;

uniform sampler3D diffuse_map ;
uniform sampler2D specular_map ;
uniform sampler3D normal_map ;

uniform mat4 lighting_mat ;
uniform float specular_factor ;
uniform float transmit_factor ;
uniform float fresnel_exp ;
uniform float eta ;
uniform float bump_factor ;

varying vec3 i_normal ;
varying vec3 i_view ;
varying vec3 uvf ;
varying vec3 tangent ;
varying vec3 binormal ;
varying float i_material ;

varying vec3 puvx1 ;
varying vec3 puvx2 ;
varying vec3 puvx3 ;
varying vec3 triangle_bary ;

// We could use NVidia's texture mode for that,
// but we prefer the more portable way.
vec3 unpack(vec3 bump) {
   return 2.0 * (bump - vec3(0.5, 0.5, 0.5)) ;
}

// Not used for the moment
vec3 rgbe_to_rgb(vec4 rgbe) {
    float scale = exp2(rgbe.a * 255.0 - 127.0);
    return scale * rgbe.rgb ;
}

// Used to lookup the lighting from the diffuse and 
// specular cubemaps.
vec3 lookup_lighting(samplerCube sampler, vec3 dir) {
    vec3 dir2 = (lighting_mat * vec4(dir,0)).xyz ;
    return textureCube(sampler, dir2).rgb ;
}

// Used for the blended lookup in material-space textures.
// The blending is used to smooth the singularities.
vec3 lookup_uvf(sampler3D cube, float I) {
   const float gamma = 1.5 ;
   vec3 color = (
               triangle_bary.x * texture3D(cube,vec3(puvx1.xy,I)) + 
               triangle_bary.y * texture3D(cube,vec3(puvx2.xy,I)) +
               triangle_bary.z * texture3D(cube,vec3(puvx3.xy,I))
              ).xyz ;
  return 2.0 * pow(color, vec3(gamma, gamma, gamma)) ;
}

float rgb_to_gray(vec3 C) {
   const vec3 weights = vec3(0.3, 0.59, 0.11) ;
   const float gamma = 0.5 ;
   return pow(dot(weights, C), gamma) ;
}

void main() {

    mat4 lighting_mat_t = lighting_mat ;

    vec3 lightVec = gl_LightSource[0].position.xyz ;

    vec3 normal = normalize(i_normal) ;
    normal = normalize(normal) ;

    // calculate diffuse component
    vec3 diffuse = lookup_lighting(lighting_map_diff, normal) ;

    // calculate specular component
    vec3 R = reflect(i_view, normal) ;    
    vec3 specular = lookup_lighting(lighting_map, R) ;

    // calculate transmission
    vec3 transmit =  lookup_lighting(lighting_map, refract(i_view, normal, eta)) ;
    float fresnel = clamp(-dot(normal, normalize(i_view)), 0.0, 1.0) ;
    fresnel = 1.0 - pow(1.0 - fresnel, fresnel_exp) ;
    
    //  combine diffuse, specular and transmitted contributions, and output final vertex color
    vec3 color = mix(diffuse, transmit, transmit_factor*fresnel) + specular_factor * specular ;

    gl_FragColor.rgb = lookup_uvf(diffuse_map, rgb_to_gray(color)) ;
    gl_FragColor.a = 1.0;
}
