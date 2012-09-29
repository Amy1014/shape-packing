#version 110

uniform samplerCube lighting_map ;
uniform samplerCube lighting_map_diff ;

uniform sampler2D diffuse_map ;
uniform sampler2D normal_map ;

uniform mat4 lighting_mat ;

uniform float eta ;
uniform float transmit_factor ;

const float specular_factor = 0.1 ;
const float bump_weight = 1.5 ;
const float fresnel_exp = 2.0 ;

varying vec3 i_normal ;
varying vec3 i_view ;
varying vec2 uv ;
varying vec3 tangent ;
varying vec3 binormal ;

vec3 unpack(vec3 bump) {
   return 2.0 * (bump - vec3(0.5, 0.5, 0.5)) ;
}

vec3 rgbe_to_rgb(vec4 rgbe) {
    float scale = exp2(rgbe.a * 255.0 - 127.0);
    return scale * rgbe.rgb ;
}

vec3 lookup_lighting(samplerCube sampler, vec3 dir) {
    vec3 dir2 = (lighting_mat * vec4(dir,0)).xyz ;
    return textureCube(sampler, dir2).rgb ;
}

void main() {

    mat4 lighting_mat_t = lighting_mat ;

    vec3 lightVec = gl_LightSource[0].position.xyz ;

    vec3 normal = normalize(i_normal) ;
    vec3 bump = texture2D(normal_map, uv).xyz ;
    bump.z = bump_weight * bump.z ;
    bump = unpack(bump) ;

    bump = bump.x * normalize(tangent) + bump.y * normalize(binormal) + bump.z * normal ;
    normal = normal + 2.0 * bump ;

    // calculate diffuse component
    vec3 diffuse = lookup_lighting(lighting_map_diff, normal) ;
    diffuse = diffuse * texture2D(diffuse_map, uv).rgb ;

    // calculate specular component
    vec3 R = reflect(i_view, normal) ;    
    vec3 specular = lookup_lighting(lighting_map, R) ;

    // calculate transmission
    vec3 transmit =  lookup_lighting(lighting_map, refract(i_view, normal, eta)) ;
    float fresnel = clamp(-dot(normal, normalize(i_view)), 0.0, 1.0) ;
    fresnel = 1.0 - pow(1.0 - fresnel, fresnel_exp) ;

    //  combine diffuse and specular contributions and output final vertex color
    gl_FragColor.rgb = mix(diffuse, transmit, transmit_factor*fresnel) + specular_factor * specular ;
    gl_FragColor.a = 1.0;
}
