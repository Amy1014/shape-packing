
uniform vec3 eye_pos ;
uniform vec3 origin ;
uniform vec3 U ;
uniform vec3 V ;

varying vec3 i_normal ;
varying vec3 i_view ;
varying vec2 uv ;
varying vec3 tangent ;
varying vec3 binormal ;


void main()
{
    gl_Position = ftransform() ;
    i_normal = (gl_NormalMatrix * gl_Normal);

// Note: we should use the 3x3 submatrix of ModelView instead
    tangent = gl_NormalMatrix * U ;
    binormal = gl_NormalMatrix * V ;

    i_view = (gl_ModelViewMatrix * gl_Vertex).xyz  - eye_pos ;

    vec3 W = gl_Vertex.xyz - origin ; W = normalize(W) ;

    uv.x = 0.5 * (1.0 + dot(W, U)) ;
    uv.y = 0.5 * (1.0 + dot(W, V)) ;
    
    gl_FrontColor = gl_Color ;
}
