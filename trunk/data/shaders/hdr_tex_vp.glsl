
uniform vec3 eye_pos ;
varying vec3 i_normal ;
varying vec3 i_normal_diff ;
varying vec3 i_view ;
varying vec2 uv ;
varying vec3 tangent ;
varying vec3 binormal ;


void main()
{
    gl_Position = ftransform();
    i_normal = (gl_NormalMatrix * gl_Normal);
    i_normal_diff = (gl_NormalMatrix * gl_MultiTexCoord0.xyz);

// Note: we should use the 3x3 submatrix of ModelView instead
    tangent = gl_NormalMatrix * gl_MultiTexCoord2.xyz ;
    binormal = gl_NormalMatrix * gl_MultiTexCoord3.xyz ;

    i_view = (gl_ModelViewMatrix * gl_Vertex).xyz  - eye_pos ;

    uv = gl_MultiTexCoord1.xy ;
    gl_FrontColor = gl_Color ;
}
