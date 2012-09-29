
uniform vec3 eye_pos ;
varying vec3 i_normal ;
varying vec3 i_normal_diff ;
varying vec3 i_view ;
varying vec2 uv ;

void main()
{
    gl_Position = ftransform();
    i_normal = - (gl_NormalMatrix * gl_Normal) ;
    i_normal_diff = - (gl_NormalMatrix * gl_MultiTexCoord0.xyz) ;
    i_view = (gl_ModelViewMatrix * gl_Vertex).xyz  - eye_pos ;
    gl_FrontColor = gl_Color ;
    uv = gl_MultiTexCoord1.xy ;
}
