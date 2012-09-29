
uniform vec3 eye_pos ;
uniform float tex_mul ;

varying vec3 i_normal ;
varying vec3 i_view ;
varying vec3 tangent ;
varying vec3 binormal ;

varying vec3 puvx1 ;
varying vec3 puvx2 ;
varying vec3 puvx3 ;
varying vec3 triangle_bary ;

void main()
{
    gl_Position = ftransform();
    i_normal = (gl_NormalMatrix * gl_Normal);

// Note: we should use the 3x3 submatrix of ModelView instead
    tangent = gl_NormalMatrix * gl_MultiTexCoord0.xyz ;
    binormal = gl_NormalMatrix * gl_MultiTexCoord1.xyz ;

    i_view = (gl_ModelViewMatrix * gl_Vertex).xyz  - eye_pos ;

    gl_FrontColor = gl_Color ;

    puvx1 = vec3(tex_mul*gl_MultiTexCoord2.xy, gl_MultiTexCoord6.x) ;
    puvx2 = vec3(tex_mul*gl_MultiTexCoord3.xy, gl_MultiTexCoord6.x) ;
    puvx3 = vec3(tex_mul*gl_MultiTexCoord4.xy, gl_MultiTexCoord6.x) ;
    triangle_bary = gl_MultiTexCoord5.xyz ;
}
