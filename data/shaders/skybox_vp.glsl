varying vec3 cubecoord ;

void main()
{
    cubecoord = gl_Vertex.xyz ;
    gl_Position = ftransform();
}
