varying vec3 xyz ;

void main()
{
    xyz = gl_Vertex.xyz;
    gl_Position = ftransform();
}
