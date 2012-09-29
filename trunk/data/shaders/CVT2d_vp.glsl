varying vec2 xy ;

void main()
{
    xy = gl_Vertex.xy;
    gl_Position = ftransform();
}
