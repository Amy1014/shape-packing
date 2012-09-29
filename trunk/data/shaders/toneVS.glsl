/*
   High Dynamic Range Rendering - Using P Buffers
   Allen Sherrod
   Article for the Game Developers Magazine.
*/


varying vec2 texCoords;

void main()
{
   gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
   texCoords = gl_MultiTexCoord0.xy;
}
