
uniform vec2 pixel_size ;

varying vec2 texCoords ;
varying vec2 texCoords2 ;

void main()
{
   gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
   texCoords = gl_MultiTexCoord0.xy - 0.5 * pixel_size ;
   texCoords2 = gl_MultiTexCoord0.xy + 0.5 * pixel_size ;
}
