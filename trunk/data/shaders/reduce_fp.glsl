/*
   High Dynamic Range Rendering - Using P Buffers
   Allen Sherrod
   Article for the Game Developers Magazine.
*/


varying vec2 texCoords ;
varying vec2 texCoords2 ;

uniform sampler2D fullTexture;

void main() {
   vec2 t00 = texCoords ;
   vec2 t01 = vec2(texCoords.x, texCoords2.y) ;
   vec2 t10 = vec2(texCoords2.x, texCoords.y) ;   
   vec2 t11 = texCoords2 ;
   
   vec3 c00 = texture2D(fullTexture, texCoords).rgb ;

   gl_FragColor = 0.25 * (
       texture2D(fullTexture, t00) +
       texture2D(fullTexture, t01) +
       texture2D(fullTexture, t10) +
       texture2D(fullTexture, t11)        
   ) ;

}
