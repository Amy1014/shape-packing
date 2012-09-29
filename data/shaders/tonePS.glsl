/*
   This file used the following reference as a starting point:
   High Dynamic Range Rendering - Using P Buffers
   Allen Sherrod
   Article for the Game Developers Magazine.
   Bruno Levy: Ported to FBO, added vignette, blur and unsharp masking.
*/


varying vec2 texCoords;

uniform sampler2D fullTexture ;
uniform sampler2D blurTexture ;
uniform sampler2D depthTexture ;

uniform float exposure;
uniform float blur_amount ;
uniform bool do_vignette ;
uniform bool do_shadows ;
uniform bool do_positive_shadows ;
uniform float gamma ;
uniform float shadows_gamma ;

void vignette(inout vec3 c, const vec2 win_bias) {
        // convert window coord to [-1, 1] range
        vec2 wpos = 2.0*(win_bias - vec2(0.5, 0.5));

        // calculate distance from origin
        float r = length(wpos.xy);
        r = 1.0 - smoothstep(0.8, 1.5, r);
        c = c * r;
}

float get_obj_z(in vec2 where, in sampler2D texture){
        float depth = texture2D(texture, where).x ;
        //si on est sur le fond, on renvoie une valeur impossible
        // if (depth >= 1.) return -10000. ;
        // Map from range 0 to 1 to range -1 to 1 => on recentre
        depth = depth * 2. - 1. ;
        float z = 
           (depth * gl_ProjectionMatrixInverse[2][2] + gl_ProjectionMatrixInverse[3][2]) / 
           (depth * gl_ProjectionMatrixInverse[2][3] + gl_ProjectionMatrixInverse[3][3]) ;
        return z ;
}


float equalize_shadow(float x) {
    float result = x ;
//  if(result > 1.0) { result = 1.0 ; }
    result = pow(x, shadows_gamma) ;
    return result ;
}

float unsharp_masking(vec2 uv) {
   float orig_depth = get_obj_z(uv, depthTexture) ;
   float smoothed_depth = get_obj_z(uv, blurTexture) ;
   float result = orig_depth - smoothed_depth ;
   if(result < 0.0) {
      result = -equalize_shadow(-result) ;
   } else {
      if(do_positive_shadows) {
         result = equalize_shadow(result) ;
      } else {
         result = 0.0 ;
      }
   }
   return result ;
}

void main() {

   vec3 scenePass = texture2D(fullTexture, texCoords).rgb ;

   if(do_shadows) {
       float depth = texture2D(depthTexture, texCoords).x ;
       if (depth >= 1.)  {
           gl_FragColor.rgb = scenePass ;
       } else {
           float shadow = -100.0 * exposure * blur_amount * unsharp_masking(texCoords) ;
             // if(shadow > 1.0) { shadow = 1.0 ; }
             scenePass = (1.0 - shadow) * scenePass  ;
             // scenePass = scenePass - shadow * vec3(1.0, 1.0, 1.0) ;
             // scenePass = (1.0 - shadow).xxx ;
       }
   } else {
       scenePass = mix(
          scenePass,
          texture2D(blurTexture, texCoords).rgb,
          blur_amount
      ).rgb ;
   }

   // Increase the color by the exposure then use a power
   // to convert to LDR range.
   scenePass = scenePass * exposure;
   scenePass =  pow(scenePass, vec3(gamma,gamma,gamma));
   if(do_vignette) {
      vignette(scenePass, texCoords) ;
   }
     
   gl_FragColor.rgb = scenePass ;
   gl_FragColor.a = 1.0 ;

}
