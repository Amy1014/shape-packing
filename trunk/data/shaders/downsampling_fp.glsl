
varying vec2 texCoords ;
uniform sampler2D depthTexture ;
uniform float size ;			// taille de la texture

vec4 downsample(in float nb) ;


void main(){	

	gl_FragColor =  downsample(9.)  ;
}

vec4 downsample(in float nb){

	float min = 2. ;
	float color = 1. ;
	float pas = (nb - 1.) / 2. ;

	for (float i = texCoords.x - pas/size ; i < texCoords.x + (pas+0.1)/size ; i += 1./size){
		for (float j = texCoords.y - pas/size ; j < texCoords.y + (pas+0.1)/size ; j += 1./size){
			color = texture2D(depthTexture, vec2(i, j)).r ;
			if ( color < min ) min = color ;
		}
	}

	return vec4(min, min, min, 1.) ;
}
