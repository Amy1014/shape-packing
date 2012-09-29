varying vec2 texCoords ;

uniform sampler2D depthTexture ;
uniform float size ;			// taille de la texture


vec3 bilinear_filtering(in int n);


void main() {
	gl_FragColor.rgb = bilinear_filtering(3);
	gl_FragColor.a = 1. ;
}

vec3 bilinear_filtering(in int n){

/*
* 	Somme ( gris_i+di_j+dj * r(di,dj))
*	_______________________________________
*
*		Somme( r(di,dj) )
*
*
*
*
*			di^2 + dj^2 + a (gris_i+di_j+dj - gris_ij)^2
* r(di,dj) = exp - ________________________________________________
*					2 s	
*
*
*/


	float color ;
	float gris ;
	float gris_ij = texture2D(depthTexture, texCoords).r ;
	float r ;
	float sum_r = 0. ;
	float di, dj ;
	float dpas = (float(n) - 1.) / 2. ;
	dpas /= size ;
	float alpha = 50. ;
	float sigma = 1. ;

	for (di = -dpas; di < dpas + (0.1 / size); di += 1. / size ){

		for (dj = -dpas; dj < dpas + (0.1 / size); dj += 1. / size ){

			gris = texture2D(depthTexture, vec2(texCoords.x+di,texCoords.y+dj)).r ;
			float d = gris - gris_ij ;
			r = exp( -( di*di + dj *dj + alpha * d * d) / (2. * sigma)) ;
			sum_r += r ;

			color += gris * r;
		}
	}

	color /= sum_r ;
 	
	return vec3(color);
}
