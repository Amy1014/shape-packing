varying vec2 texCoords ;
uniform sampler2D depthTexture ;
uniform sampler2D colorTexture ;	 // 900 * 900
uniform sampler2D occlusionTexture ;	 // 900 * 900
uniform mat4 mat_proj_inv ; 		// matrice de projection inversee (state opengl)
uniform int ao ;			// boolean
uniform float size ;			// taille de la texture

vec4 new_filtering() ;
vec3 get_obj_coords(in vec2 where) ;

void main() {
	if(texture2D(depthTexture, texCoords).z < 1. && texture2D(depthTexture, texCoords).z > 0.0){
		/*if (ao == 1) {
			gl_FragColor = (texture2D(colorTexture, texCoords)) * new_filtering();
		}
		if (ao == 0) {
			//gl_FragColor.rgb = (texture2D(colorTexture, texCoords)).rgb;
			gl_FragColor = new_filtering();
		}*/
			gl_FragColor = new_filtering();
	}
	else gl_FragColor = vec4(1.0);	

	if(texture2D(depthTexture, texCoords).z < 1. && texture2D(depthTexture, texCoords).z > 0.0){
		/*if (ao == 1) {
			gl_FragColor = (texture2D(colorTexture, texCoords)) * new_filtering();
		}
		if (ao == 0) {
			//gl_FragColor.rgb = (texture2D(colorTexture, texCoords)).rgb;
			gl_FragColor = new_filtering();
		}*/
			gl_FragColor = new_filtering();
	}
	else gl_FragColor = vec4(1.0);	

	if(texture2D(depthTexture, texCoords).z < 1. && texture2D(depthTexture, texCoords).z > 0.0){
		/*if (ao == 1) {
			gl_FragColor = (texture2D(colorTexture, texCoords)) * new_filtering();
		}
		if (ao == 0) {
			//gl_FragColor.rgb = (texture2D(colorTexture, texCoords)).rgb;
			gl_FragColor = new_filtering();
		}*/
			gl_FragColor = new_filtering();
	}
	else gl_FragColor = vec4(1.0);	

	if(texture2D(depthTexture, texCoords).z < 1. && texture2D(depthTexture, texCoords).z > 0.0){
		/*if (ao == 1) {
			gl_FragColor = (texture2D(colorTexture, texCoords)) * new_filtering();
		}
		if (ao == 0) {
			//gl_FragColor.rgb = (texture2D(colorTexture, texCoords)).rgb;
			gl_FragColor = new_filtering();
		}*/
			gl_FragColor = new_filtering();
	}
	else gl_FragColor = vec4(1.0);	
}

//----------------------------------------------------
// renvoie les coordonnees du point dans le systeme de coord de l'objet (3D) [d'apres gluUnProject]
// RQ : la matrice de projection doit etre deja inversee
//----------------------------------------------------
vec3 get_obj_coords(in vec2 where){

	vec4 p = vec4(where, texture2D(depthTexture, where).x, 1.) ;	//point (pixel traite) en coordonnees homogenes

	// Map from range 0 to 1 to range -1 to 1 => on recentre	
	p.xyz = p.xyz * 2. - 1. ;

	p = mat_proj_inv * p;

	if (p.w != 0.){
		p.xyz /= p.w ;
	}

	return vec3(p) ;
} 


vec4 new_filtering() {

	vec2 c = (floor( size * texCoords - vec2(0.5)) + vec2(0.5) ) / size; 
	vec2 d = size * (texCoords - c) ;	
	vec4 color ;
	vec4 coeff  = vec4((1. - d.x)*(1. - d.y), (d.x)*(1. - d.y), (d.x)*(d.y), (1. - d.x)*(d.y)  ) ;
	vec3 p = get_obj_coords(texCoords) ;
	vec3 p1 = get_obj_coords(c);
	vec3 p2 = get_obj_coords(c + vec2(1. / size, 0.));
	vec3 p3 = get_obj_coords(c + vec2(1. / size, 1. / size));
	vec3 p4 = get_obj_coords(c + vec2(.0, 1. / size));

	vec2 c1 = vec2(0.); 
	vec2 c2 = vec2(1. / size, 0.);
	vec2 c3 = vec2(1. / size, 1. / size);
	vec2 c4 = vec2(.0, 1. / size);


	vec4 ratio = vec4( length(p - p1) / length(d - c1),length(p - p2) / length(d - c2), length(p - p3) / length(d - c3), length(p - p4) / length(d - c4) )  ;
	float min_ratio = min (ratio.x, min (ratio.y , min(ratio.z, ratio.w) ) ) ;

	if(ratio.x > 3. * min_ratio){
		coeff.x = 0. ;
	}
	if(ratio.y > 3. * min_ratio){
		coeff.y = 0. ;
	}
	if(ratio.z > 3. * min_ratio){
		coeff.z = 0. ;
	}
	if(ratio.w > 3. * min_ratio){
		coeff.w = 0. ;
	}


	if (coeff == vec4(0., 0., 0., 0.)) coeff = vec4(min_ratio);	//TODO


	coeff = coeff / ( coeff.x + coeff.y + coeff.z + coeff.w) ; 	// la somme des poids doit etre 1


	color = coeff.x * texture2D(occlusionTexture, c) + 
		coeff.y * texture2D(occlusionTexture, c + vec2(1. / size, 0.)) + 
		coeff.z * texture2D(occlusionTexture, c + vec2(1. / size, 1. / size)) + 
		coeff.w * texture2D(occlusionTexture, c + vec2(.0, 1. / size)) ;	

	return color ;
}

