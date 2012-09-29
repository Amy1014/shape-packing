varying vec2 texCoords;

uniform sampler2D depthTexture ;	// 900 * 900
uniform sampler2D depthTexture_small ;
uniform sampler2D randomTexture ;	// 8 * 8
uniform mat4 mat_proj_inv ; 		// matrice de projection inversee (state opengl)
uniform float size ;			// largeur de la grosse texture
uniform float small ;			// largeur de la petite texture


bool dehors(in vec2 point);			// retourne vrai quand le pixel est en dehors du champs
float ambient_occlusion(in vec2 from);		//vretourne le coeff d'occlusion du fragment en cours de traitement
vec3 get_obj_coords(in vec2 where, in sampler2D texture) ;		// renvoie les coordonnées du point dans le systeme de coord de l'objet (3D)
float get_obj_z(in vec2 where, in sampler2D texture) ;		// renvoie la hauteur du point dans le systeme de coord de l'objet (3D) 

float horizon_angle(in vec2 from, in vec3 from3D, in vec2 dir, in vec3 normale) ;
vec2 horizon_point(in vec2 from, in vec2 dir) ;
vec3 normal(in vec2 where) ;
vec2 get_dir(in float alpha, mat2 base2D) ;
mat2 get_base(in vec3 normale, in float angle_normale) ;

const float PI = 3.141593;

void main() {

	if(texture2D(depthTexture, texCoords).z < 1. && texture2D(depthTexture, texCoords).z > 0.){
		gl_FragColor.rgb = vec3(ambient_occlusion(texCoords)) ;
	}
	else {
		gl_FragColor.rgb = vec3(1.0);
	}

	gl_FragColor.a = 1.0 ;
}

//----------------------------------------------------
// renvoie les coordonnées du point dans le systeme de coord de l'objet (3D) [d'après gluUnProject]
// RQ : la matrice de projection doit etre deja inversee
//----------------------------------------------------
vec3 get_obj_coords(in vec2 where, in sampler2D texture){

	vec4 p = vec4(where, texture2D(texture, where).x, 1.) ;	//point (pixel traité) en coordonnees homogenes

	// Map from range 0 to 1 to range -1 to 1 => on recentre	
	p.xyz = p.xyz * 2. - 1. ;

	p = mat_proj_inv * p;

	if (p.w != 0.){
		p.xyz /= p.w ;
	}

	return vec3(p) ;
} 

float get_obj_z(in vec2 where, in sampler2D texture){

	float depth = texture2D(texture, where).x ;
	
	//si on est sur le fond, on renvoie une valeur impossible
	if (depth >= 1.) return -10000. ;

	// Map from range 0 to 1 to range -1 to 1 => on recentre
	depth = depth * 2. - 1. ;
	float z = (depth * mat_proj_inv[2][2] + mat_proj_inv[3][2]) / ( depth * mat_proj_inv[2][3] + mat_proj_inv[3][3]) ; 
	return z ;
} 

//------------------------------------------------
// renvoie la normale au pixel, estimée d'après la carte de hauteur
//------------------------------------------------
vec3 normal(in vec2 where){

	float delta = 3. / size; 

	vec3 p = get_obj_coords(where, depthTexture);	//le point actuel
	vec3 px_p = get_obj_coords(where + vec2(delta, 0.), depthTexture);
	vec3 py_p = get_obj_coords(where + vec2(0., delta), depthTexture);
	vec3 px_t, py_t; 	//les tangentes
	
	px_t = px_p - p;
	py_t = py_p - p;

	return normalize( vec3( - px_t.z * py_t.y,
		     		- px_t.x * py_t.z,
		       		  px_t.x * py_t.y)) ;	// si on considere que px_t.y == 0 et py_t.x == 0
}

//-----------------------------------------------
// renvoie les 2 vecteurs axes de l'ellipse (plan tangent projeté)
//-----------------------------------------------
mat2 get_base(in vec3 normale, in float angle_normale){

	mat2 base2D ;
	vec2 x, y , n_p;
	float l ;

	n_p = vec2(normale) ;
	l = length(n_p) ;
	x = (n_p / l) * cos(angle_normale) ;		//length(x) = cos(angle_normale)
	y = vec2(-n_p.y, n_p.x) / l ;

	base2D = mat2(x,y) ;	
	return base2D ;
}

//-----------------------------------------------
// renvoie la direction de découpage dans la sphere en vraies coords [a transformer ensuite en coords de texture]
//-----------------------------------------------
vec2 get_dir(in float alpha, mat2 base2D){

	vec2 dir =  base2D * vec2(cos(alpha), sin(alpha));

	return normalize(dir) ;
}

//-----------------------------------------
// renvoie true si le point est en dehors de la texture
//-----------------------------------------

bool dehors(in vec2 point){
	if ( point.x > 1. || point.x < 0.|| point.y > 1. || point.y < 0.) return true;
	return false;
}



//-------------------------------------------------
// retourne le coeff d'occlusion du fragment en cours de traitement
//-------------------------------------------------

float ambient_occlusion(in vec2 from){	

	const int nb_slices = 16 ;
	float beetween_angle = 2.0 * PI / float(nb_slices) ;
	float actual_angle = 0.0 ;
	actual_angle = texture2D(randomTexture, from * 50.).r * 2. * PI  ;
	float occlusion_factor = 0.0 ;
	vec2 dir_vector = vec2(cos(actual_angle), sin(actual_angle)) ;
	vec3 from3D = get_obj_coords(from, depthTexture) ;

	for (int i=0; i < nb_slices; i++){

		float h_angle = horizon_angle(from, from3D, dir_vector, vec3(0., 0., 1.)) ;
		//float h_angle = horizon_angle(from, from3D, dir_vector, normale) ;
		actual_angle += beetween_angle ;
		dir_vector = vec2(cos(actual_angle), sin(actual_angle)) ;
		occlusion_factor += h_angle ;
	}

	occlusion_factor /= float(nb_slices) ;
	// on ajoute la contribution de chaque angle d'horizon en tant que pourcentage par rapport à PI/2
	occlusion_factor /= (PI/2.) ;
	return occlusion_factor ;

	/*
	vec3 normale = normal(from);
	float angle_normale = acos(normale.z) ;
	mat2 base = get_base(normale, angle_normale) ;
	dir_vector = get_dir(actual_angle, base);

	for (int i=0; i < nb_slices; i++){

		float h_angle = horizon_angle(from, from3D, dir_vector, normale) ;
		//h_angle = min(h_angle, PI/2.);	//TODO verifier qu'on ne peut pas avoir d'angle négatif
		clamp(h_angle, -PI/2., PI/2.);
		occlusion_factor += sin(h_angle) * sin(h_angle) ;
		actual_angle += beetween_angle ;
		dir_vector = vec2(cos(actual_angle), sin(actual_angle)) ;
	}

	occlusion_factor /= float(nb_slices) ;
	return occlusion_factor ;*/

}

//----------------------------------------
// retourne l'angle entre l'horizon et le vecteur passé en paramètres dans la direction demandée
// le vecteur doit etre normalisé
//----------------------------------------
float horizon_angle(in vec2 from, in vec3 from3D, in vec2 dir, in vec3 normale){

	vec3 horizon = get_obj_coords(horizon_point(from, dir), depthTexture);
	horizon = horizon - from3D ;	//vecteur d'horizon depuis le pixel de départ

	return acos ( dot(normale, horizon) / length(horizon) ) ;
}


//-------------------------------------------------
// retourne le point horizon pour ce fragment
// on recherche le point de l'horizon en utilisant les coord de 
// texture ainsi que la hauteur contenue dans la depth map
//-------------------------------------------------	

vec2 horizon_point(in vec2 from, in vec2 dir){

	//
	bool near = true ;

	// point en cours dans la depth map
	vec2 point = from ;

	//point de l'horizon
	vec2 point_h ;

	// stocke la (différence de) tangente la plus grande (x et y dans la depth map)
	float delta_max = -100000. ;
	
	// hauteur de référence (on veut une hauteur non pas par rapport à zero mais par rapport à la hauteur du point de depart) en coords objet
	//float height_ref = texture2D(depthTexture, from).z ;
	float height_ref = get_obj_z(from, depthTexture) ;

	// distance parcourue jusqu'a maintenant
	float length = 0.0 ;

	// largeur de la texture de profondeur utilisée pour les lookups
	float width = size ;

	// largeur du pas dans la depth map en nb de pixels 
	float pixel_step = 2. ;

	//
	vec2 dir_step = normalize(dir) * pixel_step ;

	// on calcule le pas dans la depth map la plus précise
	vec2 pas = normalize(dir) * pixel_step / size ;  // 1 pixel / pas <=> normalize(dir) / size

	int n_pas = 2 ;
	point += float(n_pas) * pas ;
	length = (float(n_pas) * pixel_step / size) ;

	// on recherche le point de l'horizon en utilisant les coord de 
	// texture ainsi que la hauteur contenue dans la depth map
	while ( !dehors(point) && length < 0.71/2. )	// 0.71 = sqrt(2)/2 pour couvrir toute la fenetre
	{
		//float height = texture2D(depthTexture, point).x ;
		float height ;
		if (near == true){
			height = get_obj_z(point, depthTexture) ;
			pas = dir_step / size ;
			//width = size ; // evident par initialisation !!
		}
		else {
			height = get_obj_z(point, depthTexture_small) ;	
			pas = dir_step / small ;
			width = small ;
		}		
		//length = (float(n_pas) * pixel_step / size) ;			
		float delta = (height - height_ref) / length;	// on rapporte sur la longueur pour utiliser comme un tangente

		if( delta > delta_max){
			delta_max = delta ;
			point_h = point ;
		}

		//if (length > 0.025) near = false; 


		point += pas ;
		//n_pas++ ;
		length += pixel_step / width ;
	}

	return point_h ;
}


