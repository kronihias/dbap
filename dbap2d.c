/*
	a-dbap2d

 Pd Port by Matthias Kronlachner, 2013
 www.matthiaskronlachner.com
 
	distance based amplitude panning in 2d spatium
	
	a-dbap2d calculates distance coefficients for sound panning 
	in custom multi-channel loudspeaker configurations. 
	
	this code product, courtesy, andré sier,
	march,may,october,november, december 2006
	agosto07

	extended from trond lossius excelent tl.dbap
	thanks to trond for ideas in extending this object (bullseye method), 
		
	thanks:
		Mathieu Chamagne

	credits: 
		theory and max patch (Trond Lossius); 
		c coding max external (André Sier);
		windows port and several ideas like gains and inverse dbap (Mathieu Chamagne);
	
	license: this external is licenced by s373 (www.s373.net) in GNU-LGPL (see www.gnu.org)

	
	public availability: this external may be found at http://www.s373.net/code/dbap 
	
	todo:
		- local / global coords for multiplanar dimension expansion on the virtual world
		- reverberated sources 
		- voronoi diagrams for fields of audio regions
		


	modifies 4th build:
		- per sound gain settable with db's
		- expansion/restructuring of speakers to struct: accomodates more params' per woofer (like on/off)
		- boxes insepts


	modifies 3rd build:
		- list method compliant with <index> <pos1> <pos2> (optioal active status)
		- region expansion
		- per sound active support
		- various methods fixes with A_GIMME
		- bullseye method


	modifies 2nd build:
		- 2nd outlet now outlet anything with more info
		- struck sound element (for angle expansion, radius & falloff factors)
		- add more sonic inputs to be used as desired; via message or num elements list m3thod
		- add damping region when outside box (rect, oval)
		- no normalize algorithm that makes sharper distances (diffuse parameter is radius)
		- add set individual snd coords and individual output
		- array of om's to set different om's on different snd's

	a-dbap2d calculates distance coefficients for sound panning 
	in custom multi-channel loudspeaker configurations

	receives a list of pairs(2d) of sound positions 
	outputs matrix~ style input messages <i o v> (<input output value(0,1)>)
	
	(i-int, f-float, l-list)
	
	io:
		args : 1 : sets num sound inputs
		args : n : sets loudspeaker configuration

		messages (please update this):
			- list 				(should be <snd-index(0,n-1)>i <posx>f <posy>f (optional active status)i )
			- speakers<l> 		(n pairs of floating point values set n speakers position)
			- speakers<l> 		(n pairs of floating point values set n speakers position)
			- num <i> 			(set num sound inputs)
			- radius <f> 		(set radius for active sound inputs and output n active sounds values)
			- radius_snd <i,f> 	(set radius for sound input i and output sound i values)
			- set <i,f,f> 		(set position for sound input i and output sound i values)
			- om <l> 			(set output mode for sound input. 1 arg sets all, n args set n)
			- falloff_mode <i> 	(set falloff mode for all sound inputs in box modes)
			- falloff <f> 		(set falloff factor for all sound inputs in box modes)

*/

#define __myVERSION__ "0.1.6 osx universal binary"


// #include "ext.h"
// #include "ext_common.h"
#include     "m_pd.h"
#include    <stdlib.h>

#include <math.h>

#define DIMEN	2
#define OUTL	10

#define SETSYM SETSYMBOL

// #define CRT_SECURE_NO_DEPRECATE  // weird sprintf deprecating on windows

#define MAXVALUE 256
#define MAXBOXES 20

// 2D
#define dot(u,v) 		((u).x * (v).x + (u).y * (v).y )
#define	len(v)			sqrt(sqr((v).x)+sqr((v).y))
#define	len2(v)			(sqr((v).x)+sqr((v).y))
#define set(u,a,b)		(u).x = (a); (u).y = (b);
#define copy(u,v)		(u).x = (v).x; (u).y = (v).y;
#define sqr(u)			((u)*(u))
#define sub(u,v)		(u).x -= (v).x; (u).y -= (v).y;
#define add(u,v)		(u).x += (v).x; (u).y += (v).y;
#define addu(u,v)		(u).x += (v); (u).y += (v);
#define scale(u,v)		(u).x *= (v).x; (u).y *= (v).y;
#define scaleu(u,a)		(u).x *= (a); (u).y *= (a);

#define CLIP(u,a,b)		(((u)<(a))?(a):((u)>(b))?(b):(u))

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

typedef struct	_vec {	double x,y; 	} vec;

typedef struct	_dbox
{	vec min,max; 	
	int active;
	double gain;
} dbox;

typedef struct	_ssnd
{	
	vec		pos; //input snd coords	
	double	a; //angle 
	double	f; //falloff factor
	double	radius; 	
	int		om;
	int		active;
	double	gain;
} ssnd; // speaker sound element
//	int					om; 			// output mode: 0 normalize output coords
										// output mode: 1 dont normalize output coords 'sharper' than om 0
										// output mode: 2 use box clipping & normalize
										// output mode: 3 use box clipping & dont normalize


typedef struct	_speak	//speakers struct
{
	vec pos;
	int active;
} speak;


typedef struct	_dbap2d
{
	t_object 	         c_ob;		
	t_outlet 		        *c_out;
	t_outlet 		        *c_out2;	
	t_atom				 out[OUTL];				// output list

	ssnd				o[MAXVALUE];			// the snds and where and how they are
	int 				num_sounds;

	speak				woofers[MAXVALUE];		// the speakers
	int					num_woofers;

	dbox				boxes[MAXBOXES];		// the damping/boosting regions
	int					num_regions;

	int 	fall_off_mode; // linear, exponential, log .. not implemented
	
	vec		boxspeakersmin, boxspeakersmax, boxspeakersdim;				// the enclosing speakers box and info
	vec		boxcoordsmin, boxcoordsmax, boxcenter, boxdim, boxradius;
	int		box_mode;
	double	box_scale;

	int		bullseye;							// a toggle to exec bullseye or not
	int		regions;
  
} dbap2d;



t_symbol *ps_boxspeakers, *ps_boxcoords, *ps_boxcenter, *ps_boxdim, *ps_boxspeakersdim, *ps_boxradius;
t_symbol *ps_distmass, *ps_numspeakers, *ps_getspeakers, *ps_version, *ps_myversion;


static void *dbap2d_new (t_symbol *msg, long argc, t_atom *argv);
void dbap2d_free (dbap2d *x);
void dbap2d_bang (dbap2d *x);
void dbap2d_assist(dbap2d *x, void *b, long m, long a, char *s);

void dbap2d_list (dbap2d *x, t_symbol *msg, short argc, t_atom *argv);	// set all snd's coords methoid
void dbap2d_speakers (dbap2d *x, t_symbol *msg, short argc, t_atom *argv);	// speakers coords
void dbap2d_setspeaker (dbap2d *x, t_symbol *msg, short argc, t_atom *argv);	// speakers coords
void dbap2d_radius(dbap2d *x, t_symbol *msg, short argc, t_atom *argv);
void dbap2d_snd_radius(dbap2d *x, t_symbol *msg, short argc, t_atom *argv);
void dbap2d_snd_diffusion(dbap2d *x, t_symbol *msg, short argc, t_atom *argv); //diffusion tris to work on 0.1. range
void dbap2d_diffusion(dbap2d *x, t_symbol *msg, short argc, t_atom *argv); 
void dbap2d_gain (dbap2d *x, t_symbol *msg, short argc, t_atom *argv);
void dbap2d_snd_gain(dbap2d *x, t_symbol *msg, short argc, t_atom *argv);
void dbap2d_db (dbap2d *x, t_symbol *msg, short argc, t_atom *argv);
void dbap2d_snd_db(dbap2d *x, t_symbol *msg, short argc, t_atom *argv);

void dbap2d_snd_active(dbap2d *x, int snd, int active);
void dbap2d_speaker_active(dbap2d *x, int speaker, int active);
void dbap2d_region_active(dbap2d *x, int box, int active);

void dbap2d_box_gain(dbap2d *x, t_symbol *msg, short argc, t_atom *argv); 


void dbap2d_num_sounds(dbap2d *x, int numsnds);
void dbap2d_om(dbap2d *x, t_symbol *msg, short argc, t_atom *argv);

void dbap2d_falloff_mode(dbap2d *x, int numsnds);
void dbap2d_falloff(dbap2d *x, t_symbol *msg, short argc, t_atom *argv);
void dbap2d_snd_falloff(dbap2d *x, t_symbol *msg, short argc, t_atom *argv);

void dbap2d_setsndc(dbap2d *x, t_symbol *msg, short argc, t_atom *argv);
//static inline dbap2d_calc_snd(dbap2d *x, int snd);
void dbap2d_calc_snd(dbap2d *x, int snd);
void dbap2d_box_scale (dbap2d *x);
void dbap2d_setboxcoords (dbap2d *x, t_symbol *msg, short argc, t_atom *argv);
void dbap2d_setboxcenter (dbap2d *x, t_symbol *msg, short argc, t_atom *argv);
void dbap2d_setboxdim (dbap2d *x, t_symbol *msg, short argc, t_atom *argv);
void dbap2d_setboxscale (dbap2d *x, t_symbol *msg, short argc, t_atom *argv);
void dbap2d_setboxmode (dbap2d *x, int x1);


void dbap2d_region (dbap2d *x, t_symbol *msg, short argc, t_atom *argv);	
void dbap2d_regions (dbap2d *x, t_symbol *msg, short argc, t_atom *argv);	
void dbap2d_setregions (dbap2d *x, int x1);
void dbap2d_num_regions(dbap2d *x, int numregions);

void dbap2d_num_speakers(dbap2d *x, int numregions);


int dbap2d_check_bullseye(dbap2d *x, int snd);
void dbap2d_exec_bullseye(dbap2d *x, int snd, int speaker);
void dbap2d_setbullseye (dbap2d *x, int x1);



 double distan(vec o, vec p);
 double distan2(vec o, vec p);
 double dbap2d_box_check (dbap2d *x, vec coords, double radius, double falloff, double distmass) ;
 double dbap2d_box_check2 (dbap2d *x, vec coords, double radius, double falloff, double distmass) ;

void dbap2d_version (dbap2d *x, t_symbol *msg, short argc, t_atom *argv);
void dbap2d_getspeakers (dbap2d *x, t_symbol *msg, short argc, t_atom *argv);


void *dbap2d_class;


void dbap2d_free(dbap2d *x)
{
  //free(x->dbap2d);
}


static void *dbap2d_new (t_symbol *msg, long argc, t_atom *argv) //input the args
{
	 int i;

  dbap2d *x = (dbap2d *)pd_new(dbap2d_class);

  x->c_out = outlet_new(&x->c_ob, NULL);
  x->c_out2 = outlet_new(&x->c_ob, NULL);
  
	if(argc){
		if(argc==1) { //num inputs
			if (argv[0].a_type == A_FLOAT) { 
				x->num_sounds = argv[0].a_w.w_float; 
			} else { error("argument is not int"); }
		} else if(argc>1) { //speakers
			dbap2d_speakers(x, msg,argc,argv);
		}
	}

	// variables init
//	x->radius = 0.001; // globals deprecated
	x->num_woofers = 1;
	x->box_scale = 1;
	x->box_mode = 0; // box is sphere
//	x->fall_off = 0.1;
	
	x->bullseye = 0;

	for(i=0;i<MAXVALUE;i++)
	{
		set (x->o[i].pos, 0.,0.);
		x->o[i].a = 0.f;
		x->o[i].f = 0.0000001f;
		x->o[i].radius = 0.001;
		x->o[i].active = 1;
		x->o[i].om = 0;	//mathieu chamagne found that multiple a-dbaps instantions would remember om. this should fix it.
		x->o[i].gain = 1.0f;
	}

	for(i=0;i<MAXVALUE;i++)
	{
//		x->speakers[i] = 0.0f; 
		x->woofers[i].pos.x = 0.0f;
		x->woofers[i].pos.y = 0.0f;

		x->woofers[i].active = 1;
	}

	for(i=0;i<MAXBOXES;i++)
	{
		x->boxes[i].min.x = 0.0f; 	
		x->boxes[i].max.x = 0.0f; 
		x->boxes[i].min.y = 0.0f; 	
		x->boxes[i].max.y = 0.0f; 
		
		x->boxes[i].active = 0;
		
		x->boxes[i].gain = 1.0f;
	}


	 return(x);	
}

// radius is spatial blur
void dbap2d_radius(dbap2d *x, t_symbol *msg, short argc, t_atom *argv)
{
	double radius;
	int i;
  
	if(argv->a_type == A_FLOAT) {
		radius = (double)argv->a_w.w_float;	
	} else if(argv->a_type == A_FLOAT) {
		//radius = (double)argv->a_w.w_float;
	}
//	x->radius = radius;

	for(i=0;i<x->num_sounds;i++) {
	
		x->o[i].radius = radius;
		dbap2d_calc_snd(x,i);
		}
}

void dbap2d_snd_radius(dbap2d *x, t_symbol *msg, short argc, t_atom *argv)
{
	if (argc) {
		double rad;
//		int i;
		int id;

		if (argc==2) { //id rad
			id = (argv+0)->a_w.w_float;
			argv++;

		if(argv->a_type == A_FLOAT) {
			rad = (double)argv->a_w.w_float;	
		}


		x->o[id].radius=rad;
		dbap2d_calc_snd(x,id);

		}
		
	}



}

void dbap2d_snd_diffusion(dbap2d *x, t_symbol *msg, short argc, t_atom *argv) //diffusion tris to work on 0.1. range
{
	if (argc) {
		double rad;
//		int i;
		int id;

		if (argc==2) { //id rad
			id = (argv+0)->a_w.w_float;
			argv++;

		if(argv->a_type == A_FLOAT) {
			rad = (double)argv->a_w.w_float;	
		}
	


		x->o[id].radius=rad*MAX(x->boxspeakersdim.x, x->boxspeakersdim.y);
//		post("snd %ld diffusion is %f, like setting the radius at %f",id,rad,x->o[id].radius);
		dbap2d_calc_snd(x,id);

		}
		
	}



}

void dbap2d_diffusion(dbap2d *x, t_symbol *msg, short argc, t_atom *argv) //diffusion tris to work on 0.1. range
{
	if (argc) {
		double rad, maxdist;
		int i;
		
		if(argv->a_type == A_FLOAT) {
			rad = (double)argv->a_w.w_float;	
		}
	
		maxdist = MAX(x->boxspeakersdim.x, x->boxspeakersdim.y);

		for(i=0;i<x->num_sounds;i++) {
			x->o[i].radius=rad*maxdist;
			dbap2d_calc_snd(x,i);
		}
		
	}



}


void dbap2d_gain(dbap2d *x, t_symbol *msg, short argc, t_atom *argv)
{
	double gain;
	int i;
	
	if(argv->a_type == A_FLOAT) {
		gain = (double)argv->a_w.w_float;	
	}
	
	

	for(i=0;i<x->num_sounds;i++) {
	
		x->o[i].gain = gain;
		dbap2d_calc_snd(x,i);
		}
}

void dbap2d_snd_gain(dbap2d *x, t_symbol *msg, short argc, t_atom *argv)
{
	if (argc) {
		double gain;
//		int i;
		int id;

		if (argc==2) { //id rad
			id = (argv+0)->a_w.w_float;
			argv++;

		if(argv->a_type == A_FLOAT) {
			gain = (double)argv->a_w.w_float;	
		}
	


		x->o[id].gain=gain;
		dbap2d_calc_snd(x,id);

		}
		
	}



}


void dbap2d_db(dbap2d *x, t_symbol *msg, short argc, t_atom *argv)
{
	double gain;
	int i;
	
	if(argv->a_type == A_FLOAT) {
		gain = (double)argv->a_w.w_float;	
	}
	
	gain = pow(10.,(gain/20.)); //db to a

	for(i=0;i<x->num_sounds;i++) {
	
		x->o[i].gain = gain;
		dbap2d_calc_snd(x,i);
		}
}

void dbap2d_snd_db(dbap2d *x, t_symbol *msg, short argc, t_atom *argv)
{
	if (argc) {
		double gain;
//		int i;
		int id;

		if (argc==2) { //id rad
			id = (argv+0)->a_w.w_float;
			argv++;

		if(argv->a_type == A_FLOAT) {
			gain = (double)argv->a_w.w_float;	
		}
	
		gain = pow(10.,(gain/20.)); //db to a


		x->o[id].gain=gain;
		dbap2d_calc_snd(x,id);

		}
		
	}



}
void dbap2d_region_active(dbap2d *x, int box, int active)
{
	int s = MAX(0,box);
	x->boxes[s].active = active;
}

void dbap2d_speaker_active(dbap2d *x, int speaker, int active)
{
	int s = MAX(0,speaker);
	x->woofers[s].active = active;
}



void dbap2d_snd_active(dbap2d *x, int snd, int active)
{
	int s = MAX(0,snd);
	x->o[s].active = active;
}


void dbap2d_num_sounds(dbap2d *x, int numsnds)
{
	x->num_sounds = MAX(0,numsnds);
}

void dbap2d_num_regions(dbap2d *x, int numregions)
{
	x->num_regions = MAX(0,numregions);
}

void dbap2d_num_speakers(dbap2d *x, int numregions)
{
	x->num_woofers = MAX(0,numregions);
}


void dbap2d_om(dbap2d *x, t_symbol *msg, short argc, t_atom *argv)
{
//	x->om=CLIP(numsnds,0,4);

	int i,num;
	int temp[MAXVALUE];
	
	num = MAX (1,argc);
  
  if (num == 1) { // if just one argument -> set output mode for all sources!
    for(i=0;i<x->num_sounds;i++) {
      if (argv[i].a_type == A_FLOAT) { temp[0] = (int) argv[0].a_w.w_float; }
      
      x->o[i].om = temp[0];
    }
    // recalculate
    for(i=0;i<x->num_sounds;i++) {
      dbap2d_calc_snd(x, i);
    }
  } else {
    
    for(i=0;i<num;i++) {
      if (argv[i].a_type == A_FLOAT) { temp[i] = (int) argv[i].a_w.w_float; }
      
      x->o[i].om = temp[i];
    }
    // recalculate
    for(i=0;i<num;i++) {
      dbap2d_calc_snd(x, i);
    }
  }

}
void dbap2d_falloff_mode(dbap2d *x, int numsnds)
{
	x->fall_off_mode = numsnds;
}

void dbap2d_falloff(dbap2d *x, t_symbol *msg, short argc, t_atom *argv)
{
	int i;
	double x1;

		if(argv->a_type == A_FLOAT) {
			x1 = (double)argv->a_w.w_float;	
		}
	
	
//	x->fall_off = x1*0.00001;
  // post ("falloff set %f", x1);
	for(i=0;i<x->num_sounds;i++) {
		x->o[i].f = x1*0.00001;
		dbap2d_calc_snd(x,i);
		}
}

void dbap2d_snd_falloff(dbap2d *x, t_symbol *msg, short argc, t_atom *argv)
{
	int snd;
	double x1;

		if(argv->a_type == A_FLOAT) {
			snd = (int)argv->a_w.w_float;//	post("should be int. truncating to %ld", snd);
		}
		
		argv++;

		if(argv->a_type == A_FLOAT) {
			x1 = (double)argv->a_w.w_float;	
		}


		x->o[snd].f = x1*0.00001;
		dbap2d_calc_snd(x,snd);
}

void dbap2d_setboxscale (dbap2d *x, t_symbol *msg, short argc, t_atom *argv)
{
	double x1;

		if(argv->a_type == A_FLOAT) {
			x1 = (double)argv->a_w.w_float;	
		}
    else { post("feed me a float please"); return; };


//	post("setting box_scale to %f",x1);
	x->box_scale = x1;
	
	dbap2d_box_scale(x);

//	for(i=0;i<x->num_sounds;i++) dbap2d_calc_snd(x,i);
		

}


 double dbap2d_box_check (dbap2d *x, vec coords, double radius, double falloff, double distmass) 
{	// returns magnitude to scale final amp
	// return 1 if snd coords inside box & factor to scale final amp
	int inside=1;
	double d1, d2,d3;
	
	

	if(!x->box_mode) { // box is sphere

		d3 = distan(x->boxcenter,coords);
		d1 = MAX(x->boxradius.x,x->boxradius.y);
		d3-= d1;

		if(d3<0.) { // inside
			return (1.0f);	
		}
		
		//outside
		if(d3>0.) {
			d2 = d3 / d1;
			d2 = d2 + falloff * distmass * 0.1;
			(d2>1.)? 1. : (d2<0.)? 0. : d2; //CLIP(d2, 0., 1.);
			d2 = 1. - d2;
			if(d2<0.)
				d2 = 0.;
			
			
			return(d2);
		}
	
	
		
		return (1.0f); 
		
	} else if(x->box_mode) {

		if ( coords.x < x->boxcoordsmin.x || coords.x > x->boxcoordsmax.x ||	
			 coords.y < x->boxcoordsmin.y || coords.y > x->boxcoordsmax.y  )
			 	inside = 0;
		
		if (inside) { return (1.0f); }

	// find the less distance to the box
	d1 = distan(x->boxcoordsmin,coords);
	d2 = distan(x->boxcoordsmax,coords);
	d3 = distan(x->boxcenter,coords);
	
	
	return (1.0f);
	}

}

 double dbap2d_box_check2 (dbap2d *x, vec coords, double radius, double falloff, double distmass) 
{	// returns magnitude to scale final amp
	// return 1 if snd coords inside box & factor to scale final amp
	int inside=1;
	double d1, d2,d3;
	
	

	if(!x->box_mode) { // box is sphere

		d3 = distan2(x->boxcenter,coords);
		d1 = MAX(x->boxradius.x,x->boxradius.y);
		d3-= d1;

		if(d3<0.) { // inside
			return (1.0f);	
		}
		
		//outside
		if(d3>0.) {
			d2 = d3 / d1;
			d2 = d2 + falloff * distmass * 0.1;
			CLIP(d2, 0., 1.);
			d2 = 1. - d2;
			if(d2<0.)
				d2 = 0.;
			
			
			return(d2);
		}
	
		
		return (1.0f); 
		
	} else if(x->box_mode) {

		if ( coords.x < x->boxcoordsmin.x || coords.x > x->boxcoordsmax.x ||	
			 coords.y < x->boxcoordsmin.y || coords.y > x->boxcoordsmax.y  )
			 	inside = 0;
		
		if (inside) { 
		//	post("inside box");post("inside box"); 
			return (1.0f); 
		}

/*
	// find the less distance to the box
	d1 = distan(x->boxcoordsmin,coords);
	d2 = distan(x->boxcoordsmax,coords);
	d3 = distan(x->boxcenter,coords);
	
	post("outside box dist 2 box1: %f", d1);
	post("outside box dist 2 box2: %f", d2);
	post("outside box dist 3 box2: %f", d3);
*/	
	return (1.0f);
	}

}

void dbap2d_box_scale (dbap2d *x)
{
	int i;
	vec min,max,s;
	t_atom *out;
	vec center;
	
	out = x->out;

	 //init locals to first speaker
//	set (min, x->speakers[0], x->speakers[1]);
//	set (max, x->speakers[0], x->speakers[1]);
	set (min, x->woofers[0].pos.x, x->woofers[0].pos.y);
	set (max, x->woofers[0].pos.x, x->woofers[0].pos.y);

	for(i=0;i<x->num_woofers;i++) {
//		set (s, x->speakers[i*DIMEN+0],x->speakers[i*DIMEN+1]);
		set (s, x->woofers[i].pos.x, x->woofers[i].pos.y);
		// min max
		if(min.x > s.x) { min.x = s.x; }	
		if(min.y > s.y) { min.y = s.y; }	
		if(max.x < s.x) { max.x = s.x; }	
		if(max.y < s.y) { max.y = s.y; }	

	}


	copy (x->boxspeakersmin,min);
	copy (x->boxspeakersmax,max);

	
	// box center
	center.x = (max.x + min.x) * 0.5f; 
	center.y = (max.y + min.y) * 0.5f; 

	copy (x->boxcenter,center);


	// box speakers dim
	s.x = (max.x - min.x) ; 
	s.y = (max.y - min.y) ; 

	copy (x->boxspeakersdim,s);


	if(x->box_scale!=1.0f) {
		//move to local axis
		sub(min,center);
		sub(max,center);

		//scale
		scaleu(min,x->box_scale);
		scaleu(max,x->box_scale);
		
		//move back to global axis
		add(min,center);
		add(max,center);
	}

	copy(x->boxcoordsmin,min);
	copy(x->boxcoordsmax,max);

	// box dim
	s.x = (max.x - min.x) ; 
	s.y = (max.y - min.y) ; 

	copy (x->boxdim,s);


	// box radius
	s.x = x->boxdim.x*0.5 ; 
	s.y = x->boxdim.y*0.5 ; 

	copy (x->boxradius,s);
	
		
	SETFLOAT(out+0, x->boxspeakersmin.x);
	SETFLOAT(out+1, x->boxspeakersmin.y);
	SETFLOAT(out+2, x->boxspeakersmax.x);
	SETFLOAT(out+3, x->boxspeakersmax.y);
		
	outlet_anything(x->c_out2, ps_boxspeakers, 2*DIMEN, out);


	SETFLOAT(out+0, x->boxcoordsmin.x);
	SETFLOAT(out+1, x->boxcoordsmin.y);
	SETFLOAT(out+2, x->boxcoordsmax.x);
	SETFLOAT(out+3, x->boxcoordsmax.y);
		
	outlet_anything(x->c_out2, ps_boxcoords, 2*DIMEN, out);


	SETFLOAT(out+0, x->boxcenter.x);
	SETFLOAT(out+1, x->boxcenter.y);

	outlet_anything(x->c_out2, ps_boxcenter, DIMEN, out);

	SETFLOAT(out+0, x->boxdim.x);
	SETFLOAT(out+1, x->boxdim.y);

	outlet_anything(x->c_out2, ps_boxdim, DIMEN, out);

	SETFLOAT(out+0, x->boxspeakersdim.x);
	SETFLOAT(out+1, x->boxspeakersdim.y);

	outlet_anything(x->c_out2, ps_boxspeakersdim, DIMEN, out);

	SETFLOAT(out+0, x->boxradius.x);
	SETFLOAT(out+1, x->boxradius.y);

	outlet_anything(x->c_out2, ps_boxradius, DIMEN, out);
	

}

void dbap2d_setboxcoords (dbap2d *x, t_symbol *msg, short argc, t_atom *argv)
{
	double temp[4];
	int i;
	t_atom *out;

	out = x->out;
	
	for(i=0;i<4;i++) {
		if (argv[i].a_type == A_FLOAT) { temp[i] = (double) argv[i].a_w.w_float; }
		
	}
	
	x->boxspeakersmin.x = temp[0];
	x->boxspeakersmin.y = temp[1];
	x->boxspeakersmax.x = temp[2];
	x->boxspeakersmax.y = temp[3];

	SETFLOAT(out+0, x->boxspeakersmin.x);
	SETFLOAT(out+1, x->boxspeakersmin.y);
	SETFLOAT(out+2, x->boxspeakersmax.x);
	SETFLOAT(out+3, x->boxspeakersmax.y);
		
	outlet_anything(x->c_out2, ps_boxspeakers, 2*DIMEN, out);

}

void dbap2d_setboxcenter (dbap2d *x, t_symbol *msg, short argc, t_atom *argv)
{
	double temp[2];
	int i;
	t_atom *out=x->out;

	for(i=0;i<2;i++) {
		if (argv[i].a_type == A_FLOAT) { temp[i] = (double) argv[i].a_w.w_float; }

		
	}
	
	x->boxcenter.x = temp[0];
	x->boxcenter.y = temp[1];

	SETFLOAT(out+0, x->boxcenter.x);
	SETFLOAT(out+1, x->boxcenter.y);

	outlet_anything(x->c_out2, ps_boxcenter, DIMEN, out);


}

void dbap2d_setboxdim (dbap2d *x, t_symbol *msg, short argc, t_atom *argv)
{
	double temp[2];
	int i;
	t_atom *out = x->out;

	for(i=0;i<2;i++) {
		if (argv[i].a_type == A_FLOAT) { temp[i] = (double) argv[i].a_w.w_float; }

		
	}
	
	x->boxdim.x = temp[0];
	x->boxdim.y = temp[1];

	SETFLOAT(out+0, x->boxdim.x);
	SETFLOAT(out+1, x->boxdim.y);

	outlet_anything(x->c_out2, ps_boxdim, DIMEN, out);

}

void dbap2d_setboxmode (dbap2d *x, int x1)
{
	x->box_mode = x1;
}

void dbap2d_setbullseye (dbap2d *x, int x1)
{
	if(x1)
		x->bullseye=1;
	else
		x->bullseye=0;
}

void dbap2d_speakers (dbap2d *x, t_symbol *msg, short argc, t_atom *argv)
{
	int i;
	double temp[MAXVALUE];
//	double oneoverdimen = 1/DIMEN;
//	vec min,max;
	t_atom *out;
	
	out = x->out;

	if(argc) {
//		x->num_woofers = MIN(MAXSPEAKERSPERSOUND,argc/DIMEN);
		x->num_woofers = MAX(0,argc/DIMEN);

		for(i=0;i<argc;i++) {
			if (argv[i].a_type == A_FLOAT) { temp[i] = (double) argv[i].a_w.w_float; }


//			x->speakers[i] = temp[i];
		}

		for(i=0;i<x->num_woofers;i++) {
			x->woofers[i].pos.x = temp[i*DIMEN+0];
			x->woofers[i].pos.y = temp[i*DIMEN+1];
		}

		
		outlet_float(x->c_out2, x->num_woofers);		
		

		if(x->box_scale) { //make box coords regarding speakers position
			
			dbap2d_box_scale(x);

		}

	} else {
		error(" no arguments via speakers method");
	}
}


void dbap2d_setspeaker (dbap2d *x, t_symbol *msg, short argc, t_atom *argv) // per speaker pos
{
	if (argc) {
		double	temp[DIMEN]={0.0f,0.0f};
		int i;
		int id;

		if (argc==3) { //id posx posy
			id = (argv+0)->a_w.w_float;
			argv++;

			for(i=0;i<DIMEN;i++) {
				if(argv->a_type == A_FLOAT) {
					temp[i] = (double)argv->a_w.w_float;	
					argv++;
				}
			}


			x->woofers[id].pos.x = temp[+0];
			x->woofers[id].pos.y = temp[+1];

//		x->speakers[id*DIMEN+0] = temp[0];
//		x->speakers[id*DIMEN+1] = temp[1];

		}
		
	}


}


void dbap2d_list (dbap2d *x, t_symbol *msg, short argc, t_atom *argv) //updated for <index> <pos> <optional active
{



	if (argc) {
	
		double	temp[DIMEN]={0.0f,0.0f};
		int i;
		int id;
		int status;

		if (argc<3) { 
			post("dbap's list method requires: snd index(int) + posx + posy (floats) + (optional active status0/1)");
		}

		if (argc==4) { //id posx posy status
			status = (argv+3)->a_w.w_float;
			x->o[id].active = status;
		}
		
		if (argc==3) { //id posx posy
			id = (argv+0)->a_w.w_float;
			argv++;

			for(i=0;i<DIMEN;i++) {
				if(argv->a_type == A_FLOAT) {
					temp[i] = (double)argv->a_w.w_float;	
					argv++;
				}
			}



		set(x->o[id].pos, temp[0], temp[1]); 
		dbap2d_calc_snd(x, id);
		}
		
	}



/*	previous method

	int i,numtriplets,j;
	double temp[MAXVALUE];
	


		
	numtriplets = argc / DIMEN;
	x->num_sounds = numtriplets;
//	post("num sounds: %ld", x->num_sounds);
	
	for(i=0; i<argc;i++) {	
			if (argv->a_type == A_FLOAT) {
				temp[i] = (double)argv->a_w.w_float;
				argv++;
			} 
			else if (argv->a_type == A_LONG) {
				temp[i] = (double)argv->a_w.w_float;
				argv++;
			} 
	}

	for(i=0;i<numtriplets;i++) {
		set(x->o[i].pos, temp[i*DIMEN+0], temp[i*DIMEN+1]); 
		dbap2d_calc_snd(x, i);
	}

*/
//	dbap2d_bang(x);
}



void dbap2d_setregions (dbap2d *x, int x1)
{
	if(x1)
		x->regions=1;
	else
		x->regions=0;
}

void dbap2d_regions (dbap2d *x, t_symbol *msg, short argc, t_atom *argv)
{
	int i;
	double temp[MAXVALUE];
//	double oneoverdimen = 1/DIMEN;
//	vec min,max;
	t_atom *out;
	
	out = x->out;

	if(argc) {
//		x->num_woofers = MIN(MAXSPEAKERSPERSOUND,argc/DIMEN);
		x->num_regions = MAX(0,argc/DIMEN*2);

		for(i=0;i<argc;i++) {
			if (argv[i].a_type == A_FLOAT) { temp[i] = (double) argv[i].a_w.w_float; }
			

//			x->speakers[i] = temp[i];
		}

		for(i=0;i<x->num_regions;i++) {
			x->boxes[i].min.x = temp[i*DIMEN*2+0];
			x->boxes[i].min.y = temp[i*DIMEN*2+1];
			x->boxes[i].max.x = temp[i*DIMEN*2+2];
			x->boxes[i].max.y = temp[i*DIMEN*2+3];
		}

		
//		outlet_float(x->c_out2, x->num_woofers);		
		

	} else {
		error(" no arguments via speakers method");
	}
}


void dbap2d_region (dbap2d *x, t_symbol *msg, short argc, t_atom *argv) // per speaker pos
{
	if (argc) {
		double	temp[DIMEN]={0.0f,0.0f};
		int i;
		int id;

		if (argc==5) { //id minx miny maxx maxy
			id = (argv+0)->a_w.w_float;
			argv++;

			for(i=0;i<DIMEN*2;i++) {
				if(argv->a_type == A_FLOAT) {
					temp[i] = (double)argv->a_w.w_float;	
					argv++;
				}
			}


			x->boxes[id].min.x = temp[+0];
			x->boxes[id].min.y = temp[+1];
			x->boxes[id].max.x = temp[+2];
			x->boxes[id].max.y = temp[+3];

//		x->speakers[id*DIMEN+0] = temp[0];
//		x->speakers[id*DIMEN+1] = temp[1];

		}
		
	}


}






int dbap2d_check_bullseye(dbap2d *x, int snd)
{
	int i;
	vec speaker,sndpt,disttemp;
	double dist;

	set(sndpt, x->o[snd].pos.x,x->o[snd].pos.y);
	for(i=0; i < x->num_woofers; i++ ) {
		set(speaker, 
			x->woofers[i].pos.x,
			x->woofers[i].pos.y);
		
	// cityblock distance
		disttemp.x = speaker.x - sndpt.x;
		disttemp.y = speaker.y - sndpt.y;
		if(disttemp.x<0.) disttemp.x *= -1;
		if(disttemp.y<0.) disttemp.y *= -1;

		dist = disttemp.x + disttemp.y;
		
		if ( dist < x->o[snd].radius ) //cityblock distance is not eucleadean, but should work for the purposes
			return i;
	}
	
	return -1;
	
}

void dbap2d_exec_bullseye(dbap2d *x, int snd, int speaker)
{
	t_atom	*out;
	int i;
	
	out = x->out;

			// output 01 levels
			for (i=0; i<x->num_woofers;i++) 
			{
				SETFLOAT(out+0, snd);
				SETFLOAT(out+1, i);
				if(i==speaker) {
//				 SETFLOAT(out+2,1.);
				 SETFLOAT(out+2,x->o[snd].gain);
				} else {
				 SETFLOAT(out+2,0.);
				}
				
				outlet_list(x->c_out, 0L, 3, out);

			}
				//output -1 on distmass
			SETFLOAT(out+0, snd);
			SETFLOAT(out+1, -1);
			outlet_anything(x->c_out2, ps_distmass, 2, out);

}

void dbap2d_exec_silence(dbap2d *x, int snd)
{
	t_atom	*out;
	int i;
	
	out = x->out;

			// output 01 levels
			for (i=0; i<x->num_woofers;i++) 
			{
				SETFLOAT(out+0, snd);
				SETFLOAT(out+1, i);
				 SETFLOAT(out+2,0.);
				
				outlet_list(x->c_out, 0L, 3, out);

			}

}


/// main method per sound
//static inline dbap2d_calc_snd(dbap2d *x, int snd)
void dbap2d_calc_snd(dbap2d *x, int snd)
{
	int 	i;
	vec		speaker;
	vec		inpt;
	vec		avec;
	double 	distance;
	double  radius;
	double  temp[MAXVALUE];
	int 	num_woofers = x->num_woofers;
	double 	sum=0., one_over_sum, val;
	int		bingo=-1;
	
	double	falloff;
	double	boxfactor=1.0f;
	
	t_atom	*out;
	
	out = x->out;

	set(inpt, x->o[snd].pos.x,x->o[snd].pos.y);
	radius = x->o[snd].radius;
	falloff = x->o[snd].f;
	
	if(x->o[snd].active == 0) {
		dbap2d_exec_silence(x,snd); // if inactive, output zero's &return
		return; 
	}

	switch(x->o[snd].om) 
	{

		case 0:
		// normalize & no box (fixed according to tl.dbap_1_over_r
		
			// pass 0 - bullseye method
		
			if(x->bullseye) {
					bingo=dbap2d_check_bullseye(x,snd);
					if(bingo!=-1) {
						dbap2d_exec_bullseye(x,snd,bingo);
						return;
					}
			}

			radius *= radius; //square radius from now on

			// pass 1 - calculate squares of distances
						
			for (i=0; i<num_woofers;i++) 
			{
			 if(x->woofers[i].active) {
				set(speaker, 
					x->woofers[i].pos.x,
					x->woofers[i].pos.y);
				
				avec.x = speaker.x-inpt.x;
				avec.y = speaker.y-inpt.y;

				avec.x *= avec.x;
				avec.y *= avec.y;

				temp[i] = 1. / (avec.x + avec.y + radius); //adding squared radius/diffusion parameter
			 }	
			}

			// pass 2 - calculate levels

			sum = 0.0f;
			for (i=0; i<num_woofers;i++) 
			{
				if(x->woofers[i].active)	sum += temp[i];
			}

			val = 1.0f/sqrt(sum);
	
			for (i=0; i<num_woofers;i++) 
			{

				if(x->woofers[i].active)	temp[i] = sqrt(temp[i]) * val; 
					
			}


			
			SETFLOAT(out+0, snd);
			SETFLOAT(out+1, val);
			outlet_anything(x->c_out2, ps_distmass, 2, out);
			
			// output
			for (i=0; i<num_woofers;i++) 
			{
				if(x->woofers[i].active) 
				{
					SETFLOAT(out+0, snd);
					SETFLOAT(out+1, i);

					if(x->o[snd].gain==1.0f) SETFLOAT(out+2, temp[i]);
					else SETFLOAT(out+2, temp[i]*x->o[snd].gain);
					
					outlet_list(x->c_out, 0L, 3, out);
	
				} else if(!x->woofers[i].active)  { //output 0's on the inactive speaker
					SETFLOAT(out+0, snd);
					SETFLOAT(out+1, i);
					SETFLOAT(out+2, 0.);
					outlet_list(x->c_out, 0L, 3, out);
				}
			}
			break;


		case 1:

			// pass 0 - bullseye method
		
			if(x->bullseye) {
					bingo=dbap2d_check_bullseye(x,snd);
					if(bingo!=-1) {
						dbap2d_exec_bullseye(x,snd,bingo);
						return;
					}
			}



			// no normalize & no box
			for (i=0; i<num_woofers;i++) 
			{
			 if(x->woofers[i].active) 
			 {
				set(speaker, 
					x->woofers[i].pos.x,
					x->woofers[i].pos.y);
				
				distance = distan2(speaker,inpt);

				distance -= radius;

				distance = MAX(radius,distance);

				distance = sqr(distance);
				distance = 1.0f/distance;
				temp[i] = sqr(distance);

				sum += temp[i];
			 }
			}
			
			one_over_sum = 1.00f/(sum);
			SETFLOAT(out+0, snd);
			SETFLOAT(out+1, one_over_sum);
			outlet_anything(x->c_out2, ps_distmass, 2, out);
			
			// output
			for (i=0; i<num_woofers;i++) 
			{
				if(x->woofers[i].active) 
				{
					temp[i] = temp[i]*one_over_sum;
				//	temp[i] = pow(temp[i],0.5);
					SETFLOAT(out+0, snd);
					SETFLOAT(out+1, i);
					if(x->o[snd].gain==1.0f) SETFLOAT(out+2, temp[i]);
					else SETFLOAT(out+2, temp[i]*x->o[snd].gain);
					
					outlet_list(x->c_out, 0L, 3, out);

				} else if(!x->woofers[i].active)  { //output 0's on the inactive speaker
					SETFLOAT(out+0, snd);
					SETFLOAT(out+1, i);
					SETFLOAT(out+2, 0.);
					outlet_list(x->c_out, 0L, 3, out);
				}

			}

			break;

		case 2:
			// pass 0 - bullseye method
		
			if(x->bullseye) {
					bingo=dbap2d_check_bullseye(x,snd);
					if(bingo!=-1) {
						dbap2d_exec_bullseye(x,snd,bingo);
						return;
					}
			}


		// normalize &  box
	
			for (i=0; i<num_woofers;i++) 
			{
				if(x->woofers[i].active)
				{
				set(speaker, 
					x->woofers[i].pos.x,
					x->woofers[i].pos.y);
				
				distance = distan(speaker,inpt);
//				post("snd1 %ld speaker %ld distance %f", snd, i, distance);

				distance -= radius;
//				post("snd2 %ld speaker %ld dist-rad %f", snd, i, distance);

				distance = MAX(radius,distance);

				distance = sqr(distance);
				distance = 1.0f/distance;
				temp[i] = sqr(distance);

				sum += temp[i];
//				sum += distance;
				}
			}
			
			one_over_sum = 1.00f/(sum);
			SETFLOAT(out+0, snd);
			SETFLOAT(out+1, one_over_sum);
			outlet_anything(x->c_out2, ps_distmass, 2, out);

//			post("normalize&box inpt %f %f", inpt.x, inpt.y);
			boxfactor = dbap2d_box_check(x, inpt,radius,falloff, one_over_sum);
//			post("snd %ld boxfactor %f", snd, boxfactor);
//      post ("falloff %f", falloff);
			
			// output
			for (i=0; i<num_woofers;i++) 
			{
				if(x->woofers[i].active)
				{
					temp[i] = temp[i]*one_over_sum;
					temp[i] = pow(temp[i],0.5);
					temp[i]*= boxfactor;		// scale box factor
					SETFLOAT(out+0, snd);
					SETFLOAT(out+1, i);
					if(x->o[snd].gain==1.0f) SETFLOAT(out+2, temp[i]);
					else SETFLOAT(out+2, temp[i]*x->o[snd].gain);
					
					outlet_list(x->c_out, 0L, 3, out);

				} else if(!x->woofers[i].active)  { //output 0's on the inactive speaker
					SETFLOAT(out+0, snd);
					SETFLOAT(out+1, i);
					SETFLOAT(out+2, 0.);
					outlet_list(x->c_out, 0L, 3, out);
				}
			}
			break;


		case 3:

			// pass 0 - bullseye method
		
			if(x->bullseye) {
					bingo=dbap2d_check_bullseye(x,snd);
					if(bingo!=-1) {
						dbap2d_exec_bullseye(x,snd,bingo);
						return;
					}
			}


			// no normalize &  box
			for (i=0; i<num_woofers;i++) 
			{
				if(x->woofers[i].active)
				{
				set(speaker, 
					x->woofers[i].pos.x,
					x->woofers[i].pos.y);
				
				distance = distan2(speaker,inpt);

				distance -= radius;

				distance = MAX(radius,distance);

				distance = sqr(distance);
				distance = 1.0f/distance;
				temp[i] = sqr(distance);

				sum += temp[i];
				}
			}
			
			one_over_sum = 1.00f/(sum);
			SETFLOAT(out+0, snd);
			SETFLOAT(out+1, one_over_sum);
			outlet_anything(x->c_out2, ps_distmass, 2, out);

	//		post("normalize&box inpt %f %f", inpt.x, inpt.y);
			boxfactor = dbap2d_box_check2(x, inpt,radius,falloff, one_over_sum);
	//		post("snd %ld boxfactor %f", snd, boxfactor);
			
			// output
			for (i=0; i<num_woofers;i++) 
			{
				if(x->woofers[i].active)
				{
					temp[i] = temp[i]*one_over_sum;
				//	temp[i] = pow(temp[i],0.5);
					temp[i]*= boxfactor;		// scale box factor
					SETFLOAT(out+0, snd);
					SETFLOAT(out+1, i);
					if(x->o[snd].gain==1.0f) SETFLOAT(out+2, temp[i]);
					else SETFLOAT(out+2, temp[i]*x->o[snd].gain);
					
					outlet_list(x->c_out, 0L, 3, out);

				} else if(!x->woofers[i].active)  { //output 0's on the inactive speaker
					SETFLOAT(out+0, snd);
					SETFLOAT(out+1, i);
					SETFLOAT(out+2, 0.);
					outlet_list(x->c_out, 0L, 3, out);
				}
			}

			break;




			case 4:
		// normalize & no box (fixed according to tl.dbap_1_over_r -- non sqrt mode
		
			// pass 0 - bullseye method
		
			if(x->bullseye) {
					bingo=dbap2d_check_bullseye(x,snd);
					if(bingo!=-1) {
						dbap2d_exec_bullseye(x,snd,bingo);
						return;
					}
			}

			radius *= radius; //square radius from now on

			// pass 1 - calculate squares of distances
						
			for (i=0; i<num_woofers;i++) 
			{
				if(x->woofers[i].active)
				{
				set(speaker, 
					x->woofers[i].pos.x,
					x->woofers[i].pos.y);
				
				avec.x = speaker.x-inpt.x;
				avec.y = speaker.y-inpt.y;

				avec.x *= avec.x;
				avec.y *= avec.y;

				temp[i] = 1. / (avec.x + avec.y + radius); //adding squared radius/diffusion parameter
				}
			}

			// pass 2 - calculate levels

			sum = 0.0f;
			for (i=0; i<num_woofers;i++) 
			{
				if(x->woofers[i].active) sum += temp[i];
			}

	//		val = 1.0f/sqrt(sum);
			val = 1.0f/sum;

	
			for (i=0; i<num_woofers;i++) 
			{

				if(x->woofers[i].active) temp[i] = temp[i] * val; 
					
			}


			
			SETFLOAT(out+0, snd);
			SETFLOAT(out+1, val);
			outlet_anything(x->c_out2, ps_distmass, 2, out);
			
			// output
			for (i=0; i<num_woofers;i++) 
			{
				if(x->woofers[i].active) 
				{
					SETFLOAT(out+0, snd);
					SETFLOAT(out+1, i);

					if(x->o[snd].gain==1.0f) SETFLOAT(out+2, temp[i]);
					else SETFLOAT(out+2, temp[i]*x->o[snd].gain);
					
					outlet_list(x->c_out, 0L, 3, out);

				} else if(!x->woofers[i].active)  { //output 0's on the inactive speaker
					SETFLOAT(out+0, snd);
					SETFLOAT(out+1, i);
					SETFLOAT(out+2, 0.);
					outlet_list(x->c_out, 0L, 3, out);
				}
			}
			break;

	
	
	
	}


}



void dbap2d_setsndc(dbap2d *x, t_symbol *msg, short argc, t_atom *argv)
{
	
	if (argc) {
		double temp[DIMEN] = {0.,0.};
		int i;
		int id;

		if (argc==3) { //id posx posy
			id = (argv+0)->a_w.w_float;
			argv++;

			for(i=0;i<DIMEN;i++) {
				if(argv->a_type == A_FLOAT) {
					temp[i] = (double)argv->a_w.w_float;	
					argv++;
				}
			}


		x->o[id].pos.x=temp[0];
		x->o[id].pos.y=temp[1];

		dbap2d_calc_snd(x, id);
		}
		
	}



//	set(x->o[snd].pos, (double)x1, (double)x2);
//	
//	dbap2d_calc_snd(x, snd);

	
}

 double distan(vec o, vec p)
{
	vec dist;
	
	dist.x = p.x - o.x;
	dist.y = p.y - o.y;
	
	return(len(dist));
	
}

 double distan2(vec o, vec p)
{
	vec dist;
	
	dist.x = p.x - o.x;
	dist.y = p.y - o.y;
	
	return(len2(dist));
}





void dbap2d_bang (dbap2d *x)
{


//	dbap2d_calc_snd(x,0);
}


void dbap2d_getspeakers (dbap2d *x, t_symbol *msg, short argc, t_atom *argv)
{

	t_atom  out[MAXVALUE];
	int i;
	

	for(i=0;i<x->num_woofers;i++) {
		SETFLOAT(out+i*DIMEN+0,x->woofers[i].pos.x);
		SETFLOAT(out+i*DIMEN+1,x->woofers[i].pos.y);
	}
	outlet_anything(x->c_out2, ps_getspeakers, x->num_woofers*DIMEN, out);

}


void dbap2d_version (dbap2d *x, t_symbol *msg, short argc, t_atom *argv)
{
	t_atom * out;
	
	out = x->out;
	
	SETSYM(out+0,ps_myversion);
	outlet_anything(x->c_out2, ps_version, 1, out);

}

void dbap2d_assist(dbap2d *x, void *b, long m, long a, char *s)
{
    if (m==1) { post("list of coords, speakers"); }
    else if (m==2&&a==0) { post("(list) to matrix~"); }
    else if (m==2&&a==1) { post("info"); }
}

void dbap2d_setup (void)
{
  //long int tick = gettime();
  
  dbap2d_class = class_new(gensym("dbap2d"), (t_newmethod)dbap2d_new, (t_method)dbap2d_free, sizeof(dbap2d), 0, A_GIMME, 0);
  
  class_addmethod(dbap2d_class, (t_method)dbap2d_list,					gensym("list"),				A_GIMME,0);
  class_addmethod(dbap2d_class, (t_method)dbap2d_speakers, gensym("speakers"),			A_GIMME,0);
  class_addmethod(dbap2d_class, (t_method)dbap2d_getspeakers,gensym("getspeakers"),		A_GIMME,0);
  class_addmethod(dbap2d_class, (t_method)dbap2d_setspeaker,				gensym("speaker"),			A_GIMME,0);
  class_addmethod(dbap2d_class, (t_method)dbap2d_num_sounds,				gensym("num"),				A_DEFFLOAT,0);
  class_addmethod(dbap2d_class, (t_method)dbap2d_radius,					gensym("radius"),			A_GIMME,0);
  class_addmethod(dbap2d_class, (t_method)dbap2d_snd_radius,				gensym("radius_snd"),			A_GIMME,0);
  class_addmethod(dbap2d_class, (t_method)dbap2d_snd_radius,				gensym("radius_input"),			A_GIMME,0);
  class_addmethod(dbap2d_class, (t_method)dbap2d_diffusion,				gensym("diffusion"),		A_GIMME,0);
  class_addmethod(dbap2d_class, (t_method)dbap2d_snd_diffusion,			gensym("diffusion_snd"),		A_GIMME,0);
  class_addmethod(dbap2d_class, (t_method)dbap2d_snd_diffusion,			gensym("diffusion_input"),		A_GIMME,0);
  class_addmethod(dbap2d_class, (t_method)dbap2d_gain,					gensym("gain"),				A_GIMME,0);
  class_addmethod(dbap2d_class, (t_method)dbap2d_snd_gain,				gensym("gain_snd"),				A_GIMME,0);
  class_addmethod(dbap2d_class, (t_method)dbap2d_snd_gain,				gensym("gain_input"),			A_GIMME,0);
  class_addmethod(dbap2d_class, (t_method)dbap2d_db,						gensym("db"),				A_GIMME,0);
  class_addmethod(dbap2d_class, (t_method)dbap2d_snd_db,					gensym("db_snd"),				A_GIMME,0);
  class_addmethod(dbap2d_class, (t_method)dbap2d_snd_db,					gensym("db_input"),				A_GIMME,0);
  class_addmethod(dbap2d_class, (t_method)dbap2d_setsndc,				gensym("set"),				A_GIMME,0);
  class_addmethod(dbap2d_class, (t_method)dbap2d_om,						gensym("om"),				A_GIMME,0);
  class_addmethod(dbap2d_class, (t_method)dbap2d_falloff_mode,			gensym("falloff_mode"),		A_FLOAT,0);
  class_addmethod(dbap2d_class, (t_method)dbap2d_falloff,				gensym("falloff"),			A_GIMME,0);
  class_addmethod(dbap2d_class, (t_method)dbap2d_snd_falloff,			gensym("falloff_snd"),			A_GIMME,0);
  class_addmethod(dbap2d_class, (t_method)dbap2d_snd_falloff,			gensym("falloff_input"),		A_GIMME,0);
  
  class_addmethod(dbap2d_class, (t_method)dbap2d_setboxcoords,			gensym("boxcoords"),		A_GIMME,0);
  class_addmethod(dbap2d_class, (t_method)dbap2d_setboxcenter,			gensym("boxcenter"),		A_GIMME,0);
  class_addmethod(dbap2d_class, (t_method)dbap2d_setboxdim,				gensym("boxdim"),			A_GIMME,0);
  class_addmethod(dbap2d_class, (t_method)dbap2d_setboxscale,			gensym("boxscale"),			A_GIMME,0);
  class_addmethod(dbap2d_class, (t_method)dbap2d_setboxmode,				gensym("boxmode"),			A_FLOAT,0);
  class_addmethod(dbap2d_class, (t_method)dbap2d_setbullseye,			gensym("bullseye"),			A_FLOAT,0);
  
  class_addmethod(dbap2d_class, (t_method)dbap2d_snd_active,				gensym("active"),			A_DEFFLOAT,A_DEFFLOAT,0);
  class_addmethod(dbap2d_class, (t_method)dbap2d_snd_active,				gensym("active_snd"),			A_DEFFLOAT,A_DEFFLOAT,0);
  class_addmethod(dbap2d_class, (t_method)dbap2d_snd_active,				gensym("active_input"),			A_DEFFLOAT,A_DEFFLOAT,0);
  class_addmethod(dbap2d_class, (t_method)dbap2d_speaker_active,			gensym("active_speaker"),			A_DEFFLOAT,A_DEFFLOAT,0);
  class_addmethod(dbap2d_class, (t_method)dbap2d_region_active,				gensym("active_regions"),			A_DEFFLOAT,A_DEFFLOAT,0);
  
  class_addmethod(dbap2d_class, (t_method)dbap2d_setregions,				gensym("regions"),			A_FLOAT,0);
  class_addmethod(dbap2d_class, (t_method)dbap2d_regions,				gensym("regions_all"),			A_GIMME,0);
  class_addmethod(dbap2d_class, (t_method)dbap2d_region,					gensym("region"),			A_GIMME,0);
  
  class_addmethod(dbap2d_class, (t_method)dbap2d_num_sounds,				gensym("num_sounds"),				A_DEFFLOAT,0);
  class_addmethod(dbap2d_class, (t_method)dbap2d_num_regions,			gensym("num_regions"),				A_DEFFLOAT,0);
  class_addmethod(dbap2d_class, (t_method)dbap2d_num_speakers,			gensym("num_speakers"),				A_DEFFLOAT,0);
  
  
	
  class_addmethod(dbap2d_class, (t_method)dbap2d_assist,					gensym("assist"), 			A_CANT, 0);
  class_addmethod(dbap2d_class, (t_method)dbap2d_version,				gensym("version"),				A_GIMME,0);
  
  post("dbap2d by andre sier, inspired by trond lossius");
  
  ps_boxspeakers = gensym("boxspeakers");
  ps_boxcoords = gensym("boxcoords");
  ps_boxcenter =  gensym("boxcenter");
  ps_boxdim = gensym("boxdim");
  ps_boxspeakersdim = gensym("boxspeakersdim");
  ps_boxradius = gensym("boxradius");
  ps_distmass = gensym("distmass");
  ps_numspeakers = gensym("numspeakers");
  ps_getspeakers = gensym("speakers");
  ps_version = gensym("version");
  ps_myversion = gensym(__myVERSION__);
}

/*
old main method

		case 0:
		// normalize & no box
			for (i=0; i<num_woofers;i++) 
			{
				set(speaker, 
					x->speakers[i*DIMEN+0],
					x->speakers[i*DIMEN+1]);
				
				distance = distan(speaker,inpt);
//				post("snd1 %ld speaker %ld distance %f", snd, i, distance);

				distance -= radius;
//				post("snd2 %ld speaker %ld dist-rad %f", snd, i, distance);

				distance = MAX(radius,distance);

				distance = sqr(distance);
				distance = 1.0f/distance;
				temp[i] = sqr(distance);

				sum += temp[i];
//				sum += distance;
			}
			
			one_over_sum = 1.00f/sum;
			SETFLOAT(out+0, snd);
			SETFLOAT(out+1, one_over_sum);
			outlet_anything(x->c_out2, ps_distmass, 2, out);
			
			// output
			for (i=0; i<num_woofers;i++) 
			{
				temp[i] = temp[i]*one_over_sum;
				temp[i] = pow(temp[i],0.5);
				SETFLOAT(out+0, snd);
				SETFLOAT(out+1, i);
				SETFLOAT(out+2, temp[i]);
				
				outlet_list(x->c_out, 0L, 3, out);

			}
			break;



*/