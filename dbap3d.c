/*
	a-dbap3d
 
 Pd Port by Matthias Kronlachner, 2013
 www.matthiaskronlachner.com
 
	distance based amplitude panning in 2d spatium
	courtesy, andrŽ sier, march,may 2006
	«
	extended from trond lossius excelent tl.dbap
	thanks to trond for the ideas in extending this object

	modifies 2nd build:
		- 2nd outlet now outlet anything with more info
		- struck sound element (for angle expansion, radius & falloff factors)
		- add more sonic inputs to be used as desired; via message or num elements list m3thod
		- add damping region when outside box (rect, oval)
		- no normalize algorithm that makes sharper distances (diffuse parameter is radius)
		- add set individual snd coords and individual output
		- array of om's to set different om's on different snd's

	a-dbap3d calculates distance coefficients for sound panning 
	in custom multi-channel loudspeaker configurations

	receives a list of pairs(2d) of sound positions 
	outputs matrix~ style input messages <i o v> (<input output value(0,1)>)
	
	(i-int, f-float, l-list)
	
	io:
		args : 1 : sets num sound inputs
		args : n : sets loudspeaker configuration

		messages:
			- list 				(n pairs of floating point values set n sound position output values)
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

#define DIMEN	3
#define OUTL	10

#define MAXVALUE 256
#define MAXBOXES 20
//#define MAXSPEAKERSPERSOUND 50

#define A_DEFLONG A_DEFFLOAT
#define A_LONG A_FLOAT

#define SETLONG SETFLOAT
#define SETSYM SETSYMBOL

#define outlet_int outlet_float

#define w_long w_float

// 3D
#define dot(u,v)ÊÊ 		((u).x * (v).x + (u).y * (v).y + (u).z * (v).z)
#define	len(v)			sqrt(sqr((v).x)+sqr((v).z)+sqr((v).z))
#define	len2(v)			(sqr((v).x)+sqr((v).y)+sqr((v).z))
#define set(u,a,b,c)	(u).x = (a); (u).y = (b); (u).z = (c);
#define copy(u,v)		(u).x = (v).x; (u).y = (v).y; (u).z = (v).z;
#define sqr(u)			((u)*(u))
#define sub(u,v)		(u).x -= (v).x; (u).y -= (v).y; (u).z -= (v).z;
#define add(u,v)		(u).x += (v).x; (u).y += (v).y; (u).z += (v).z;
#define scale(u,v)		(u).x *= (v).x; (u).y *= (v).y; (u).z *= (v).z;
#define scaleu(u,a)		(u).x *= (a); (u).y *= (a); (u).z *= (a);

#define CLIP(u,a,b)		(((u)<(a))?(a):((u)>(b))?(b):(u))


#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

typedef struct	{	double x,y,z; 	} vec;

typedef struct	_dbox
{	vec min,max; 	
	int active;
	double gain;
} dbox;



typedef struct	
{	
	vec		pos; //input snd coords	
	double	a; //angle 
	double	f; //falloff factor
	double	radius; 	
	int		om;
	int		active;
	double	gain;
} ssnd; // speaker sound element


typedef struct	_speak	//speakers struct
{
	vec pos;
	int active;
} speak;


typedef struct	
{
	t_object 	         c_ob;			
	t_outlet 		        *c_out;		
	t_outlet 		        *c_out2;	
	t_atom				 out[OUTL]; // output list

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
} dbap3d;


t_symbol *ps_boxspeakers, *ps_boxcoords, *ps_boxcenter, *ps_boxdim, *ps_boxspeakersdim, *ps_boxradius;
t_symbol *ps_distmass, *ps_numspeakers, *ps_getspeakers, *ps_version, *ps_myversion;



void *dbap3d_new (t_symbol *msg, short argc, t_atom *argv);
void dbap3d_bang (dbap3d *x);
void dbap3d_assist(dbap3d *x, void *b, long m, long a, char *s);

void dbap3d_list (dbap3d *x, t_symbol *msg, short argc, t_atom *argv);	// set all snd's coords methoid
void dbap3d_speakers (dbap3d *x, t_symbol *msg, short argc, t_atom *argv);	// speakers coords
void dbap3d_setspeaker (dbap3d *x, t_symbol *msg, short argc, t_atom *argv);	// speakers coords
void dbap3d_radius(dbap3d *x, t_symbol *msg, short argc, t_atom *argv);
void dbap3d_snd_radius(dbap3d *x, t_symbol *msg, short argc, t_atom *argv);
void dbap3d_snd_diffusion(dbap3d *x, t_symbol *msg, short argc, t_atom *argv); //diffusion tris to work on 0.1. range
void dbap3d_diffusion(dbap3d *x, t_symbol *msg, short argc, t_atom *argv); 
void dbap3d_gain (dbap3d *x, t_symbol *msg, short argc, t_atom *argv);
void dbap3d_snd_gain(dbap3d *x, t_symbol *msg, short argc, t_atom *argv);
void dbap3d_db (dbap3d *x, t_symbol *msg, short argc, t_atom *argv);
void dbap3d_snd_db(dbap3d *x, t_symbol *msg, short argc, t_atom *argv);

void dbap3d_snd_active(dbap3d *x, int snd, int active);
void dbap3d_speaker_active(dbap3d *x, int speaker, int active);
void dbap3d_region_active(dbap3d *x, int box, int active);

void dbap3d_box_gain(dbap3d *x, t_symbol *msg, short argc, t_atom *argv); 


void dbap3d_num_sounds(dbap3d *x, int numsnds);
void dbap3d_om(dbap3d *x, t_symbol *msg, short argc, t_atom *argv);

void dbap3d_falloff_mode(dbap3d *x, int numsnds);
void dbap3d_falloff(dbap3d *x, t_symbol *msg, short argc, t_atom *argv);
void dbap3d_snd_falloff(dbap3d *x, t_symbol *msg, short argc, t_atom *argv);

void dbap3d_setsndc(dbap3d *x, t_symbol *msg, short argc, t_atom *argv);
//static inline dbap3d_calc_snd(dbap3d *x, int snd);
void dbap3d_calc_snd(dbap3d *x, int snd);
void dbap3d_box_scale (dbap3d *x);
void dbap3d_setboxcoords (dbap3d *x, t_symbol *msg, short argc, t_atom *argv);
void dbap3d_setboxcenter (dbap3d *x, t_symbol *msg, short argc, t_atom *argv);
void dbap3d_setboxdim (dbap3d *x, t_symbol *msg, short argc, t_atom *argv);
void dbap3d_setboxscale (dbap3d *x, t_symbol *msg, short argc, t_atom *argv);
void dbap3d_setboxmode (dbap3d *x, int x1);


void dbap3d_region (dbap3d *x, t_symbol *msg, short argc, t_atom *argv);	
void dbap3d_regions (dbap3d *x, t_symbol *msg, short argc, t_atom *argv);	
void dbap3d_setregions (dbap3d *x, int x1);
void dbap3d_num_regions(dbap3d *x, int numregions);

void dbap3d_num_speakers(dbap3d *x, int numregions);


int dbap3d_check_bullseye(dbap3d *x, int snd);
void dbap3d_exec_bullseye(dbap3d *x, int snd, int speaker);
void dbap3d_setbullseye (dbap3d *x, int x1);



 double distan(vec o, vec p);
 double distan2(vec o, vec p);
 double dbap3d_box_check (dbap3d *x, vec coords, double radius, double falloff, double distmass) ;
 double dbap3d_box_check2 (dbap3d *x, vec coords, double radius, double falloff, double distmass) ;

void dbap3d_version (dbap3d *x, t_symbol *msg, short argc, t_atom *argv);
void dbap3d_getspeakers (dbap3d *x, t_symbol *msg, short argc, t_atom *argv);


void *dbap3d_class;

t_symbol *ps_nothing;

//

void dbap3d_free(dbap3d *x)
{
  
}



void *dbap3d_new (t_symbol *msg, short argc, t_atom *argv) //input the args 
{
  int i,j;

  dbap3d *x = (dbap3d *)pd_new(dbap3d_class);

  x->c_out = outlet_new(&x->c_ob, NULL);
  x->c_out2 = outlet_new(&x->c_ob, NULL);
  
	if(argc){
		if(argc==1) { //num inputs
			if (argv[0].a_type == A_LONG) { 
				x->num_sounds = argv[0].a_w.w_long; 
			} else { error("argument is not int"); }
		} else if(argc>1) { //speakers
			dbap3d_speakers(x, msg,argc,argv);
		}
	}

	// variables init
//	x->radius = 0.001;
	x->num_woofers = 1;
	x->box_scale = 1;
	x->box_mode = 0;
//	x->fall_off = 0.1;

	for(i=0;i<MAXVALUE;i++)
	{
		set (x->o[i].pos, 0.,0.,0.);
		x->o[i].a = 0.f;
		x->o[i].f = 0.0000001f;
		x->o[i].radius = 0.001;
		x->o[i].active = 1;
		x->o[i].om = 0;	//mathieu chamagne found that multiple a-dbaps instantions would remember om. this should fix it.
		x->o[i].gain = 1.0f;
	}

	for(i=0;i<MAXVALUE;i++)
	{
//		x->speakerlist[i] = 0.0f; 
		x->woofers[i].pos.x = 0.0f;
		x->woofers[i].pos.y = 0.0f;

		x->woofers[i].active = 1;
	}


	 return(x);	
}

void dbap3d_radius(dbap3d *x, t_symbol *msg, short argc, t_atom *argv)
{
	double radius;
	int i;
	
	if(argv->a_type == A_FLOAT) {
		radius = (double)argv->a_w.w_float;	
	}
	
	
//	x->radius = radius;

	for(i=0;i<x->num_sounds;i++) {
	
		x->o[i].radius = radius;
		dbap3d_calc_snd(x,i);
		}
}

void dbap3d_snd_radius(dbap3d *x, t_symbol *msg, short argc, t_atom *argv)
{
	if (argc) {
		double rad;
//		int i;
		int id;

		if (argc==2) { //id rad
			id = (argv+0)->a_w.w_long;
			argv++;

		if(argv->a_type == A_FLOAT) {
			rad = (double)argv->a_w.w_float;	
		}
	


		x->o[id].radius=rad;
		dbap3d_calc_snd(x,id);

		}
		
	}

}

void dbap3d_snd_diffusion(dbap3d *x, t_symbol *msg, short argc, t_atom *argv) //diffusion tris to work on 0.1. range
{
	if (argc) {
		double rad;
//		int i;
		int id;

		if (argc==2) { //id rad
			id = (argv+0)->a_w.w_long;
			argv++;

		if(argv->a_type == A_FLOAT) {
			rad = (double)argv->a_w.w_float;	
		}
	


		x->o[id].radius=rad*MAX(x->boxspeakersdim.x, x->boxspeakersdim.y);
//		post("snd %ld diffusion is %f, like setting the radius at %f",id,rad,x->o[id].radius);
		dbap3d_calc_snd(x,id);

		}
		
	}



}

void dbap3d_diffusion(dbap3d *x, t_symbol *msg, short argc, t_atom *argv) //diffusion tris to work on 0.1. range
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
			dbap3d_calc_snd(x,i);
		}
		
	}



}


void dbap3d_gain(dbap3d *x, t_symbol *msg, short argc, t_atom *argv)
{
	double gain;
	int i;
	
	if(argv->a_type == A_FLOAT) {
		gain = (double)argv->a_w.w_float;	
	}
	
	

	for(i=0;i<x->num_sounds;i++) {
	
		x->o[i].gain = gain;
		dbap3d_calc_snd(x,i);
		}
}

void dbap3d_snd_gain(dbap3d *x, t_symbol *msg, short argc, t_atom *argv)
{
	if (argc) {
		double gain;
//		int i;
		int id;

		if (argc==2) { //id rad
			id = (argv+0)->a_w.w_long;
			argv++;

		if(argv->a_type == A_FLOAT) {
			gain = (double)argv->a_w.w_float;	
		}
	


		x->o[id].gain=gain;
		dbap3d_calc_snd(x,id);

		}
		
	}



}


void dbap3d_db(dbap3d *x, t_symbol *msg, short argc, t_atom *argv)
{
	double gain;
	int i;
	
	if(argv->a_type == A_FLOAT) {
		gain = (double)argv->a_w.w_float;	
	}
	
	gain = pow(10.,(gain/20.)); //db to a

	for(i=0;i<x->num_sounds;i++) {
	
		x->o[i].gain = gain;
		dbap3d_calc_snd(x,i);
		}
}

void dbap3d_snd_db(dbap3d *x, t_symbol *msg, short argc, t_atom *argv)
{
	if (argc) {
		double gain;
//		int i;
		int id;

		if (argc==2) { //id rad
			id = (argv+0)->a_w.w_long;
			argv++;

		if(argv->a_type == A_FLOAT) {
			gain = (double)argv->a_w.w_float;	
		}
	
		gain = pow(10.,(gain/20.)); //db to a


		x->o[id].gain=gain;
		dbap3d_calc_snd(x,id);

		}
		
	}



}
void dbap3d_region_active(dbap3d *x, int box, int active)
{
	int s = MAX(0,box);
	x->boxes[s].active = active;
}

void dbap3d_speaker_active(dbap3d *x, int speaker, int active)
{
	int s = MAX(0,speaker);
	x->woofers[s].active = active;
}



void dbap3d_snd_active(dbap3d *x, int snd, int active)
{
	int s = MAX(0,snd);
	x->o[s].active = active;
}


void dbap3d_num_sounds(dbap3d *x, int numsnds)
{
	x->num_sounds = MAX(0,numsnds);
}

void dbap3d_num_regions(dbap3d *x, int numregions)
{
	x->num_regions = MAX(0,numregions);
}

void dbap3d_num_speakers(dbap3d *x, int numregions)
{
	x->num_woofers = MAX(0,numregions);
}


void dbap3d_om(dbap3d *x, t_symbol *msg, short argc, t_atom *argv)
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
      dbap3d_calc_snd(x, i);
    }
  } else {
    
    for(i=0;i<num;i++) {
      if (argv[i].a_type == A_FLOAT) { temp[i] = (int) argv[i].a_w.w_float; }
      
      x->o[i].om = temp[i];
    }
    // recalculate
    for(i=0;i<num;i++) {
      dbap3d_calc_snd(x, i);
    }
  }
}
void dbap3d_falloff_mode(dbap3d *x, int numsnds)
{
	x->fall_off_mode = numsnds;
}

void dbap3d_falloff(dbap3d *x, t_symbol *msg, short argc, t_atom *argv)
{
	int i;
	double x1;

		if(argv->a_type == A_FLOAT) {
			x1 = (double)argv->a_w.w_float;	
		}
	
	
//	x->fall_off = x1*0.00001;

	for(i=0;i<x->num_sounds;i++) {
		x->o[i].f = x1*0.00001;
//		dbap3d_calc_snd(x,i);
		}
}

void dbap3d_snd_falloff(dbap3d *x, t_symbol *msg, short argc, t_atom *argv)
{
	int snd;
	double x1;

		if(argv->a_type == A_FLOAT) {
			snd = (int)argv->a_w.w_float;	post("should be int. truncating to %ld", snd);
		}
		
		argv++;

		if(argv->a_type == A_FLOAT) {
			x1 = (double)argv->a_w.w_float;	
		}


		x->o[snd].f = x1*0.00001;
		dbap3d_calc_snd(x,snd);
}

void dbap3d_setboxscale (dbap3d *x, t_symbol *msg, short argc, t_atom *argv)
{
	double x1;

		if(argv->a_type == A_FLOAT) {
			x1 = (double)argv->a_w.w_float;	
		}  else { post("feed me a float please"); return; };


//	post("setting box_scale to %f",x1);
	x->box_scale = x1;
	
	dbap3d_box_scale(x);

//	for(i=0;i<x->num_sounds;i++) dbap3d_calc_snd(x,i);
		

}


 double dbap3d_box_check (dbap3d *x, vec coords, double radius, double falloff, double distmass) 
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

 double dbap3d_box_check2 (dbap3d *x, vec coords, double radius, double falloff, double distmass) 
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

void dbap3d_box_scale (dbap3d *x)
{
	int i;
	vec min,max,s;
	t_atom *out;
	vec center;
	
	out = x->out;

	 //init locals to first speaker
//	set (min, x->speakers[0], x->speakers[1]);
//	set (max, x->speakers[0], x->speakers[1]);
	set (min, x->woofers[0].pos.x, x->woofers[0].pos.y,x->woofers[0].pos.z);
	set (max, x->woofers[0].pos.x, x->woofers[0].pos.y,x->woofers[0].pos.z);

	for(i=0;i<x->num_woofers;i++) {
//		set (s, x->speakers[i*DIMEN+0],x->speakers[i*DIMEN+1]);
		set (s, x->woofers[i].pos.x, x->woofers[i].pos.y,x->woofers[0].pos.z);
		// min max
		if(min.x > s.x) { min.x = s.x; }	
		if(min.y > s.y) { min.y = s.y; }	
		if(min.z > s.z) { min.z = s.z; }	
		if(max.x < s.x) { max.x = s.x; }	
		if(max.y < s.y) { max.y = s.y; }	
		if(max.z < s.z) { max.z = s.z; }	

	}


	copy (x->boxspeakersmin,min);
	copy (x->boxspeakersmax,max);

	
	// box center
	center.x = (max.x + min.x) * 0.5f; 
	center.y = (max.y + min.y) * 0.5f; 
	center.z = (max.z + min.z) * 0.5f; 

	copy (x->boxcenter,center);


	// box speakers dim
	s.x = (max.x - min.x) ; 
	s.y = (max.y - min.y) ; 
	s.z = (max.z - min.z) ; 

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
	s.z = (max.z - min.z) ; 

	copy (x->boxdim,s);


	// box radius
	s.x = x->boxdim.x*0.5 ; 
	s.y = x->boxdim.y*0.5 ; 
	s.z = x->boxdim.z*0.5 ; 

	copy (x->boxradius,s);
	
		
	SETFLOAT(out+0, x->boxspeakersmin.x);
	SETFLOAT(out+1, x->boxspeakersmin.y);
	SETFLOAT(out+2, x->boxspeakersmin.z);
	SETFLOAT(out+3, x->boxspeakersmax.x);
	SETFLOAT(out+4, x->boxspeakersmax.y);
	SETFLOAT(out+5, x->boxspeakersmax.z);
		
	outlet_anything(x->c_out2, ps_boxspeakers, 2*DIMEN, out);


	SETFLOAT(out+0, x->boxcoordsmin.x);
	SETFLOAT(out+1, x->boxcoordsmin.y);
	SETFLOAT(out+2, x->boxcoordsmin.z);
	SETFLOAT(out+3, x->boxcoordsmax.x);
	SETFLOAT(out+4, x->boxcoordsmax.y);
	SETFLOAT(out+5, x->boxcoordsmax.z);
		
	outlet_anything(x->c_out2, ps_boxcoords, 2*DIMEN, out);


	SETFLOAT(out+0, x->boxcenter.x);
	SETFLOAT(out+1, x->boxcenter.y);
	SETFLOAT(out+2, x->boxcenter.z);

	outlet_anything(x->c_out2, ps_boxcenter, DIMEN, out);

	SETFLOAT(out+0, x->boxdim.x);
	SETFLOAT(out+1, x->boxdim.y);
	SETFLOAT(out+2, x->boxdim.z);

	outlet_anything(x->c_out2, ps_boxdim, DIMEN, out);

	SETFLOAT(out+0, x->boxspeakersdim.x);
	SETFLOAT(out+1, x->boxspeakersdim.y);
	SETFLOAT(out+2, x->boxspeakersdim.z);

	outlet_anything(x->c_out2, ps_boxspeakersdim, DIMEN, out);

	SETFLOAT(out+0, x->boxradius.x);
	SETFLOAT(out+1, x->boxradius.y);
	SETFLOAT(out+2, x->boxradius.z);

	outlet_anything(x->c_out2, ps_boxradius, DIMEN, out);
	

}

void dbap3d_setboxcoords (dbap3d *x, t_symbol *msg, short argc, t_atom *argv)
{
	double temp[DIMEN*2];
	int i;
	t_atom *out;

	out = x->out;
	
	for(i=0;i<DIMEN*2;i++) {
		if (argv[i].a_type == A_FLOAT) { temp[i] = (double) argv[i].a_w.w_float; }
		
	}
	
	x->boxspeakersmin.x = temp[0];
	x->boxspeakersmin.y = temp[1];
	x->boxspeakersmin.z = temp[2];
	x->boxspeakersmax.x = temp[3];
	x->boxspeakersmax.y = temp[4];
	x->boxspeakersmax.z = temp[5];
	SETFLOAT(out+0, x->boxspeakersmin.x);
	SETFLOAT(out+1, x->boxspeakersmin.y);
	SETFLOAT(out+2, x->boxspeakersmin.z);
	SETFLOAT(out+3, x->boxspeakersmax.x);
	SETFLOAT(out+4, x->boxspeakersmax.y);
	SETFLOAT(out+5, x->boxspeakersmax.z);
		
	outlet_anything(x->c_out2, ps_boxspeakers, 2*DIMEN, out);

}

void dbap3d_setboxcenter (dbap3d *x, t_symbol *msg, short argc, t_atom *argv)
{
	double temp[DIMEN];
	int i;
	t_atom *out=x->out;

	for(i=0;i<DIMEN;i++) {
		if (argv[i].a_type == A_FLOAT) { temp[i] = (double) argv[i].a_w.w_float; }
		
	}
	
	x->boxcenter.x = temp[0];
	x->boxcenter.y = temp[1];
	x->boxcenter.z = temp[2];

	SETFLOAT(out+0, x->boxcenter.x);
	SETFLOAT(out+1, x->boxcenter.y);
	SETFLOAT(out+2, x->boxcenter.z);

	outlet_anything(x->c_out2, ps_boxcenter, DIMEN, out);


}

void dbap3d_setboxdim (dbap3d *x, t_symbol *msg, short argc, t_atom *argv)
{
	double temp[DIMEN];
	int i;
	t_atom *out = x->out;

	for(i=0;i<DIMEN;i++) {
		if (argv[i].a_type == A_FLOAT) { temp[i] = (double) argv[i].a_w.w_float; }
		
	}
	
	x->boxdim.x = temp[0];
	x->boxdim.y = temp[1];
	x->boxdim.z = temp[2];

	SETFLOAT(out+0, x->boxdim.x);
	SETFLOAT(out+1, x->boxdim.y);
	SETFLOAT(out+2, x->boxdim.z);

	outlet_anything(x->c_out2, ps_boxdim, DIMEN, out);

}

void dbap3d_setboxmode (dbap3d *x, int x1)
{
	x->box_mode = x1;
}

void dbap3d_setbullseye (dbap3d *x, int x1)
{
	if(x1)
		x->bullseye=1;
	else
		x->bullseye=0;
}

void dbap3d_speakers (dbap3d *x, t_symbol *msg, short argc, t_atom *argv)
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
			x->woofers[i].pos.z = temp[i*DIMEN+2];
		}

		
		outlet_int(x->c_out2, x->num_woofers);		
		

		if(x->box_scale) { //make box coords regarding speakers position
			
			dbap3d_box_scale(x);

		}

	} else {
		error(" no arguments via speakers method");
	}
}


void dbap3d_setspeaker (dbap3d *x, t_symbol *msg, short argc, t_atom *argv) // per speaker pos
{
	if (argc) {
		double	temp[DIMEN]={0.0f,0.0f,0.0f};
		int i;
		int id;

		if (argc==(DIMEN+1)) { //id posx posy posz
			id = (argv+0)->a_w.w_long;
			argv++;

			for(i=0;i<DIMEN;i++) {
				if(argv->a_type == A_FLOAT) {
					temp[i] = (double)argv->a_w.w_float;	
					argv++;
				}
			}


			x->woofers[id].pos.x = temp[+0];
			x->woofers[id].pos.y = temp[+1];
			x->woofers[id].pos.z = temp[+2];

//		x->speakers[id*DIMEN+0] = temp[0];
//		x->speakers[id*DIMEN+1] = temp[1];

		}
		
	}


}


void dbap3d_list (dbap3d *x, t_symbol *msg, short argc, t_atom *argv) //updated for <index> <pos> <optional active
{



	if (argc) {
	
		double	temp[DIMEN]={0.0f,0.0f,0.0f};
		int i;
		int id;
		int status;

		if (argc<(DIMEN+1)) { 
			post("dbap's list method requires: snd index(int) + posx + posy + posz(floats) + (optional active status0/1)");
			return;
		}

		if (argc==(DIMEN+2)) { //id posx posy status
			status = (argv+4)->a_w.w_long;
			x->o[id].active = status;
		}
		
		//here only if
		//if (argc==(DIMEN+1)){ //id posx posy
			
		id = (argv+0)->a_w.w_long;
		argv++;

		for(i=0;i<DIMEN;i++) {
			if(argv->a_type == A_FLOAT) {
				temp[i] = (double)argv->a_w.w_float;	
				argv++;
			}
		}



		set(x->o[id].pos, temp[0], temp[1],temp[2]); 
		dbap3d_calc_snd(x, id);
		
		
	}




}



void dbap3d_setregions (dbap3d *x, int x1)
{
	if(x1)
		x->regions=1;
	else
		x->regions=0;
}

void dbap3d_regions (dbap3d *x, t_symbol *msg, short argc, t_atom *argv)
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
			x->boxes[i].min.z = temp[i*DIMEN*2+2];
			x->boxes[i].max.x = temp[i*DIMEN*2+3];
			x->boxes[i].max.y = temp[i*DIMEN*2+4];
			x->boxes[i].max.z = temp[i*DIMEN*2+5];
		}

		
//		outlet_int(x->c_out2, x->num_woofers);		
		

	} else {
		error(" no arguments via speakers method");
	}
}


void dbap3d_region (dbap3d *x, t_symbol *msg, short argc, t_atom *argv) // per speaker pos
{
	if (argc) {
		double	temp[DIMEN]={0.0f,0.0f,0.0f};
		int i;
		int id;

		if (argc==7) { //id minx miny minz maxx maxy maxz
			id = (argv+0)->a_w.w_long;
			argv++;

			for(i=0;i<DIMEN*2;i++) {
				if(argv->a_type == A_FLOAT) {
					temp[i] = (double)argv->a_w.w_float;	
					argv++;
				}
			}


			x->boxes[id].min.x = temp[+0];
			x->boxes[id].min.y = temp[+1];
			x->boxes[id].min.z = temp[+2];
			x->boxes[id].max.x = temp[+3];
			x->boxes[id].max.y = temp[+4];
			x->boxes[id].max.z = temp[+5];

//		x->speakers[id*DIMEN+0] = temp[0];
//		x->speakers[id*DIMEN+1] = temp[1];

		}
		
	}


}






int dbap3d_check_bullseye(dbap3d *x, int snd)
{
	int i;
	vec speaker,sndpt,disttemp;
	double dist;

	set(sndpt, x->o[snd].pos.x,x->o[snd].pos.y,x->o[snd].pos.z);
	for(i=0; i < x->num_woofers; i++ ) {
		set(speaker, 
			x->woofers[i].pos.x,
			x->woofers[i].pos.y,
			x->woofers[i].pos.z);
		
	// cityblock distance
		disttemp.x = speaker.x - sndpt.x;
		disttemp.y = speaker.y - sndpt.y;
		disttemp.z = speaker.z - sndpt.z;
		if(disttemp.x<0.) disttemp.x *= -1;
		if(disttemp.y<0.) disttemp.y *= -1;
		if(disttemp.z<0.) disttemp.z *= -1;

		dist = disttemp.x + disttemp.y + disttemp.z;
		
		if ( dist < x->o[snd].radius ) //cityblock distance is not eucleadean, but should work for the purposes
			return i;
	}
	
	return -1;
	
}

void dbap3d_exec_bullseye(dbap3d *x, int snd, int speaker)
{
	t_atom	*out;
	int i;
	
	out = x->out;

			// output 01 levels
			for (i=0; i<x->num_woofers;i++) 
			{
				SETLONG(out+0, snd);
				SETLONG(out+1, i);
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

void dbap3d_exec_silence(dbap3d *x, int snd)
{
	t_atom	*out;
	int i;
	
	out = x->out;

			// output 01 levels
			for (i=0; i<x->num_woofers;i++) 
			{
				SETLONG(out+0, snd);
				SETLONG(out+1, i);
				 SETFLOAT(out+2,0.);
				
				outlet_list(x->c_out, 0L, 3, out);

			}

}


/// main method per sound
//static inline dbap3d_calc_snd(dbap3d *x, int snd)
void dbap3d_calc_snd(dbap3d *x, int snd)
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

	set(inpt, x->o[snd].pos.x,x->o[snd].pos.y,x->o[snd].pos.z);
	radius = x->o[snd].radius;
	falloff = x->o[snd].f;
	
	if(x->o[snd].active == 0) {
		dbap3d_exec_silence(x,snd); // if inactive, output zero's &return
		return; 
	}

	switch(x->o[snd].om) 
	{

		case 0:
		// normalize & no box (fixed according to tl.dbap_1_over_r
		
			// pass 0 - bullseye method
		
			if(x->bullseye) {
					bingo=dbap3d_check_bullseye(x,snd);
					if(bingo!=-1) {
						dbap3d_exec_bullseye(x,snd,bingo);
						return;
					}
			}

			radius *= radius; //square radius from now on

			// pass 1 - calculate squares of distances
						
			for (i=0; i<num_woofers;i++) 
			{
			 if(x->woofers[i].active) {
				set(speaker,x->woofers[i].pos.x,x->woofers[i].pos.y,x->woofers[i].pos.z);
				
				avec.x = speaker.x-inpt.x;
				avec.y = speaker.y-inpt.y;
				avec.z = speaker.z-inpt.z;

				avec.x *= avec.x;
				avec.y *= avec.y;
				avec.z *= avec.z;

				temp[i] = 1. / (avec.x + avec.y + avec.z + radius); //adding squared radius/diffusion parameter
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
					SETLONG(out+0, snd);
					SETLONG(out+1, i);

					if(x->o[snd].gain==1.0f) SETFLOAT(out+2, temp[i]);
					else SETFLOAT(out+2, temp[i]*x->o[snd].gain);
					
					outlet_list(x->c_out, 0L, 3, out);
	
				} else if(!x->woofers[i].active)  { //output 0's on the inactive speaker
					SETLONG(out+0, snd);
					SETLONG(out+1, i);
					SETFLOAT(out+2, 0.);
					outlet_list(x->c_out, 0L, 3, out);
				}
			}
			break;


		case 1:

			// pass 0 - bullseye method
		
			if(x->bullseye) {
					bingo=dbap3d_check_bullseye(x,snd);
					if(bingo!=-1) {
						dbap3d_exec_bullseye(x,snd,bingo);
						return;
					}
			}



			// no normalize & no box
			for (i=0; i<num_woofers;i++) 
			{
			 if(x->woofers[i].active) 
			 {
				set(speaker,x->woofers[i].pos.x,x->woofers[i].pos.y,x->woofers[i].pos.z);
				
				distance = distan2(speaker,inpt);

				distance -= radius;

				distance = MAX(radius,distance);

				distance = sqr(distance);
				distance = 1.0f/distance;
				temp[i] = sqr(distance);

				sum += temp[i];
			 }
			}
			
			one_over_sum = 1.00f/sum;
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
					SETLONG(out+0, snd);
					SETLONG(out+1, i);
					if(x->o[snd].gain==1.0f) SETFLOAT(out+2, temp[i]);
					else SETFLOAT(out+2, temp[i]*x->o[snd].gain);
					
					outlet_list(x->c_out, 0L, 3, out);

				} else if(!x->woofers[i].active)  { //output 0's on the inactive speaker
					SETLONG(out+0, snd);
					SETLONG(out+1, i);
					SETFLOAT(out+2, 0.);
					outlet_list(x->c_out, 0L, 3, out);
				}

			}

			break;

		case 2:
			// pass 0 - bullseye method
		
			if(x->bullseye) {
					bingo=dbap3d_check_bullseye(x,snd);
					if(bingo!=-1) {
						dbap3d_exec_bullseye(x,snd,bingo);
						return;
					}
			}


		// normalize &  box
	
			for (i=0; i<num_woofers;i++) 
			{
				if(x->woofers[i].active)
				{
				set(speaker,x->woofers[i].pos.x,x->woofers[i].pos.y,x->woofers[i].pos.z);
				
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
			
			one_over_sum = 1.00f/sum;
			SETFLOAT(out+0, snd);
			SETFLOAT(out+1, one_over_sum);
			outlet_anything(x->c_out2, ps_distmass, 2, out);

//			post("normalize&box inpt %f %f", inpt.x, inpt.y);
			boxfactor = dbap3d_box_check(x, inpt, radius, falloff, one_over_sum);
//			post("snd %ld boxfactor %f", snd, boxfactor);

			
			// output
			for (i=0; i<num_woofers;i++) 
			{
				if(x->woofers[i].active)
				{
					temp[i] = temp[i]*one_over_sum;
					temp[i] = pow(temp[i],0.5);
					temp[i]*= boxfactor;		// scale box factor
					SETLONG(out+0, snd);
					SETLONG(out+1, i);
          
					if(x->o[snd].gain==1.0f) SETFLOAT(out+2, temp[i]); // what for??
					else SETFLOAT(out+2, temp[i]*x->o[snd].gain);
					
					outlet_list(x->c_out, 0L, 3, out);

				} else if(!x->woofers[i].active)  { //output 0's on the inactive speaker
					SETLONG(out+0, snd);
					SETLONG(out+1, i);
					SETFLOAT(out+2, 0.);
					outlet_list(x->c_out, 0L, 3, out);
				}
			}
			break;


		case 3:

			// pass 0 - bullseye method
		
			if(x->bullseye) {
					bingo=dbap3d_check_bullseye(x,snd);
					if(bingo!=-1) {
						dbap3d_exec_bullseye(x,snd,bingo);
						return;
					}
			}


			// no normalize &  box
			for (i=0; i<num_woofers;i++) 
			{
				if(x->woofers[i].active)
				{
				set(speaker,x->woofers[i].pos.x,x->woofers[i].pos.y,x->woofers[i].pos.z);
				
				distance = distan2(speaker,inpt);

				distance -= radius;

				distance = MAX(radius,distance);

				distance = sqr(distance);
				distance = 1.0f/distance;
				temp[i] = sqr(distance);

				sum += temp[i];
				}
			}
			
			one_over_sum = 1.00f/sum;
			SETFLOAT(out+0, snd);
			SETFLOAT(out+1, one_over_sum);
			outlet_anything(x->c_out2, ps_distmass, 2, out);

	//		post("normalize&box inpt %f %f", inpt.x, inpt.y);
			boxfactor = dbap3d_box_check2(x, inpt,radius,falloff, one_over_sum);
	//		post("snd %ld boxfactor %f", snd, boxfactor);
			
			// output
			for (i=0; i<num_woofers;i++) 
			{
				if(x->woofers[i].active)
				{
					temp[i] = temp[i]*one_over_sum;
				//	temp[i] = pow(temp[i],0.5);
					temp[i]*= boxfactor;		// scale box factor
					SETLONG(out+0, snd);
					SETLONG(out+1, i);
					if(x->o[snd].gain==1.0f) SETFLOAT(out+2, temp[i]);
					else SETFLOAT(out+2, temp[i]*x->o[snd].gain);
					
					outlet_list(x->c_out, 0L, 3, out);

				} else if(!x->woofers[i].active)  { //output 0's on the inactive speaker
					SETLONG(out+0, snd);
					SETLONG(out+1, i);
					SETFLOAT(out+2, 0.);
					outlet_list(x->c_out, 0L, 3, out);
				}
			}

			break;




			case 4:
		// normalize & no box (fixed according to tl.dbap_1_over_r -- non sqrt mode
		
			// pass 0 - bullseye method
		
			if(x->bullseye) {
					bingo=dbap3d_check_bullseye(x,snd);
					if(bingo!=-1) {
						dbap3d_exec_bullseye(x,snd,bingo);
						return;
					}
			}

			radius *= radius; //square radius from now on

			// pass 1 - calculate squares of distances
						
			for (i=0; i<num_woofers;i++) 
			{
				if(x->woofers[i].active)
				{
				set(speaker,x->woofers[i].pos.x,x->woofers[i].pos.y,x->woofers[i].pos.z);
				
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
					SETLONG(out+0, snd);
					SETLONG(out+1, i);

					if(x->o[snd].gain==1.0f) SETFLOAT(out+2, temp[i]);
					else SETFLOAT(out+2, temp[i]*x->o[snd].gain);
					
					outlet_list(x->c_out, 0L, 3, out);

				} else if(!x->woofers[i].active)  { //output 0's on the inactive speaker
					SETLONG(out+0, snd);
					SETLONG(out+1, i);
					SETFLOAT(out+2, 0.);
					outlet_list(x->c_out, 0L, 3, out);
				}
			}
			break;

	
	
	
	}


}



void dbap3d_setsndc(dbap3d *x, t_symbol *msg, short argc, t_atom *argv)
{
	
	if (argc) {
		double temp[DIMEN] = {0.,0.,0.};
		int i;
		int id;

		if (argc==(DIMEN+1)) { //id posx posy
			id = (argv+0)->a_w.w_long;
			argv++;

			for(i=0;i<DIMEN;i++) {
				if(argv->a_type == A_FLOAT) {
					temp[i] = (double)argv->a_w.w_float;	
					argv++;
				}
			}


		x->o[id].pos.x=temp[0];
		x->o[id].pos.y=temp[1];
		x->o[id].pos.z=temp[2];
		
		dbap3d_calc_snd(x, id);
		}
		
	}



//	set(x->o[snd].pos, (double)x1, (double)x2);
//	
//	dbap3d_calc_snd(x, snd);

	
}

 double distan(vec o, vec p)
{
	vec dist;
	
	dist.x = p.x - o.x;
	dist.y = p.y - o.y;
	dist.z = p.z - o.z;
	
	return(len(dist));
	
}

 double distan2(vec o, vec p)
{
	vec dist;
	
	dist.x = p.x - o.x;
	dist.y = p.y - o.y;
	dist.z = p.z - o.z;
		
	return(len2(dist));	
}





void dbap3d_bang (dbap3d *x)
{


//	dbap3d_calc_snd(x,0);
}


void dbap3d_getspeakers (dbap3d *x, t_symbol *msg, short argc, t_atom *argv)
{

	t_atom  out[MAXVALUE];
	int i;
	

	for(i=0;i<x->num_woofers;i++) {
		SETFLOAT(out+i*DIMEN+0,x->woofers[i].pos.x);
		SETFLOAT(out+i*DIMEN+1,x->woofers[i].pos.y);
		SETFLOAT(out+i*DIMEN+2,x->woofers[i].pos.z);

	}
	outlet_anything(x->c_out2, ps_getspeakers, x->num_woofers*DIMEN, out);

}


void dbap3d_version (dbap3d *x, t_symbol *msg, short argc, t_atom *argv)
{
	t_atom * out;
	
	out = x->out;
	
	SETSYM(out+0,ps_myversion);
	outlet_anything(x->c_out2, ps_version, 1, out);

}

void dbap3d_assist(dbap3d *x, void *b, long m, long a, char *s)
{
    if (m==1) { post("list of coords, speakers"); }
    else if (m==2&&a==0) { post("(list) to matrix~"); }
    else if (m==2&&a==1) { post("info"); }
}

void setup(void)
{
  //long int tick = gettime();
  dbap3d_class = class_new(gensym("dbap3d"), (t_newmethod)dbap3d_new, (t_method)dbap3d_free, sizeof(dbap3d), 0, A_GIMME, 0);
  
  //setup((t_messlist**)&dbap3d_class,(method)dbap3d_new,0L,(short)sizeof(dbap3d),0L,
  //A_GIMME,0);
  
  // addbang(dbap3d_class, (t_method)dbap3d_bang);
  class_addmethod(dbap3d_class, (t_method)dbap3d_list,	gensym("list"),	A_GIMME,0);
  class_addmethod(dbap3d_class, (t_method)dbap3d_speakers,				gensym("speakers"),			A_GIMME,0);
  class_addmethod(dbap3d_class, (t_method)dbap3d_getspeakers,			gensym("getspeakers"),		A_GIMME,0);
  class_addmethod(dbap3d_class, (t_method)dbap3d_setspeaker,				gensym("speaker"),			A_GIMME,0);
  class_addmethod(dbap3d_class, (t_method)dbap3d_num_sounds,				gensym("num"),				A_DEFLONG,0);
  class_addmethod(dbap3d_class, (t_method)dbap3d_radius,					gensym("radius"),			A_GIMME,0);
  class_addmethod(dbap3d_class, (t_method)dbap3d_snd_radius,				gensym("radius_snd"),			A_GIMME,0);
  class_addmethod(dbap3d_class, (t_method)dbap3d_snd_radius,				gensym("radius_input"),			A_GIMME,0);
  class_addmethod(dbap3d_class, (t_method)dbap3d_diffusion,				gensym("diffusion"),		A_GIMME,0);
  class_addmethod(dbap3d_class, (t_method)dbap3d_snd_diffusion,			gensym("diffusion_snd"),		A_GIMME,0);
  class_addmethod(dbap3d_class, (t_method)dbap3d_snd_diffusion,			gensym("diffusion_input"),		A_GIMME,0);
  class_addmethod(dbap3d_class, (t_method)dbap3d_gain,					gensym("gain"),				A_GIMME,0);
  class_addmethod(dbap3d_class, (t_method)dbap3d_snd_gain,				gensym("gain_snd"),				A_GIMME,0);
  class_addmethod(dbap3d_class, (t_method)dbap3d_snd_gain,				gensym("gain_input"),			A_GIMME,0);
  class_addmethod(dbap3d_class, (t_method)dbap3d_db,						gensym("db"),				A_GIMME,0);
  class_addmethod(dbap3d_class, (t_method)dbap3d_snd_db,					gensym("db_snd"),				A_GIMME,0);
  class_addmethod(dbap3d_class, (t_method)dbap3d_snd_db,					gensym("db_input"),				A_GIMME,0);
  class_addmethod(dbap3d_class, (t_method)dbap3d_setsndc,				gensym("set"),				A_GIMME,0);
  class_addmethod(dbap3d_class, (t_method)dbap3d_om,						gensym("om"),				A_GIMME,0);
  class_addmethod(dbap3d_class, (t_method)dbap3d_falloff_mode,			gensym("falloff_mode"),		A_LONG,0);
  class_addmethod(dbap3d_class, (t_method)dbap3d_falloff,				gensym("falloff"),			A_GIMME,0);
  class_addmethod(dbap3d_class, (t_method)dbap3d_snd_falloff,			gensym("falloff_snd"),			A_GIMME,0);
  class_addmethod(dbap3d_class, (t_method)dbap3d_snd_falloff,			gensym("falloff_input"),		A_GIMME,0);
  
  class_addmethod(dbap3d_class, (t_method)dbap3d_setboxcoords,			gensym("boxcoords"),		A_GIMME,0);
  class_addmethod(dbap3d_class, (t_method)dbap3d_setboxcenter,			gensym("boxcenter"),		A_GIMME,0);
  class_addmethod(dbap3d_class, (t_method)dbap3d_setboxdim,				gensym("boxdim"),			A_GIMME,0);
  class_addmethod(dbap3d_class, (t_method)dbap3d_setboxscale,			gensym("boxscale"),			A_GIMME,0);
  class_addmethod(dbap3d_class, (t_method)dbap3d_setboxmode,				gensym("boxmode"),			A_LONG,0);
  class_addmethod(dbap3d_class, (t_method)dbap3d_setbullseye,			gensym("bullseye"),			A_LONG,0);
  
  class_addmethod(dbap3d_class, (t_method)dbap3d_snd_active,				gensym("active"),			A_DEFLONG,A_DEFLONG,0);
  class_addmethod(dbap3d_class, (t_method)dbap3d_snd_active,				gensym("active_snd"),			A_DEFLONG,A_DEFLONG,0);
  class_addmethod(dbap3d_class, (t_method)dbap3d_snd_active,				gensym("active_input"),			A_DEFLONG,A_DEFLONG,0);
  class_addmethod(dbap3d_class, (t_method)dbap3d_speaker_active,			gensym("gactive_speaker"),			A_DEFLONG,A_DEFLONG,0);
  class_addmethod(dbap3d_class, (t_method)dbap3d_region_active,				gensym("gactive_regions"),			A_DEFLONG,A_DEFLONG,0);
  
  class_addmethod(dbap3d_class, (t_method)dbap3d_setregions,				gensym("gregions"),			A_LONG,0);
  class_addmethod(dbap3d_class, (t_method)dbap3d_regions,				gensym("gregions_all"),			A_GIMME,0);
  class_addmethod(dbap3d_class, (t_method)dbap3d_region,					gensym("gregion"),			A_GIMME,0);
  
  class_addmethod(dbap3d_class, (t_method)dbap3d_num_sounds,				gensym("gnum_sounds"),				A_DEFLONG,0);
  class_addmethod(dbap3d_class, (t_method)dbap3d_num_regions,			gensym("gnum_regions"),				A_DEFLONG,0);
  class_addmethod(dbap3d_class, (t_method)dbap3d_num_speakers,			gensym("gnum_speakers"),				A_DEFLONG,0);
  
  
	
  class_addmethod(dbap3d_class, (t_method)dbap3d_assist,					gensym("gassistgensym"), 			A_CANT, 0);
  
  class_addmethod(dbap3d_class, (t_method)dbap3d_version,				gensym("gversiongensym"),				A_GIMME,0);
  
  post("a-dbap3d by andre sier, inspired by trond lossius");
  
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
