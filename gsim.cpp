#include <math.h>
#include <iostream>

#include "gsim.h"

#define N_ITERATIONS 4096
#define N_BODYS_MAX  100

Vector3 Acceleration12(
	Vector3 b1,
	Vector3 b2,
	double m1,
	double m2)
{
	Vector3 l_r21 = b2 - b1;
	double l_r = l_r21.normalize();
	l_r*=l_r;
	double l_a = GM*m2/l_r;
	return l_a*l_r21;
}

Vector3 Body::Acceleration(Body* p_body)
{
	return Acceleration12(m_position,p_body->m_position,
												m_mass,p_body->m_mass);
}

Vector3 Body::rkAcceleration(Body* p_body)
{
	return Acceleration12(m_rk_position,p_body->m_rk_position,
												m_mass,p_body->m_mass);
}

void System::AddBody(Body* p_body)
{
	m_bodies.push_back(p_body);
}

void System::rkAccelerations( void )
{
	list<Body*>::iterator l_it_b1;
	list<Body*>::iterator l_it_b2;
	// set all accelerations to zero
	for( l_it_b1=m_bodies.begin(); l_it_b1!=m_bodies.end(); l_it_b1++){
		(*l_it_b1)->m_rk_acceleration = Vector3(0,0,0);
	}
	// calculate accelerations for all pairs of bodies
	for( l_it_b1=m_bodies.begin(); l_it_b1!=m_bodies.end(); l_it_b1++){
		l_it_b2 = l_it_b1;
		l_it_b2++;
		for( ; l_it_b2!=m_bodies.end(); l_it_b2++){
			// compute the normal vector of body 1 relative to 2
			Vector3 l_n12 = (*l_it_b1)->m_rk_position - (*l_it_b2)->m_rk_position;
			double l_magr = l_n12.normalize();
			l_magr*=l_magr;
			// normal vector of body 2 relative to 1
			Vector3 l_n21 = -l_n12;
			double l_g=GM/l_magr;
			double l_a1 = (*l_it_b2)->m_mass * l_g;
			double l_a2 = (*l_it_b1)->m_mass * l_g;
			
			(*l_it_b1)->m_rk_acceleration += l_a1*l_n21;
			(*l_it_b2)->m_rk_acceleration += l_a2*l_n12;
		}
	}
}

void System::rkPhase1Positions( void )
{
	list<Body*>::iterator l_it;
	for(l_it=m_bodies.begin(); l_it!=m_bodies.end(); l_it++){
		(*l_it)->m_rk_position = (*l_it)->m_position;
	}
}

void System::rkPhase2Positions( double p_dt )
{
	list<Body*>::iterator l_it;
	for(l_it=m_bodies.begin(); l_it!=m_bodies.end(); l_it++){
		(*l_it)->m_rk_position = (*l_it)->m_position + (*l_it)->m_kr1*p_dt*0.5;
	}
}

void System::rkPhase3Positions( double p_dt )
{
	list<Body*>::iterator l_it;
	for(l_it=m_bodies.begin(); l_it!=m_bodies.end(); l_it++){
		(*l_it)->m_rk_position = (*l_it)->m_position + (*l_it)->m_kr2*p_dt*0.5;
	}
}

void System::rkPhase4Positions( double p_dt )
{
	list<Body*>::iterator l_it;
	for(l_it=m_bodies.begin(); l_it!=m_bodies.end(); l_it++){
		(*l_it)->m_rk_position = (*l_it)->m_position + (*l_it)->m_kr3*p_dt;
	}
}

void System::rkPhase1Integrate( double p_dt )
{
	list<Body*>::iterator l_it;
	for(l_it=m_bodies.begin(); l_it!=m_bodies.end(); l_it++){
		(*l_it)->m_kr1 = (*l_it)->m_velocity;
		(*l_it)->m_kv1 = (*l_it)->m_rk_acceleration;
	}
}

void System::rkPhase2Integrate( double p_dt )
{
	list<Body*>::iterator l_it;
	for(l_it=m_bodies.begin(); l_it!=m_bodies.end(); l_it++){
		(*l_it)->m_kr2 = (*l_it)->m_velocity + (*l_it)->m_kv1*p_dt*0.5;
		(*l_it)->m_kv2 = (*l_it)->m_rk_acceleration;
	}
}

void System::rkPhase3Integrate( double p_dt )
{
	list<Body*>::iterator l_it;
	for(l_it=m_bodies.begin(); l_it!=m_bodies.end(); l_it++){
		(*l_it)->m_kr3 = (*l_it)->m_velocity + (*l_it)->m_kv2*p_dt*0.5;
		(*l_it)->m_kv3 = (*l_it)->m_rk_acceleration;
	}
}

void System::rkPhase4Integrate( double p_dt )
{
	list<Body*>::iterator l_it;
	for(l_it=m_bodies.begin(); l_it!=m_bodies.end(); l_it++){
		(*l_it)->m_kr4 = (*l_it)->m_velocity + (*l_it)->m_kv3*p_dt;
		(*l_it)->m_kv4 = (*l_it)->m_rk_acceleration;
		(*l_it)->m_velocity += p_dt/6.0*((*l_it)->m_kv1 + 2.0* (*l_it)->m_kv2 + 2.0* (*l_it)->m_kv3 + (*l_it)->m_kv4);
		(*l_it)->m_position += p_dt/6.0*((*l_it)->m_kr1 + 2.0* (*l_it)->m_kr2 + 2.0* (*l_it)->m_kr3 + (*l_it)->m_kr4);
	}
}

void System::rkIntegrate( double p_dt )
{
	rkPhase1Positions();
	rkAccelerations();
	rkPhase1Integrate(p_dt);

	rkPhase2Positions(p_dt);
	rkAccelerations();
	rkPhase2Integrate(p_dt);

	rkPhase3Positions(p_dt);
	rkAccelerations();
	rkPhase3Integrate(p_dt);
	
	rkPhase4Positions(p_dt);
	rkAccelerations();
	rkPhase4Integrate(p_dt);
}

void SimulateNBodies( double p_radius, double p_mass, unsigned long N_bodies );


int main(int p_narg, char** p_argv)
{
	Body l_b1;
	Body l_b2;

	double l_magr12 = 384400e3;

	l_b1.m_mass = 5.97219e24;
	l_b2.m_mass = 7.3477e22;

	l_b1.m_position.x = -l_b2.m_mass*l_magr12/(l_b1.m_mass+l_b2.m_mass);
	l_b1.m_position.y = 0.0;
	l_b1.m_position.z = 0.0;

	l_b2.m_position.x = l_b1.m_mass*l_magr12/(l_b1.m_mass+l_b2.m_mass);
	l_b2.m_position.y = 0.0;
	l_b2.m_position.z = 0.0;
	
	double l_omega = sqrt(GM*(l_b1.m_mass+l_b2.m_mass)/pow(l_magr12,3.0));
	
	double l_v1 = l_b2.m_mass/(l_b1.m_mass+l_b2.m_mass)*l_magr12*l_omega;
	double l_v2 = l_b1.m_mass/(l_b1.m_mass+l_b2.m_mass)*l_magr12*l_omega;

	double l_f = l_omega/2.0/M_PI;
	double T = 1.0/l_f;
	double l_dt = T/N_ITERATIONS;
	
	l_b1.m_velocity.x = 0.0;
	l_b1.m_velocity.y = l_v1;
	l_b1.m_velocity.z = 0.0;
	
	l_b2.m_velocity.x = 0.0;
	l_b2.m_velocity.y = -l_v2;
	l_b2.m_velocity.z = 0.0;
	
	System l_system;
	l_system.AddBody(&l_b1);
	l_system.AddBody(&l_b2);

	std::cout.precision(15);
	
	l_b1.m_position.print("b1 start");
	l_b2.m_position.print("b2 start");
	
	for(long n=N_ITERATIONS;n;--n){
		l_system.rkIntegrate( l_dt );
	}
	l_b1.m_position.print("b1 finish");
	l_b2.m_position.print("b2 finish");
	Vector3 l_r12 = l_b1.m_position - l_b2.m_position;
	l_magr12 = l_r12.mag();
	std::cout << "r12:" << l_magr12 << "\n";
	
    //SimulateNBodies( 1.0, 1.0, 4 );
}

void SimulateNBodies( double p_radius, double p_mass, unsigned long N_bodies )
{
	if(N_bodies<2)
		return;
	if(N_bodies>N_BODYS_MAX)
		return;
	
	Body *l_bodies = new Body[N_bodies];
	long n;
	for(n=0;n<N_bodies;n++){
		double theta = (double)n*2*M_PI/N_bodies;
		l_bodies[n].m_position.x = p_radius*cos(theta);
		l_bodies[n].m_position.y = p_radius*sin(theta);
		l_bodies[n].m_position.z = 0.0;
		l_bodies[n].m_mass = p_mass;
	}
	l_bodies[0].m_position.print("body 1 position start");
	Vector3 l_a1(0.0,0.0,0.0);
	for(n=1;n<N_bodies;n++){
		Vector3 l_nk1 = l_bodies[n].m_position - l_bodies[0].m_position;
		double l_magrk1 = l_nk1.normalize();
		l_magrk1*=l_magrk1;
		double l_ak1 = GM*l_bodies[n].m_mass/l_magrk1;
		l_a1 += l_ak1 * l_nk1;
	}
		
	Vector3 l_omega( 0.0, 0.0, sqrt(l_a1.mag()/l_bodies[0].m_position.mag()) );
	for(n=0;n<N_bodies;n++){
		l_bodies[n].m_velocity = l_omega^l_bodies[n].m_position;
	}
	
	System l_system;
	for(n=0;n<N_bodies;n++){
		l_system.AddBody(&l_bodies[n]);
	}

	double l_f = l_omega.z/2.0/M_PI;
	double l_T = 1.0/l_f;
	double l_dT = l_T/N_ITERATIONS;
	
	for(n=N_ITERATIONS;n!=0;n--){
		l_system.rkIntegrate( l_dT );
	}
	
	l_bodies[0].m_position.print("body 1 position finish");

}
