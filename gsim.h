#include "vector.h"
// inserting into a list
#include <list>
using namespace std;

#define GM 6.67384e-11

class System;

class Body
{
public:
	double m_mass;
	Vector3 m_position;
	Vector3 m_velocity;
private:
	friend class System;
	Vector3 m_rk_position;
	Vector3 m_rk_acceleration;
	Vector3 m_kr1;
	Vector3 m_kr2;
	Vector3 m_kr3;
	Vector3 m_kr4;
	Vector3 m_kv1;
	Vector3 m_kv2;
	Vector3 m_kv3;
	Vector3 m_kv4;
public:
	Vector3 Acceleration(Body* p_body);
	Vector3 rkAcceleration(Body* p_body);
};

class System
{
private:
	list<Body*> m_bodies;
public:
	void AddBody(Body* p_body);
	void rkIntegrate( double p_dt );
private:
	void rkAccelerations( void );
	void rkPhase1Positions( void );
	void rkPhase2Positions( double p_dt );
	void rkPhase3Positions( double p_dt );
	void rkPhase4Positions( double p_dt );
	void rkPhase1Integrate( double p_dt );
	void rkPhase2Integrate( double p_dt );
	void rkPhase3Integrate( double p_dt );
	void rkPhase4Integrate( double p_dt );
};
