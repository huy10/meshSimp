#pragma once
#include "stdafx.h"
#include "SimpleObject.h"
#include "Vec3f.h"
#include <boost/numeric/ublas/matrix.hpp> 
#include <boost/numeric/ublas/vector.hpp> 

#define Debug
#define MAX_POINT_DIS 200
#define MAX_ITR 10
//using namespace std;
//using namespace SimpleOBJ;

class _plane{
public: _plane(){}
		~_plane(void){}
private:
	float a,b,c,d;

};

class  _pair
{
public :
	_pair(){}
	~_pair(){}

public :
	int v1;
	int v2;
	float cost;
	//Vec3f vbar;
	std::vector <float>  vbar ;
};
std::vector <_plane> planes;
char * fileName = ".\\test_data\\fixed.perfect.dragon.100K.0.07.obj";