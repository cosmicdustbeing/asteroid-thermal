// This header file contains procedures that read in shape files for use in a TPM
//
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
//#include <math>

#define PI   3.1415926

//FILE *f_shape;

using namespace std;

void cross_product(vector<float> a, vector<float> b, vector<float> &c){

  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = -a[0] * b[2] + a[2] * b[0];
  c[2] = a[0] * b[1] - a[1] * b[0];
}

void outer_product(vector<float> a, vector<float> b, vector< vector<float> > &c){

  //  vector< vector<float> > c(a.size(),vector<float>(b.size()));

  for (int u=0;u<a.size();u++){
    for (int v=0;v<b.size();v++){
      c[u][v] = a[u]*b[v];
    }
  }
}

void vec_add(vector<float> a, vector<float> b, vector<float> &c){

  c[0] = a[0] + b[0];
  c[1] = a[1] + b[1];
  c[2] = a[2] + b[2];
}

void vec_sub(vector<float> a, vector<float> b, vector<float> &c){

  c[0] = a[0] - b[0];
  c[1] = a[1] - b[1];
  c[2] = a[2] - b[2];
}

void vec_scale(vector<float> a, const float fac, vector<float> &b){

  b[0] = a[0]*fac;
  b[1] = a[1]*fac;
  b[2] = a[2]*fac;
}

float dot_product(vector<float> &a, vector<float> &b){

  float adotb = (a[0]*b[0]) + (a[1]*b[1]) + (a[2]*b[2]);

  return adotb;
}

bool raytri(vector<float> rayOrigin, vector<float> rayVector, 
	    vector<float> vertex0, vector<float> vertex1, vector<float> vertex2){

  const float EPSILON = 1e-10; // small number
  vector<float> edge1(3),edge2(3),h(3),s(3),q(3);
  float a,f,u,v;

  vec_sub(vertex1,vertex0,edge1);
  vec_sub(vertex2,vertex0,edge2);

  cross_product(rayVector,edge2,h);
  a = dot_product(edge1,h);

  if (a > -EPSILON && a < EPSILON) return false; // ray is parallel to this triangle

  f = 1.0/a;
  vec_sub(rayOrigin,vertex0,s);
  u = f*(dot_product(s,h));

  if (u < 0.0 || u > 1.0) return false;

  cross_product(s,edge1,q);
  v = f*(dot_product(rayVector,q));

  if (v < 0.0 || u + v > 1.0) return false;

  float t = f*dot_product(edge2,q);

  if (t > EPSILON && t < 1/EPSILON) return true;
  else return false;
}

void rayplane(vector<float> jedge, vector <float> jvertex, vector<float> inorm,
	      vector<float> avgivert, vector<float> &point){

  vector<float> conn_vec(3);

  vec_sub(jvertex,avgivert,conn_vec);
  double prod1 = dot_product(conn_vec,inorm);
  double prod2 = dot_product(jedge,inorm);
  double prod3 = prod1/prod2;

  vec_scale(jedge,prod3,jedge);
  vec_sub(jvertex,jedge,point);

}
