////////////////////////////////////////////////////////////////////////////////////
// Copyright © Charalambos "Charis" Poullis, charalambos@poullis.org    	  //
// This work can only be used under an exclusive license of the author.           //
//////////////////////////////////////////////////////////////////////////////////// 
#include "ONB.h"

#define ONB_EPSILON 0.01f 

void ONB::initFromU( const Vector3f& u ) {
   Vector3f n(1.0f, 0.0f, 0.0f);
   Vector3f m(0.0f, 1.0f, 0.0f);
   
   U = u;
   U.normalize();
   V = U.cross(n);
   if (V.norm() < ONB_EPSILON)
      V = U.cross(m);
   W = U.cross(V);
}

void  ONB::initFromV( const Vector3f& v ) {
   Vector3f n(1.0f, 0.0f, 0.0f);
   Vector3f m(0.0f, 1.0f, 0.0f);
   
   V = v;
   V.normalize();
   U = V.cross(n);
   if (U.norm()*U.norm() < ONB_EPSILON)
      U = V.cross(m);
   W = U.cross(V);
}

void  ONB::initFromW( const Vector3f& w ) {
   Vector3f n(1.0f, 0.0f, 0.0f);
   Vector3f m(0.0f, 1.0f, 0.0f);
   
   W = w;
   W.normalize();
   U = W.cross(n);
   if (U.norm() < ONB_EPSILON)
      U = W.cross(m);
   V = W.cross(U);
}

void ONB::initFromUV( const Vector3f& u, const Vector3f& v ) {
   U = u;
   U.normalize();
   W = u.cross(v);
   W.normalize();
   V = W.cross(U);
}

void ONB::initFromVU( const Vector3f& v, const Vector3f& u ) {
   V = v;
   V.normalize();
   W = u.cross(v);
   W.normalize();
   U = V.cross(W);
}

void ONB::initFromUW( const Vector3f& u, const Vector3f& w ) {
   U = u;
   U.normalize();
   V = w.cross(u);
   V.normalize();
   W = U.cross(V);
}

void  ONB::initFromWU( const Vector3f& w, const Vector3f& u ) {
   W = w;
   W.normalize();
   V = w.cross(u);
   V.normalize();
   U = V.cross(W);
}

void  ONB::initFromVW( const Vector3f& v, const Vector3f& w ) {
   V = v;
   V.normalize();
   U = v.cross(w);
   U.normalize();
   W = U.cross(V);
}

void ONB::initFromWV( const Vector3f& w, const Vector3f& v ) {
   W = w;
   W.normalize();
   U = v.cross(w);
   U.normalize();
   V = W.cross(U);
}

bool  operator==( const ONB & o1, const ONB & o2 )
{ return( o1.u() == o2.u() && o1.v() == o2.v() && o1.w() == o2.w()); }

/*
std::istream & operator>>( std::istream & is, ONB & t ) {
   Vector3f new_u, new_v, new_w;
   is >> new_u >> new_v >> new_w;
   t.initFromUV( new_u, new_v );

   return is;
}
*/

std::ostream & operator<<( std::ostream & os, const ONB & t ) {
   os << t.u() << "\n" << t.v() << "\n" << t.w() << "\n";
   return os;
}


