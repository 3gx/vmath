#pragma once

#include <cstdlib>
#include <vector>
#include "vector3.h"

template<typename REAL>
struct ParticleRef
{
  typedef REAL real;
  typedef Vector3ref<real> vec3ref;
  vec3ref pos;
  vec3ref vel;
  real &mass;
  ParticleRef(
      real &posx, real &posy, real &posz,
      real &velx, real &vely, real &velz,
      real &_mass) : 
    pos(posx, posy, posz), 
    vel(velx, vely, velz), 
    mass(_mass) {};
  ParticleRef(
      const real &posx, const real &posy, const real &posz,
      const real &velx, const real &vely, const real &velz,
      const real &_mass) : 
    pos(posx, posy, posz), 
    vel(velx, vely, velz), 
    mass(_mass) {};
};

template<typename REAL>
struct Particle
{
  typedef REAL real;
  typedef Vector3<real> vec3;

  vec3 pos, vel;
  real mass;
  real iPad;

  Particle() { static_assert(sizeof(Particle() == 8*sizeof(REAL)), "sizeof(Particle) is not correct!"); }
  Particle(const vec3 &_pos, const vec3 &_vel, const real _mass) : 
    pos(_pos), vel(_vel), mass(_mass) { Particle(); }
  Particle(const ParticleRef<REAL> p) : 
    pos(p.pos), vel(p.vel), mass(p.mass) { Particle(); }
};

template<typename REAL>
struct ParticleData
{
  typedef REAL real;
  typedef Vector3ref<real> vec3ref;

  private:
  std::vector<real> posx, posy, posz;
  std::vector<real> velx, vely, velz;
  std::vector<real> mass;
  public:

  size_t size() const { return posx.size(); }
  void resize(const size_t size) 
  {
    posx.resize(size); posy.resize(size); posz.resize(size);
    velx.resize(size); vely.resize(size); velz.resize(size);
    mass.resize(size); 
  }
  void reserve(const size_t size) 
  {
    posx.reserve(size); posy.reserve(size); posz.reserve(size);
    velx.reserve(size); vely.reserve(size); velz.reserve(size);
    mass.reserve(size); 
  }
  ParticleRef<real> operator[] (const int i)
  {
    return ParticleRef<real>(
        posx[i], posy[i], posz[i],
        velx[i], vely[i], velz[i],
        mass[i]);
  }
  const ParticleRef<real> operator[] (const int i) const
  {
    return ParticleRef<real>(
        posx[i], posy[i], posz[i],
        velx[i], vely[i], velz[i],
        mass[i]);
  }
};


