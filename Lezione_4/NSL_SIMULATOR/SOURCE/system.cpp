/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <cmath>
#include <cstdlib>
#include <string>
#include "system.h"

using namespace std;
using namespace arma;

void System :: step(){ // Perform a simulation step
  if(_sim_type == 0) this->Verlet();  // Perform a MD step
  else for(int i=0; i<_npart; i++) this->move(int(_rnd.Rannyu()*_npart)); // Perform a MC step on a randomly choosen particle
  _nattempts += _npart; //update number of attempts performed on the system
  return;
}

void System :: Verlet(){
  double xnew, ynew, znew;
  for(int i=0; i<_npart; i++){ //Force acting on particle i
    _fx(i) = this->Force(i,0);
    _fy(i) = this->Force(i,1);
    _fz(i) = this->Force(i,2);
  }
  for(int i=0; i<_npart; i++){ //Verlet integration scheme
    xnew = this->pbc( 2.0 * _particle(i).getposition(0,true) - _particle(i).getposition(0,false) + _fx(i) * pow(_delta,2), 0);
    ynew = this->pbc( 2.0 * _particle(i).getposition(1,true) - _particle(i).getposition(1,false) + _fy(i) * pow(_delta,2), 1);
    znew = this->pbc( 2.0 * _particle(i).getposition(2,true) - _particle(i).getposition(2,false) + _fz(i) * pow(_delta,2), 2);
    _particle(i).setvelocity(0, this->pbc(xnew - _particle(i).getposition(0,false), 0)/(2.0 * _delta));
    _particle(i).setvelocity(1, this->pbc(ynew - _particle(i).getposition(1,false), 1)/(2.0 * _delta));
    _particle(i).setvelocity(2, this->pbc(znew - _particle(i).getposition(2,false), 2)/(2.0 * _delta));
    _particle(i).acceptmove(); // xold = xnew
    _particle(i).setposition(0, xnew);
    _particle(i).setposition(1, ynew);
    _particle(i).setposition(2, znew);
  }
  _naccepted += _npart;
  return;
}

double System :: Force(int i, int dim){
  double f=0.0, dr;
  vec distance;
  distance.resize(_ndim);
  for (int j=0; j<_npart; j++){
    if(i != j){
      distance(0) = this->pbc( _particle(i).getposition(0,true) - _particle(j).getposition(0,true), 0);
      distance(1) = this->pbc( _particle(i).getposition(1,true) - _particle(j).getposition(1,true), 1);
      distance(2) = this->pbc( _particle(i).getposition(2,true) - _particle(j).getposition(2,true), 2);
      dr = sqrt( dot(distance,distance) );
      if(dr < _r_cut){
        f += distance(dim) * (48.0/pow(dr,14) - 24.0/pow(dr,8));
      }
    }
  }
  return f;
}

void System :: move(int i){ // Propose a MC move for particle i
  if(_sim_type == 3){ //Gibbs sampler for Ising
    
    // random spin
    //int spin = arma::sign(_rnd.Rannyu()-0.5); //ok but not the best 6.732918e-03
    //int spin =1; //unexpectedly better 6.121719e-03
    //int spin = _particle(i).getspin(); //terrible, precise but not accurate
    int spin = std::pow(-1,i); //the best 5.501383e-03
    double delta_E = 2.0 * spin * 
                     ( _J * (_particle(this->pbc(i-1)).getspin() + _particle(this->pbc(i+1)).getspin() ) + _H );
    double acceptance = 1./(1. + exp(-_beta*delta_E));
    if(_rnd.Rannyu() < acceptance){
      _particle(i).setspin(spin);
    }else{
      _particle(i).setspin(-spin);
    }
      _naccepted++;
  } else {           // M(RT)^2
    if(_sim_type == 1){       // LJ system
      vec shift(_ndim);       // Will store the proposed translation
      for(int j=0; j<_ndim; j++){
        shift(j) = _rnd.Rannyu(-1.0,1.0) * _delta; // uniform distribution in [-_delta;_delta)
      }
      _particle(i).translate(shift, _side);  //Call the function Particle::translate
      if(this->metro(i)){ //Metropolis acceptance evaluation
        _particle(i).acceptmove();
        _naccepted++;
      } else _particle(i).moveback(); //If translation is rejected, restore the old configuration
    } else {                  // Ising 1D sim_type = 2
      if(this->metro(i)){     //Metropolis acceptance evaluation for a spin flip involving spin i
        _particle(i).flip();  //If accepted, the spin i is flipped
        _naccepted++;
      }
    }
  }
  return;
}

bool System :: metro(int i){ // Metropolis algorithm
  bool decision = false;
  double delta_E, acceptance;
  if(_sim_type == 1) delta_E = this->Boltzmann(i,true) - this->Boltzmann(i,false);
  else delta_E = 2.0 * _particle(i).getspin() * 
                 ( _J * (_particle(this->pbc(i-1)).getspin() + _particle(this->pbc(i+1)).getspin() ) + _H );
  acceptance = exp(-_beta*delta_E);
  if(_rnd.Rannyu() < acceptance ) decision = true; //Metropolis acceptance step
  return decision;
}

double System :: Boltzmann(int i, bool xnew){
  double energy_i=0.0;
  double dx, dy, dz, dr;
  for (int j=0; j<_npart; j++){
    if(j != i){
      dx = this->pbc(_particle(i).getposition(0,xnew) - _particle(j).getposition(0,1), 0);
      dy = this->pbc(_particle(i).getposition(1,xnew) - _particle(j).getposition(1,1), 1);
      dz = this->pbc(_particle(i).getposition(2,xnew) - _particle(j).getposition(2,1), 2);
      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);
      if(dr < _r_cut){
        energy_i += 1.0/pow(dr,12) - 1.0/pow(dr,6);
      }
    }
  }
  return 4.0 * energy_i;
}

double System :: pbc(double position, int i){ // Enforce periodic boundary conditions
  return position - _side(i) * rint(position / _side(i));
}

int System :: pbc(int i){ // Enforce periodic boundary conditions for spins
  if(i >= _npart) i = i - _npart;
  else if(i < 0)  i = i + _npart;
  return i;
} 

void System :: initialize(int phase){ // Initialize the System object according to the content of the input files in the ../INPUT/ directory

  int p1, p2; // Read from ../INPUT/Primes a pair of numbers to be used to initialize the RNG
  ifstream Primes("../INPUT/Primes");
  Primes >> p1 >> p2 ;
  Primes.close();
  int seed[4]; // Read the seed of the RNG
  ifstream Seed("../INPUT/seed.in");
  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  _rnd.SetRandom(seed,p1,p2);

  ofstream couta("../OUTPUT/acceptance.dat"); // Set the heading line in file ../OUTPUT/acceptance.dat
  couta << "#   N_BLOCK:  ACCEPTANCE:" << '\n';
  couta.close();

  ifstream input{};
  ofstream coutf;
  coutf.open("../OUTPUT/output.dat");

  switch(phase){ // Initialize the system according to the simulation type
    case 0:
      input.open("../INPUT/input.gas");
      coutf << "LJ GAS SIMULATION\n";
      break;
    case 1:
      input.open("../INPUT/input.liquid");
      coutf << "LJ LIQUID SIMULATION\n";
      break;
    case 2:
      input.open("../INPUT/input.solid");
      coutf << "LJ SOLID SIMULATION\n";
      break;
    case 3: 
      input.open("../INPUT/input.ising");
      coutf << "ISING 1D SIMULATION\n";
      break;
    default:
      cerr << "PROBLEM: unknown phase\n";
      exit(EXIT_FAILURE);
  }

  string property;
  double delta;
  while ( !input.eof() ){
    input >> property;
    if( property == "SIMULATION_TYPE" ){
      input >> _sim_type;

      switch (_sim_type) { // Print the simulation type
        case 0:
          coutf << "LJ MOLECULAR DYNAMICS (NVE) SIMULATION\n"; //N = number of particles, V = volume, E = energy
          break;
        case 1:
          coutf << "LJ MONTE CARLO (NVT) SIMULATION\n"; //N = number of particles, V = volume, T = temperature
          break;
        case 2:
          coutf << "ISING 1D MONTE CARLO (M(RT)^2) SIMULATION\n"; // M(RT)^2 = Metropolis, Rosenbluth, Teller 
          break;
        case 3:
          coutf << "ISING 1D MONTE CARLO (GIBBS) SIMULATION\n";
          break;
        default:
          cerr << "PROBLEM: unknown simulation type\n"; 
          exit(EXIT_FAILURE);
        }

      if(_sim_type > 1){
        input >> _J;
        input >> _H;
      }
      
    } else if( property == "RESTART" ){
      input >> _restart;
    } else if( property == "TEMP" ){
      input >> _temp;
      _beta = 1.0/_temp;
      coutf << "TEMPERATURE= " << _temp << '\n';
    } else if( property == "NPART" ){
      input >> _npart;
      _fx.resize(_npart);
      _fy.resize(_npart);
      _fz.resize(_npart);
      _particle.set_size(_npart);
      for(int i=0; i<_npart; i++){ 
        _particle(i).initialize();
        if(_rnd.Rannyu() > 0.5) _particle(i).flip(); // to randomize the spin configuration (high temperature)
      }
      coutf << "NPART= " << _npart << '\n';
    } else if( property == "RHO" ){
      input >> _rho;
      _volume = _npart/_rho;
      _side.resize(_ndim);
      _halfside.resize(_ndim);
      double side = pow(_volume, 1.0/3.0);
      for(int i=0; i<_ndim; i++) _side(i) = side;
      _halfside=0.5*_side;
      coutf << "SIDE= ";
      for(int i=0; i<_ndim; i++){
        coutf << setw(12) << _side[i];
      }
      coutf << '\n';
    } else if( property == "R_CUT" ){
      input >> _r_cut;
      coutf << "R_CUT= " << _r_cut << '\n';
    } else if( property == "DELTA" ){
      input >> delta;
      coutf << "DELTA= " << delta << '\n';
      _delta = delta;
    } else if( property == "NBLOCKS" ){
      input >> _nblocks;
      coutf << "NBLOCKS= " << _nblocks << '\n';
    } else if( property == "NSTEPS" ){
      input >> _nsteps;
      coutf << "NSTEPS= " << _nsteps << '\n';
    } else if( property == "ENDINPUT" ){
      coutf << "Reading input completed!" << '\n';
      break;
    } else cerr << "PROBLEM: unknown input" << '\n';
  }
  input.close();
  this->read_configuration();
  this->initialize_velocities();
  coutf << "System initialized!" << '\n';
  coutf.close();
  return;
}

void System :: initialize_velocities(){
  if(_restart and _sim_type==0){
    ifstream cinf;
    cinf.open("../INPUT/CONFIG/velocities.in");
    if(cinf.is_open()){
      double vx, vy, vz;
      for(int i=0; i<_npart; i++){
        cinf >> vx >> vy >> vz;
        _particle(i).setvelocity(0,vx);
        _particle(i).setvelocity(1,vy);
        _particle(i).setvelocity(2,vz);
      }
    } else cerr << "PROBLEM: Unable to open INPUT file velocities.in"<< '\n';
    cinf.close();
  } else {
    vec vx(_npart), vy(_npart), vz(_npart);
    vec sumv(_ndim);
    sumv.zeros();
    for (int i=0; i<_npart; i++){
      vx(i) = _rnd.Gauss(0.,sqrt(_temp));
      vy(i) = _rnd.Gauss(0.,sqrt(_temp));
      vz(i) = _rnd.Gauss(0.,sqrt(_temp));
      sumv(0) += vx(i);
      sumv(1) += vy(i);
      sumv(2) += vz(i);
    }
    for (int idim=0; idim<_ndim; idim++) sumv(idim) = sumv(idim)/double(_npart);
    double sumv2 = 0.0, scalef;
    for (int i=0; i<_npart; i++){
      vx(i) = vx(i) - sumv(0);
      vy(i) = vy(i) - sumv(1);
      vz(i) = vz(i) - sumv(2);
      sumv2 += vx(i) * vx(i) + vy(i) * vy(i) + vz(i) * vz(i);
    }
    sumv2 /= double(_npart);
    scalef = sqrt(3.0 * _temp / sumv2);   // velocity scale factor 
    for (int i=0; i<_npart; i++){
      _particle(i).setvelocity(0, vx(i)*scalef);
      _particle(i).setvelocity(1, vy(i)*scalef);
      _particle(i).setvelocity(2, vz(i)*scalef);
    }
  }
  if(_sim_type == 0){
  double xold, yold, zold;
    for (int i=0; i<_npart; i++){
      xold = this->pbc( _particle(i).getposition(0,true) - _particle(i).getvelocity(0)*_delta, 0);
      yold = this->pbc( _particle(i).getposition(1,true) - _particle(i).getvelocity(1)*_delta, 1);
      zold = this->pbc( _particle(i).getposition(2,true) - _particle(i).getvelocity(2)*_delta, 2);
      _particle(i).setpositold(0, xold);
      _particle(i).setpositold(1, yold);
      _particle(i).setpositold(2, zold);
    }
  }
  return;
}

void System :: initialize_properties(){ // Initialize data members used for measurement of properties

  string property;
  int index_property = 0;
  _nprop = 0;

  _measure_penergy  = false; //Defining which properties will be measured
  _measure_kenergy  = false;
  _measure_tenergy  = false;
  _measure_pressure = false;
  _measure_gofr     = false;
  _measure_magnet   = false;
  _measure_cv       = false;
  _measure_chi      = false;

  ifstream input;
  if(_sim_type < 2) input.open("../INPUT/properties.dat");
  else input.open("../INPUT/properties.ising");

  if (input.is_open()){
    while ( !input.eof() ){
      input >> property;
      if( property == "POTENTIAL_ENERGY" ){
        ofstream coutp("../OUTPUT/potential_energy.dat");
        coutp << "#     BLOCK:  ACTUAL_PE:     PE_AVE:      ERROR:" << '\n';
        coutp.close();
        _nprop++;
        _index_penergy = index_property;
        _measure_penergy = true;
        index_property++;
        _vtail = 0.0; // TO BE FIXED IN EXERCISE 7
      } else if( property == "KINETIC_ENERGY" ){
        ofstream coutk("../OUTPUT/kinetic_energy.dat");
        coutk << "#     BLOCK:   ACTUAL_KE:    KE_AVE:      ERROR:" << '\n';
        coutk.close();
        _nprop++;
        _measure_kenergy = true;
        _index_kenergy = index_property;
        index_property++;
      } else if( property == "TOTAL_ENERGY" ){
        ofstream coutt("../OUTPUT/total_energy.dat");
        coutt << "#     BLOCK:   ACTUAL_TE:    TE_AVE:      ERROR:" << '\n';
        coutt.close();
        _nprop++;
        _measure_tenergy = true;
        _index_tenergy = index_property;
        index_property++;
      } else if( property == "TEMPERATURE" ){
        ofstream coutte("../OUTPUT/temperature.dat");
        coutte << "#     BLOCK:   ACTUAL_T:     T_AVE:       ERROR:" << '\n';
        coutte.close();
        _nprop++;
        _measure_temp = true;
        _index_temp = index_property;
        index_property++;
      } else if( property == "PRESSURE" ){
        ofstream coutpr("../OUTPUT/pressure.dat");
        coutpr << "#     BLOCK:   ACTUAL_P:     P_AVE:       ERROR:" << '\n';
        coutpr.close();
        _nprop++;
        _measure_pressure = true;
        _index_pressure = index_property;
        index_property++;
        _ptail = 0.0; // TO BE FIXED IN EXERCISE 7
      } else if( property == "GOFR" ){
        ofstream coutgr("../OUTPUT/gofr.dat");
        coutgr << "# DISTANCE:     AVE_GOFR:        ERROR:" << '\n';
        coutgr.close();
        input>>_n_bins;
        _nprop+=_n_bins;
        _bin_size = (_halfside.min() )/(double)_n_bins;
        _measure_gofr = true;
        _index_gofr = index_property;
        index_property+= _n_bins;
      } else if( property == "MAGNETIZATION" ){
        ofstream coutpr("../OUTPUT/magnetization.dat");
        coutpr << "#     BLOCK:   ACTUAL_M:     M_AVE:       ERROR:" << '\n';
        coutpr.close();
        _nprop++;
        _measure_magnet = true;
        _index_magnet = index_property;
        index_property++;
      } else if( property == "SPECIFIC_HEAT" ){
        ofstream coutpr("../OUTPUT/specific_heat.dat");
        coutpr << "#     BLOCK:   ACTUAL_CV:    CV_AVE:      ERROR:" << '\n';
        coutpr.close();
        _nprop++;
        _measure_cv = true;
        _index_cv = index_property;
        index_property++;
      } else if( property == "SUSCEPTIBILITY" ){
        ofstream coutpr("../OUTPUT/susceptibility.dat");
        coutpr << "#     BLOCK:   ACTUAL_X:     X_AVE:       ERROR:" << '\n';
        coutpr.close();
        _nprop++;
        _measure_chi = true;
        _index_chi = index_property;
        index_property++;
      } else if( property == "ENDPROPERTIES" ){
        ofstream coutf;
        coutf.open("../OUTPUT/output.dat",ios::app);
        coutf << "Reading properties completed!" << '\n';
        coutf.close();
        break;
      } else cerr << "PROBLEM: unknown property" << '\n';
    }
    input.close();
  } else cerr << "PROBLEM: Unable to open properties.dat" << '\n';

  // according to the number of properties, resize the vectors _measurement,_average,_block_av,_global_av,_global_av2
  _measurement.resize(_nprop);
  _average.resize(_nprop);
  _block_av.resize(_nprop);
  _global_av.resize(_nprop);
  _global_av2.resize(_nprop);
  _average.zeros();
  _global_av.zeros();
  _global_av2.zeros();
  _nattempts = 0;
  _naccepted = 0;
  return;
}

void System :: finalize(){
  this->write_configuration();
  _rnd.SaveSeed();
  ofstream coutf;
  coutf.open("../OUTPUT/output.dat",ios::app);
  coutf << "Simulation completed!" << '\n';
  coutf.close();
  return;
}

// Write current configuration as a .xyz file in directory ../OUTPUT/CONFIG/
void System :: write_configuration(){
  ofstream coutf;
  if(_sim_type < 2){
    coutf.open("../OUTPUT/CONFIG/config.xyz");
    ofstream coutg{"../OUTPUT/CONFIG/config_old.xyz"};
    if(coutf.is_open()){
      coutf << _npart << '\n';
      coutf << "#Comment!" << '\n';
      coutg << _npart << '\n';
      coutg << "#Comment!" << '\n';
      for(int i=0; i<_npart; i++){
        coutf << "LJ" << "  " 
              << setw(16) << _particle(i).getposition(0,true)/_side(0)          // x
              << setw(16) << _particle(i).getposition(1,true)/_side(1)          // y
              << setw(16) << _particle(i).getposition(2,true)/_side(2) << '\n'; // z
        coutg << "LJ" << "  "
              << setw(16) << _particle(i).getposition(0,false)/_side(0)          // x
              << setw(16) << _particle(i).getposition(1,false)/_side(1)          // y
              << setw(16) << _particle(i).getposition(2,false)/_side(2) << '\n'; // z
      }
    } else cerr << "PROBLEM: Unable to open config.xyz" << '\n';
    coutf.close();
    this->write_velocities();
  } else {
    coutf.open("../OUTPUT/CONFIG/config.spin");
    for(int i=0; i<_npart; i++) coutf << _particle(i).getspin() << " ";
    coutf.close();
  }
  return;
}

// Write configuration nconf as a .xyz file in directory ../OUTPUT/CONFIG/
void System :: write_XYZ(int nconf){
  ofstream coutf;
  coutf.open("../OUTPUT/CONFIG/config_" + to_string(nconf) + ".xyz");
  if(coutf.is_open()){
    coutf << _npart << '\n';
    coutf << "#Comment!" << '\n';
    for(int i=0; i<_npart; i++){
      coutf << "LJ" << "  " 
            << setw(16) << _particle(i).getposition(0,true)          // x
            << setw(16) << _particle(i).getposition(1,true)          // y
            << setw(16) << _particle(i).getposition(2,true) << '\n'; // z
    }
  } else cerr << "PROBLEM: Unable to open config.xyz" << '\n';
  coutf.close();
  return;
}

void System :: write_velocities(){
  ofstream coutf;
  coutf.open("../OUTPUT/CONFIG/velocities.out");
  if(coutf.is_open()){
    for(int i=0; i<_npart; i++){
      coutf << setw(16) << _particle(i).getvelocity(0)          // vx
            << setw(16) << _particle(i).getvelocity(1)          // vy
            << setw(16) << _particle(i).getvelocity(2) << '\n'; // vz
    }
  } else cerr << "PROBLEM: Unable to open velocities.dat" << '\n';
  coutf.close();
  return;
}

// Read configuration from a .xyz file in directory ../OUTPUT/CONFIG/
void System :: read_configuration(){

  string file;
  if(_sim_type < 2) file = "xyz";
  else file = "ising";

  ifstream cinf;
  cinf.open("../INPUT/CONFIG/config."+file);

  if(cinf.is_open()){
    string comment;
    string particle;
    double x, y, z;
    int ncoord;
    cinf >> ncoord;
    
    if (ncoord != _npart){
      cerr << "PROBLEM: conflicting number of coordinates in input.dat & config."+ file +"not match!" << '\n';
      exit(EXIT_FAILURE);
    }

    cinf >> comment;

    if(_restart and _sim_type==0){//Check if the system is restarted from a previous run
      ifstream cing{"../INPUT/CONFIG/config_old.xyz"}; //Read the previous configuration
      double xold{}, yold{}, zold{};
      int ncoord_old;
      if(cing.is_open()){
        cing >> ncoord_old;
        if (ncoord_old != _npart){
          cerr << "PROBLEM: conflicting number of coordinates in input.dat & config_old.xyz not match!\n";
          exit(EXIT_FAILURE);
        }
        ofstream coutf{"../OUTPUT/output.dat",ios::app};
        coutf << "Restarting from previous run!\n";
        coutf.close();

        cing >> comment;
        for(int i{}; i<_npart; i++){ // Read the coordinates of the particles
          cinf >> particle >> x >> y >> z; // coordinates in conf.xyz are in units of _side
          cing >> particle >> xold >> yold >> zold; // coordinates in conf_old.xyz
          _particle(i).setposition(0, this->pbc(_side(0)*x, 0) ); //rescale the coordinates to the physical dimensions
          _particle(i).setposition(1, this->pbc(_side(1)*y, 1) );
          _particle(i).setposition(2, this->pbc(_side(2)*z, 2) );
          //_particle(i).acceptmove(); // _x_old = _x_new, don't need this because the old configuration is already stored in config_old.xyz
          _particle(i).setpositold(0, this->pbc(_side(0)*xold, 0) );
          _particle(i).setpositold(1, this->pbc(_side(1)*yold, 1) );
          _particle(i).setpositold(2, this->pbc(_side(2)*zold, 2) );
        }
      } else {
        cerr << "PROBLEM: Unable to open INPUT file config_old.xyz\n";
        exit(EXIT_FAILURE);
      }
      cinf.close();
      cing.close();
    } else {
      
      for(int i{}; i<_npart; i++){ // Read the coordinates of the particles
        cinf >> particle >> x >> y >> z; // units of coordinates in conf.xyz is _side
        _particle(i).setposition(0, this->pbc(_side(0)*x, 0) ); //rescale the coordinates to the physical dimensions
        _particle(i).setposition(1, this->pbc(_side(1)*y, 1) );
        _particle(i).setposition(2, this->pbc(_side(2)*z, 2) );
        _particle(i).acceptmove(); // _x_old = _x_new
      }
      cinf.close();
    }
  } else {
    cerr << "PROBLEM: Unable to open INPUT file config.xyz\n";
    exit(EXIT_FAILURE);
  }
  
  cinf.close();
  if(_restart and _sim_type > 1){
    int spin;
    cinf.open("../INPUT/CONFIG/config.spin");
    for(int i=0; i<_npart; i++){
      cinf >> spin;
      _particle(i).setspin(spin);
    }
    cinf.close();
  }
  return;
}

void System :: block_reset(int blk){ // Reset block accumulators to zero
  ofstream coutf;
  if(blk>0){
    coutf.open("../OUTPUT/output.dat",ios::app);
    coutf << "Block completed: " << blk << '\n';
    coutf.close();
  }
  _block_av.zeros();
  return;
}

void System :: measure(){ // Measure properties
  _measurement.zeros();
  // POTENTIAL ENERGY, VIRIAL, GOFR ///////////////////////////////////////////
  int bin{};
  vec distance{};
  distance.resize(_ndim);
  double penergy_temp{}, dr{}; // temporary accumulator for potential energy
  double kenergy_temp{}; // temporary accumulator for kinetic energy
  double tenergy_temp{}; //temporary accumulator for total energy
  double magnetization{};// temporary accumulator for magnetization
  double virial{};       // temporary accumulator for virial
  double vij{};         // LJ potential
  double appo{};
  if (_measure_penergy or _measure_pressure or _measure_gofr) {
    for (int i=0; i<_npart-1; i++){
      for (int j=i+1; j<_npart; j++){
        distance(0) = this->pbc( _particle(i).getposition(0,true) - _particle(j).getposition(0,true), 0);
        distance(1) = this->pbc( _particle(i).getposition(1,true) - _particle(j).getposition(1,true), 1);
        distance(2) = this->pbc( _particle(i).getposition(2,true) - _particle(j).getposition(2,true), 2);
        dr = sqrt( dot(distance,distance) );
        // GOFR ... TO BE FIXED IN EXERCISE 7
        if(dr < _r_cut){
          appo = 1.0/pow(dr,6); // to avoid multiple calculations
          vij = 1.0/pow(dr,12) - appo;
          if(_measure_penergy)  penergy_temp += vij; // LJ potential
          if(_measure_pressure) virial += vij + 0.5*appo; // VIRIAL
        }
      }
    }
  }
  // POTENTIAL ENERGY //////////////////////////////////////////////////////////
  if (_measure_penergy){
    penergy_temp = _vtail + 4.0 * penergy_temp / double(_npart);
    _measurement(_index_penergy) = penergy_temp;
  }
  // KINETIC ENERGY ////////////////////////////////////////////////////////////
  if (_measure_kenergy){
    for (int i=0; i<_npart; i++) kenergy_temp += 0.5 * dot( _particle(i).getvelocity() , _particle(i).getvelocity() ); 
    kenergy_temp /= double(_npart);
    _measurement(_index_kenergy) = kenergy_temp;
  }
  // TOTAL ENERGY (kinetic+potential) //////////////////////////////////////////
  if (_measure_tenergy){
    if (_sim_type < 2) _measurement(_index_tenergy) = kenergy_temp + penergy_temp;
    else { 
      double s_i, s_j;
      for (int i=0; i<_npart; i++){
        s_i = double(_particle(i).getspin());
        s_j = double(_particle(this->pbc(i+1)).getspin());
        tenergy_temp += - _J * s_i * s_j - 0.5 * _H * (s_i + s_j);
        magnetization += s_i; //sum of spins
      }
      _measurement(_index_tenergy) = tenergy_temp / double(_npart); // energy per particle
    }
  }
  // TEMPERATURE ///////////////////////////////////////////////////////////////
  if (_measure_temp and _measure_kenergy) _measurement(_index_temp) = (2.0/3.0) * kenergy_temp;

  // PRESSURE //////////////////////////////////////////////////////////////////
  if (_measure_pressure) _measurement(_index_pressure) = _rho * _measurement(_index_temp) + 48.0* virial/(3.0*_volume);
  
  // MAGNETIZATION /////////////////////////////////////////////////////////////
  if (_measure_magnet) _measurement(_index_magnet) = magnetization / double(_npart);
  
  // SPECIFIC HEAT /////////////////////////////////////////////////////////////
  if (_measure_cv) _measurement(_index_cv) = tenergy_temp * tenergy_temp; 
  
  // SUSCEPTIBILITY ////////////////////////////////////////////////////////////
  if (_measure_chi) _measurement(_index_chi) = _beta*(magnetization*magnetization) / double(_npart);
  
  _block_av += _measurement; //Update block accumulators

  return;
}

void System :: averages(int blk){

  ofstream coutf;
  double average, sum_average, sum_ave2;

  if(_measure_cv){
    //block_av(cv) = <(E)^2> - <E>^2
    _block_av(_index_cv) -=  std::pow(_block_av(_index_tenergy)*double(_npart),2) / double(_nsteps);
    _block_av(_index_cv) *= _beta*_beta;
    _block_av(_index_cv) /= double(_npart); // specific heat per particle
  }

  _average     = _block_av / double(_nsteps);
  _global_av  += _average;
  _global_av2 += _average % _average; // % -> element-wise multiplication

  // POTENTIAL ENERGY //////////////////////////////////////////////////////////
  if (_measure_penergy){
    coutf.open("../OUTPUT/potential_energy.dat",ios::app);
    average  = _average(_index_penergy);
    sum_average = _global_av(_index_penergy);
    sum_ave2 = _global_av2(_index_penergy);
    coutf << scientific
          << setw(12) << blk 
          << setw(18) << average
          << setw(18) << sum_average/double(blk)
          << setw(18) << this->error(sum_average, sum_ave2, blk) << '\n';
    coutf.close();
  }
  // KINETIC ENERGY ////////////////////////////////////////////////////////////
  if (_measure_kenergy){
    coutf.open("../OUTPUT/kinetic_energy.dat",ios::app);
    average  = _average(_index_kenergy);
    sum_average = _global_av(_index_kenergy);
    sum_ave2 = _global_av2(_index_kenergy);
    coutf << scientific
          << setw(12) << blk
          << setw(18) << average
          << setw(18) << sum_average/double(blk)
          << setw(18) << this->error(sum_average, sum_ave2, blk) << '\n';
    coutf.close();
  }
  // TOTAL ENERGY //////////////////////////////////////////////////////////////
  if (_measure_tenergy){
    coutf.open("../OUTPUT/total_energy.dat",ios::app);
    average  = _average(_index_tenergy);
    sum_average = _global_av(_index_tenergy);
    sum_ave2 = _global_av2(_index_tenergy);
    coutf << scientific
          << setw(12) << blk
          << setw(18) << average
          << setw(18) << sum_average/double(blk)
          << setw(18) << this->error(sum_average, sum_ave2, blk) << '\n';
    coutf.close();
  }
  // TEMPERATURE ///////////////////////////////////////////////////////////////
  if (_measure_temp){
    coutf.open("../OUTPUT/temperature.dat",ios::app);
    average  = _average(_index_temp);
    sum_average = _global_av(_index_temp);
    sum_ave2 = _global_av2(_index_temp);
    coutf << scientific
          << setw(12) << blk
          << setw(18) << average
          << setw(18) << sum_average/double(blk)
          << setw(18) << this->error(sum_average, sum_ave2, blk) << '\n';
    coutf.close();
  }
  // PRESSURE //////////////////////////////////////////////////////////////////
  if (_measure_pressure){
    coutf.open("../OUTPUT/pressure.dat",ios::app);
    average  = _average(_index_pressure);
    sum_average = _global_av(_index_pressure);
    sum_ave2 = _global_av2(_index_pressure);
    coutf << scientific
          << setw(12) << blk
          << setw(18) << average
          << setw(18) << sum_average/double(blk)
          << setw(18) << this->error(sum_average, sum_ave2, blk) << '\n';
    coutf.close();
  }
  // GOFR //////////////////////////////////////////////////////////////////////
  // TO BE FIXED IN EXERCISE 7
  // MAGNETIZATION /////////////////////////////////////////////////////////////
  if (_measure_magnet){
    coutf.open("../OUTPUT/magnetization.dat",ios::app);
    average  = _average(_index_magnet);
    sum_average = _global_av(_index_magnet);
    sum_ave2 = _global_av2(_index_magnet);
    coutf << scientific
          << setw(12) << blk
          << setw(18) << average
          << setw(18) << sum_average/double(blk)
          << setw(18) << this->error(sum_average, sum_ave2, blk) << '\n';
    coutf.close();
  }
  // SPECIFIC HEAT /////////////////////////////////////////////////////////////
  if (_measure_cv){
    coutf.open("../OUTPUT/specific_heat.dat",ios::app);
    average  = _average(_index_cv);
    sum_average = _global_av(_index_cv);
    sum_ave2 = _global_av2(_index_cv);
    coutf << scientific
          << setw(12) << blk
          << setw(18) << average
          << setw(18) << sum_average/double(blk)
          << setw(18) << this->error(sum_average, sum_ave2, blk) << '\n';
    coutf.close();
  }
  // SUSCEPTIBILITY ////////////////////////////////////////////////////////////
  if (_measure_chi){
    coutf.open("../OUTPUT/susceptibility.dat",ios::app);
    average  = _average(_index_chi);
    sum_average = _global_av(_index_chi);
    sum_ave2 = _global_av2(_index_chi);
    coutf << scientific
          << setw(12) << blk
          << setw(18) << average
          << setw(18) << sum_average/double(blk)
          << setw(18) << this->error(sum_average, sum_ave2, blk) << '\n';
    coutf.close();
  }
  // ACCEPTANCE ////////////////////////////////////////////////////////////////
  double fraction;
  coutf.open("../OUTPUT/acceptance.dat",ios::app);
  if(_nattempts > 0) fraction = double(_naccepted)/double(_nattempts);
  else fraction = 0.0; 
  coutf << setw(12) << blk << setw(18) << fraction << '\n';
  coutf.close();
  
  return;
}

double System :: error(double acc, double acc2, int blk){
  if(blk <= 1) return 0.0;
  else return sqrt( fabs(acc2/double(blk) - pow( acc/double(blk) ,2) )/double(blk) );
}

int System :: get_nbl(){
  return _nblocks;
}

int System :: get_nsteps(){
  return _nsteps;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
