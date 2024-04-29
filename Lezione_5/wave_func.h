#pragma once

#include <fstream>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>

#include "RndGen/random.h"

class WaveFunctionSimulation {
private:
    Random rnd;
    int nstep;
    int nblk;
    int steps;
    double L;
    bool gauss;
    int np_rend;
    int mod;
    int contastep;
    std::string tr_pr;

    // Variables for Ground State
    double d;
    double var_prog_d;
    double mean_prog_d;
    double r[3]; //initial position
    double y[3]; //position before move
    double x[3]; //position after move
    double dr_g[3]; //proposed displacement
    double q;
    double A;
    int attempted_gs;
    int accepted_gs;

    // Variables for Excited State
    double de;
    double var_prog_de;
    double mean_prog_de;
    double re[3]; //initial position
    double ye[3]; //position before move
    double xe[3]; //position after move
    double dr_e[3]; //proposed displacement
    double qe;
    double Ae;
    int attempted_es;
    int accepted_es;

    // Private helper functions
    void ProposeStep();
    double prob_gs(double x[3]);
    double prob_2P0(double x[3]);
    double min(double s, double t);

public:
    WaveFunctionSimulation();
    ~WaveFunctionSimulation();
    void Input();
    void ResetAll();
    void SavePos();
    void Move();
    void Accumulate();
    void PrintAccRate(int blknum);
    void BlockAverages();
    void SaveDist(int blknum);
    double error(double av, double av2, int n);
    int getnblk() { return nblk; }
    int getnstep() { return nstep; }
};

