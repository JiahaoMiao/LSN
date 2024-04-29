#include "wave_func.h"

WaveFunctionSimulation::WaveFunctionSimulation() : rnd(), nstep(0), nblk(0), steps(0), L(0), gauss(false), np_rend(0), mod(0), contastep(0),
                                                   tr_pr(""), d(0), var_prog_d(0), mean_prog_d(0),r{},y{},x{},dr_g{}, q(0), A(0), attempted_gs(0), accepted_gs(0),
                                                   de(0), var_prog_de(0), mean_prog_de(0), re{},ye{},xe{},dr_e{}, qe(0), Ae(0), attempted_es(0), accepted_es(0) {}
WaveFunctionSimulation::~WaveFunctionSimulation() {}

void WaveFunctionSimulation::Input() {
    // Read input informations
    std::ifstream ReadInput{"input.dat"};
    std::ofstream fout{"output.dat"};
    fout << ",-------------------------------------," << '\n';
    fout << "|            Hydrogen atom            |" << '\n';
    fout << "| Monte Carlo simulation (Metropolis) |" << '\n';
    fout << "'-------------------------------------'" << '\n';

    // Read seed for random numbers
    rnd.initialize();

    ReadInput >> nblk;  // number of blocks
    ReadInput >> nstep; // number of steps per block
    ReadInput >> L;     // step length
    ReadInput >> r[0];  // initial position
    ReadInput >> re[0]; // initial position
    ReadInput >> gauss; // transition probability
    ReadInput >> mod;  // save to file every mod steps

    if (!gauss)
        tr_pr = "unif";
    else if (gauss)
        tr_pr = "gauss";

    steps = nstep * nblk; // total number of steps = steps per block * number of blocks
    
    fout << "Total number of steps = " << steps << '\n';
    fout << "Number of blocks = " << nblk << '\n';
    fout << "Number of steps in one block = " << nstep << '\n';
    fout << "Step length = " << L << '\n';
    fout << "Initial position (GS) = " << r[0] << '\n';
    fout << "Initial position (ExSt) = " << re[0] << '\n';

    ReadInput.close();
    fout.close();
}

void WaveFunctionSimulation::ResetAll() {
    // reset distance
    d = 0;
    de = 0;
    // reset attempted moves
    attempted_gs = 0; 
    attempted_es = 0;
    // reset accepted moves
    accepted_gs = 0;
    accepted_es = 0;
}

void WaveFunctionSimulation::SavePos() {
    std::ofstream WritePosGS{"GS/pos_" + tr_pr + ".dat", std::ios::app};
    std::ofstream WritePosES{"ES/pos_" + tr_pr + ".dat", std::ios::app};

    // Save current position
    for (int k = 0; k < 3; k++) {
        y[k] = r[k];   // GS
        ye[k] = re[k]; // ES

        // Also save to file every mod steps
        if (contastep % mod == 0) {
            WritePosGS << y[k] << " ";
            WritePosES << ye[k] << " ";
        }
    }
    // newline every mod steps
    if (contastep % mod == 0) {
        WritePosGS << '\n';
        WritePosES << '\n';
    }

    WritePosGS.close();
    WritePosES.close();
}

void WaveFunctionSimulation::ProposeStep() {
    if (!gauss) {
        for (int k = 0; k < 3; k++) {
            dr_g[k] = rnd.Rannyu(-L, L); // Uniform transition probability in a cube of side 2L
            dr_e[k] = rnd.Rannyu(-L, L); // Same for excited state
        }
    } else if (gauss) {
        for (int k = 0; k < 3; k++) {
            dr_g[k] = rnd.Gauss(0, L); // Gaussian transition probability with sigma = L
            dr_e[k] = rnd.Gauss(0, L); // Same for excited state
        }
    }
}

void WaveFunctionSimulation::Move() {
    ProposeStep();

    // Try to move
    contastep++;
    for (int k = 0; k < 3; k++) {
        x[k] = y[k] + dr_g[k];   // GS
        xe[k] = ye[k] + dr_e[k]; // ES
    }

    // Calculate probability and accept according to Metropolis algorithm

    // GS
    attempted_gs++;
    q = prob_gs(x) / prob_gs(y);
    A = min(1, q); // Metropolis acceptance probability

    if (A == 1) {
        for (int k = 0; k < 3; k++) r[k] = x[k];
        accepted_gs++;
    } else {
        if (rnd.Rannyu() < A) { //accept with probability A
            for (int k = 0; k < 3; k++) r[k] = x[k];
            accepted_gs++;
        }
    }

    // ES
    attempted_es++;
    qe = prob_2P0(xe) / prob_2P0(ye);
    Ae = min(1, qe);

    if (Ae == 1) {
        for (int k = 0; k < 3; k++) re[k] = xe[k];
        accepted_es++;
    } else {
        if (rnd.Rannyu() < Ae) {
            for (int k = 0; k < 3; k++) re[k] = xe[k];
            accepted_es++;
        }
    }
}

void WaveFunctionSimulation::Accumulate() {
    d += sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
    de += sqrt(re[0] * re[0] + re[1] * re[1] + re[2] * re[2]);
}

void WaveFunctionSimulation::PrintAccRate(int i) {
    std::ofstream fout{"output.dat", std::ios::app};
    fout << "Block # " << i + 1 << '\n';
    fout << "                    (GrSt)  (ExSt)  " << '\n';
    fout << "Acceptance rates:   " << (double)accepted_gs / attempted_gs
         << "  " << (double)accepted_es / attempted_es << '\n';
    fout << "-----------------------------------" << '\n';
}

void WaveFunctionSimulation::BlockAverages() {
    d /= nstep;
    mean_prog_d += d;
    var_prog_d += d * d;

    de /= nstep;
    mean_prog_de += de;
    var_prog_de += de * de;
}

void WaveFunctionSimulation::SaveDist(int i) {
    std::ofstream WriteResultGS;
    WriteResultGS.open("GS/dist_" + tr_pr + ".dat", std::ios::app);
    std::ofstream WriteResultES;
    WriteResultES.open("ES/dist_" + tr_pr + ".dat", std::ios::app);

    // save to file
    if (WriteResultGS.is_open()) {
        if (i == 0)
            WriteResultGS << d << " " << mean_prog_d / (i + 1) << " " << 0 << '\n';
        else
            WriteResultGS << d << " " << mean_prog_d / (i + 1) << " "
                          << error(mean_prog_d / (i + 1), var_prog_d / (i + 1), i)
                          << " " << '\n';

    } else {
        std::cerr << "PROBLEM: Unable to open resultGS.dat" << '\n';
    }

    if (WriteResultES.is_open()) {
        if (i == 0)
            WriteResultES << de << " " << mean_prog_de / (i + 1) << " " << 0 << '\n';
        else
            WriteResultES << de << " " << mean_prog_de / (i + 1) << " "
                          << error(mean_prog_de / (i + 1), var_prog_de / (i + 1), i)
                          << " " << '\n';

    } else {
        std::cerr << "PROBLEM: Unable to open resultES.dat" << '\n';
    }

    WriteResultGS.close();
    WriteResultES.close();
}

double WaveFunctionSimulation::error(double av, double av2, int n) {
    return n == 0 ? 0.0 : sqrt((av2 - av * av) / n);
}

double WaveFunctionSimulation::prob_gs(double x[3]) {
    double d = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
    double psi = pow(M_E, -d) / sqrt(M_PI);
    return psi * psi;
}

double WaveFunctionSimulation::prob_2P0(double x[3]) {
    double d = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
    double costheta = x[2] / d;
    double psi = 1. / 8. * sqrt(2. / M_PI) * d * pow(M_E, -d / 2) * costheta;
    return psi * psi;
}

double WaveFunctionSimulation::min(double s, double t) {
    return s <= t ? s : t;
}
