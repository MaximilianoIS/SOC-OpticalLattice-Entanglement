#include <armadillo>
#include <fstream>
using namespace arma;

struct OATPars { double Jpar, Jperp, JOmg; };

OATPars map_xxz_to_oat(double J, double U, double Omega, double phi){
    const double c  = std::cos(phi);
    const double c2 = std::cos(0.5*phi);             
    const double U2 = U*U;
    const double O2 = Omega*Omega;
    const double den = U2 - O2;

    OATPars p;
    p.Jpar  = (J*J/U) * ( (U2*c - O2*c2*c2) / den );          
    p.Jperp = (J*J/U) * ( (U2     - O2*c2*c2) / den );        
    p.JOmg  = 0.5*Omega;                                      
    return p;
}

struct DickeOps { cx_mat Sx, Sy, Sz; };

DickeOps build_dicke_ops(int L){
    const int dim = L + 1;
    const double S = 0.5 * L;
    cx_mat Sp(dim,dim,fill::zeros), Sm(dim,dim,fill::zeros), Sz(dim,dim,fill::zeros);

    for(int m=-int(S); m<int(S); ++m){
        double c = std::sqrt(S*(S+1) - m*(m+1));
        Sp(m+1+int(S), m+int(S)) = c;
        Sm(m+int(S),   m+1+int(S)) = c;
    }
    for(int m=-int(S); m<=int(S); ++m) Sz(m+int(S), m+int(S)) = double(m);

    DickeOps ops;
    ops.Sx = 0.5*(Sp + Sm);
    ops.Sy = cx_double(0.0, -0.5) * (Sp - Sm);
    ops.Sz = Sz;
    return ops;
}

cx_mat build_H_OAT(int L, double Jpar, double Jperp, double JOmega, const DickeOps& ops){
    const double iso = 2.0*(Jpar + Jperp)/(L - 1);   
    const double tw  = 2.0*(Jpar - Jperp)/(L - 1);   

    cx_mat S2 = ops.Sx*ops.Sx + ops.Sy*ops.Sy + ops.Sz*ops.Sz;
    cx_mat H(ops.Sx.n_rows, ops.Sx.n_cols, fill::zeros);
    H += iso * S2;
    H += tw * (ops.Sz * ops.Sz);
    H += 2.0 * JOmega * ops.Sx;
    return H;
}

// Spins up along z
cx_vec dicke_all_up(int L) {
    int dim = L + 1;           
    cx_vec psi(dim);           
    for (int i = 0; i < dim; ++i)
        psi(i) = cx_double(0.0, 0.0);  
    psi(dim - 1) = cx_double(1.0, 0.0); //last basis vector
    return psi;
}

// All spins in the +x direction
cx_vec dicke_all_right(int L, const DickeOps& ops) {
    cx_vec psi_up = dicke_all_up(L);

    double angle = -datum::pi / 2.0;     
    cx_mat Ry = expmat(cx_double(0, -angle) * ops.Sy);

    cx_vec psi_right = Ry * psi_up;

    return psi_right;
}

int main(){
    int    L     = 8;            
    double J     = 1.0;
    double U     = 50.0;         
    double phi   = 0.5;
    double Omega = 15.0 * J; 
    double dt    = 0.001;

    auto pars = map_xxz_to_oat(J, U, Omega, phi);
    double tmax = 8.0 / fabs(pars.Jpar - pars.Jperp);
    DickeOps ops = build_dicke_ops(L);
    cx_mat H = build_H_OAT(L, pars.Jpar, pars.Jperp, pars.JOmg, ops);

    cx_vec psi0 = dicke_all_up(L);

    vec evals; cx_mat evecs;
    eig_sym(evals, evecs, H);
    cx_vec coeffs = evecs.t() * psi0;

    std::ofstream fout("contrast_oat.csv");
    fout << "t,Cx,Cz\n";
    const cx_double I(0,1);

    for(double t=0.0; t<=tmax+1e-12; t+=dt){
        cx_vec phases = exp(-I * evals * t);
        cx_vec psi_t  = evecs * (coeffs % phases);

        double Ex = real(cdot(psi_t, ops.Sx*psi_t));
        double Ey = real(cdot(psi_t, ops.Sy*psi_t));
        double Ez = real(cdot(psi_t, ops.Sz*psi_t));

        double Cx = (2.0/double(L)) * Ex;
        double Cz = (2.0/double(L)) * std::hypot(Ey, Ez);

        fout << t << "," << Cx << "," << Cz << "\n";
    }
    fout.close();
    return 0;
}

