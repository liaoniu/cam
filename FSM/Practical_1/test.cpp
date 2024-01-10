#include <iostream>
#include <vector>
#include <array>
#include <fstream>
#include <sstream>
#include <cmath>
using namespace std;
using mat = vector<array<double, 4>>;   // u[i][var] = var[i]

array<double, 4> indi_primitiveToConserved(array<double, 4> &u_pri, string EoS);
array<double, 4> indi_conservedToPrimitive(array<double, 4> &u_con, string EoS);
mat conservedToPrimitive(mat &u, string EoS);
mat primitiveToConserved(mat &u, string EoS);
double cs(double &gamma, double &p, double &rho, string EoS);
double calc_dt(mat &u_con, double &dx, string EoS);
array<double, 4> f(array<double, 4> &u);
tuple<double, double> fK_and_dfK(double p_star, array<double, 4> &uK);
tuple<double, double> f_and_df(double p_star, array<double, 4> &uL, array<double, 4> &uR);
double compute_p_star(array<double, 4> &uL, array<double, 4> &uR, double epsilon);
double zeta(double delta_minus, double delta_plus, string limiter);
tuple<mat, mat> dataReconstruct(mat &u, string direction);
array<double, 4> getXFlux(array<double, 4> &uL, array<double, 4> &uR, double dx, double dt, string method, string EoS); //u_left and u_right
mat split(mat &u, mat &flux, double dt, double &dx, int &nxCells, string method, string EoS);
void transmissiveBC(mat& u, int nxCells);

void transmissiveBC(mat &u, int nxCells){
    for (int var = 0; var != 4; var++){
        u[0][var] = u[1][var];
        u[nxCells+1][var] = u[nxCells][var];
    }
}

array<double, 4> indi_primitiveToConserved(array<double, 4> &u_pri, string EoS){
    array<double, 4> u_con = u_pri;
    double rho, v, p, gamma;
    tie(rho, v, p, gamma) = tie(u_pri[0], u_pri[1], u_pri[2], u_pri[3]);
    double rhoV, epsilon, E, rhoGamma;
    double p_ref, e_ref, Gamma;
    if (EoS == "ideal_gas"){
        p_ref = 0;
        e_ref = 0;
        Gamma = gamma - 1;
    }
    else if (EoS == "stiffened_gas"){
        double p_infty = 3*1e8;
        p_ref = -gamma*p_infty;
        e_ref = 0;
        Gamma = gamma - 1;
    }

    rhoV = rho*v;
    epsilon = (p - p_ref)/Gamma/rho + e_ref;
    E = rho*epsilon + rho*(v*v)/2;
    rhoGamma = rho*gamma;
    tie(u_con[0], u_con[1], u_con[2], u_con[3]) = tie(rho, rhoV, E, rhoGamma);
    return u_con;
}

array<double, 4> indi_conservedToPrimitive(array<double, 4> &u_con, string EoS){
    array<double, 4> u_pri = u_con;
    double rho, rhoV, E, rhoGamma;
    tie(rho, rhoV, E, rhoGamma) = tie(u_con[0], u_con[1], u_con[2], u_con[3]);
    double v, epsilon, p;
    double p_ref, e_ref, Gamma;
    double gamma = rhoGamma/rho;
    if (EoS == "ideal_gas"){
        p_ref = 0;
        e_ref = 0;
        Gamma = gamma - 1;
    }
    else if (EoS == "stiffened_gas"){
        double p_infty = 3*1e8;
        p_ref = -gamma*p_infty;
        e_ref = 0;
        Gamma = gamma - 1;
    }

    v = rhoV/rho;
    epsilon = E/rho - 0.5*(v*v);
    p = p_ref + Gamma*rho*(epsilon - e_ref);
    tie(u_pri[0], u_pri[1], u_pri[2], u_pri[3]) = tie(rho, v, p, gamma);
    return u_pri;
}

mat conservedToPrimitive(mat &u, string EoS="ideal_gas"){
    //u would be changed
    for (int i = 0; i != u.size(); i++){
        u[i] = indi_conservedToPrimitive(u[i], EoS);
    }
    return u;
}

mat primitiveToConserved(mat &u, string EoS="ideal_gas"){
    //u would be changed
    for (int i = 0; i != u.size(); i++){
        u[i] = indi_primitiveToConserved(u[i], EoS);
    }
    return u;
}


double cs(double &gamma, double &p, double &rho, string EoS="ideal_gas"){
    if (EoS == "ideal_gas"){
        return sqrt(gamma*p/rho);
    }
    else if (EoS == "stiffened_gas"){
        double p_infty = 3*1e8;
        return sqrt(gamma*(p + p_infty)/rho);
    }
    return 0;
}

double calc_dt(mat &u_con, double &dx, string EoS="ideal_gas"){
    double delta = dx;
    double C = 0.8;
    double a_max = 0;
    int m = u_con.size();
    for (int i = 0; i != m; i++){
            double rho, v, p, gamma;
            array<double, 4> u_pri_indi = indi_conservedToPrimitive(u_con[i], EoS);
            tie(rho, v, p, gamma) = tie(u_pri_indi[0], u_pri_indi[1], u_pri_indi[2], u_pri_indi[3]);
            double a = sqrt(v*v) + cs(gamma, p, rho, EoS);
            a_max = ((a_max < a) ? a : a_max);
    }
    double dt = C*delta/a_max;
    return dt;
}

array<double, 4> f(array<double, 4> &u){
    double rho, v, p, gamma;
    double rhoV, E, rhoGamma;
    tie(rho, rhoV, E, rhoGamma) = tie(u[0], u[1], u[2], u[3]);
    v = rhoV/rho;
    gamma = rhoGamma/rho;
    p = (gamma - 1)*(E - 0.5*(rhoV*v));
    array<double, 4> f_res;
    f_res[0] = rhoV;
    f_res[1] = rhoV*v + p;
    f_res[2] = (E + p)*v;
    f_res[3] = rhoGamma*v;
    return f_res;
}



//Godunov Functions
tuple<double, double> fK_and_dfK(double p_star, array<double, 4> &uK){
    double rhoK, vK, pK, gamma;
    tie(rhoK, vK, pK, gamma) = tie(uK[0], uK[1], uK[2], uK[3]);
    double AK = 2/((gamma + 1)*rhoK);
    double BK = (gamma - 1)/(gamma + 1)*pK;
    double fK, dfK;
    if (p_star > pK){
        fK = (p_star - pK)*sqrt(AK/(p_star + BK));
        dfK = sqrt(AK/(BK + p_star))*(1 - (p_star - pK)/2/(BK + p_star));
    }
    else{
        double cs = sqrt(gamma*pK/rhoK);
        fK = 2*cs/(gamma - 1)*(pow(p_star/pK, (gamma - 1)/2/gamma) - 1);
        dfK = 1/rhoK/cs*pow(p_star/pK, -(gamma+1)/2/gamma);
    }
    return tie(fK, dfK);
}

tuple<double, double> f_and_df(double p_star, array<double, 4> &uL, array<double, 4> &uR){
    double f, df;
    double vL, vR, dv;
    tie(ignore, vL, ignore, ignore) = tie(uL[0], uL[1], uL[2], uL[3]);
    tie(ignore, vR, ignore, ignore) = tie(uR[0], uR[1], uR[2], uR[3]);
    dv = vR - vL;
    double fL, fR, dfL, dfR;
    tie(fL, dfL) = fK_and_dfK(p_star, uL);
    tie(fR, dfR) = fK_and_dfK(p_star, uR);
    f = fL + fR + dv;
    df = dfL + dfR;
    return tie(f, df);
}

double compute_p_star(array<double, 4> &uL, array<double, 4> &uR, double epsilon = 1e-3){
    double rhoL, vL, pL, rhoR, vR, pR, gammaL, gammaR;
    tie(rhoL, vL, pL, gammaL) = tie(uL[0], uL[1], uL[2], uL[3]);
    tie(rhoR, vR, pR, gammaR) = tie(uR[0], uR[1], uR[2], uR[3]);
    double csL = sqrt(gammaL*pL/rhoL);
    double csR = sqrt(gammaR*pR/rhoR);
    double p_init, p_PV;
    p_PV = (pR + pL)/2 - 1/8*(vR - vL)*(rhoL + rhoR)*(csL + csR);
    p_init = pow((csL+ csR - (gammaL - 1)*(vR - vL)/2)/(csL/pow(pL, (gammaL - 1)/2/gammaL) + csR/pow(pR, (gammaL - 1)/2/gammaL)), 2*gammaL/(gammaL - 1));
    double p_new = p_init;
    double p_old;
    double f, df;
    double err;
    do{
        p_old = p_new;
        tie(f, df) = f_and_df(p_new, uL, uR);
        p_new = p_old - f/df;
        err = fabs(p_new - p_old)/2/fabs(p_new + p_old);
    }while(err > epsilon);
    return p_new;
}

tuple<double, double, double, double, double, bool, bool> compute_stars_shock(array<double, 4> uL, array<double, 4> uR, double p_star){
    double v_star, rhoL_star, rhoR_star, SL = 0, SR = 0;
    double rhoL, vL, pL, rhoR, vR, pR, gammaL, gammaR;
    bool left_shock = false, right_shock = false;
    tie(rhoL, vL, pL, gammaL) = tie(uL[0], uL[1], uL[2], uL[3]);
    tie(rhoR, vR, pR, gammaR) = tie(uR[0], uR[1], uR[2], uR[3]);
    double csL, csR;
    csL = sqrt(gammaL*pL/rhoL);
    csR = sqrt(gammaR*pR/rhoR);
    if (p_star > pL){
        double gamma = gammaL;
        SL = vL - csL*sqrt((gamma+1)/2/gamma*p_star/pL + (gamma-1)/2/gamma);
        rhoL_star = rhoL*((p_star/pL + (gamma-1)/(gamma+1))/((gamma - 1)/(gamma+1)*(p_star/pL)+1));
        v_star = (1 - rhoL/rhoL_star)*SL + vL*rhoL/rhoL_star;
        left_shock = true;
    }
    else{
        double gamma = gammaL;
        rhoL_star = rhoL*pow(p_star/pL, 1/gamma);
        v_star = vL - 2*csL/(gamma - 1)*(pow(p_star/pL, (gamma - 1)/2/gamma) - 1);
    }
    if (p_star > pR){
        double gamma = gammaR;
        SR = vR + csR*sqrt((gamma+1)/2/gamma*p_star/pR + (gamma-1)/2/gamma);
        rhoR_star = rhoR*((p_star/pR + (gamma-1)/(gamma+1))/((gamma - 1)/(gamma+1)*(p_star/pR)+1));
        v_star = (1 - rhoR/rhoR_star)*SR + vR*rhoR/rhoR_star;
        right_shock = true;
    }
    else{
        double gamma = gammaR;
        rhoR_star = rhoR*pow(p_star/pR, 1/gamma);
        v_star = vR + 2*csR/(gamma - 1)*(pow(p_star/pR, (gamma - 1)/2/gamma) - 1);
    }
    return tie(v_star, rhoL_star, rhoR_star, SL, SR, left_shock, right_shock);
}

array<double, 4> primitiveInRarefaction(array<double, 4> uK, bool left){
    double rhoK, vK, pK, gamma;
    double rhoK_raref, vK_raref, pK_raref;
    tie(rhoK, vK, pK, gamma) = tie(uK[0], uK[1], uK[2], uK[3]);
    double csK = sqrt(gamma*pK/rhoK);
    if (left){
        rhoK_raref = rhoK*pow(2/(gamma+1) + (gamma-1)/(gamma+1)/csK*(vK), 2/(gamma - 1));
        vK_raref = 2/(gamma + 1)*(csK + (gamma - 1)/2*vK);
        pK_raref = pK*pow(2/(gamma+1) + (gamma-1)/(gamma+1)/csK*(vK), 2*gamma/(gamma - 1));
    }
    else{
        rhoK_raref = rhoK*pow(2/(gamma+1) - (gamma-1)/(gamma+1)/csK*(vK), 2/(gamma - 1));
        vK_raref = 2/(gamma + 1)*(-csK + (gamma - 1)/2*vK);
        pK_raref = pK*pow(2/(gamma+1) - (gamma-1)/(gamma+1)/csK*(vK), 2*gamma/(gamma - 1));
    }
    return array<double, 4> {rhoK_raref, vK_raref, pK_raref, gamma};
}

array<double, 4> calc_central_u(array<double, 4> uL, array<double, 4> uR, double p_star){
    double rhoL, vL, pL, rhoR, vR, pR, v_star, rhoL_star, rhoR_star, SL, SR, gammaL, gammaR;
    bool left_shock, right_shock;
    tie(rhoL, vL, pL, gammaL) = tie(uL[0], uL[1], uL[2], uL[3]);
    tie(rhoR, vR, pR, gammaR) = tie(uR[0], uR[1], uR[2], uR[3]);
    double csL = sqrt(gammaL*pL/rhoL);
    double csR = sqrt(gammaR*pR/rhoR);
    tie(v_star, rhoL_star, rhoR_star, SL, SR, left_shock, right_shock) = compute_stars_shock(uL, uR, p_star);
    double csL_star = sqrt(gammaL*p_star/rhoL_star);
    double csR_star = sqrt(gammaR*p_star/rhoR_star);
    //Compute gamma through Godunov
    double gamma;
    if (vL > vR){
        double S = (vR*gammaR - vL*gammaL)/(gammaR - gammaL);
        gamma = ((S > 0) ? gammaL : gammaR);
    }
    else{
        if (vL > 0)
            gamma = gammaL;
        else if (vR < 0)
            gamma = gammaR;
        else
            gamma = 0;
    }
    array<double, 4> u_mid;
    if (left_shock){
        if (SL > 0){
            u_mid = uL;
            u_mid[3] = gamma;
            return u_mid;
        }
    }
    else{
        if (vL - csL > 0){
            u_mid = uL;
            u_mid[3] = gamma;
            return u_mid;
        }
    }
    if (right_shock){
        if (SR < 0){
            u_mid = uR;
            u_mid[3] = gamma;
            return u_mid;
        }
    }
    else{
        if (vR + csR < 0){
            u_mid = uR;
            u_mid[3] = gamma;
            return u_mid;
        }
    }
    if (!left_shock && (v_star - csL_star > 0)){
        u_mid = primitiveInRarefaction(uL, true);
        u_mid[3] = gamma;
        return u_mid;
    }
    if (!right_shock && (v_star + csR_star < 0)){
        u_mid = primitiveInRarefaction(uR, false);
        u_mid[3] = gamma;
        return u_mid;
    }
    else{
        if (v_star < 0){
            return array<double, 4> {rhoR_star, v_star, p_star, gamma};
        }
        else{
            return array<double, 4> {rhoL_star, v_star, p_star, gamma};
        }
    }
    
}



array<double, 4> getXFlux(array<double, 4> &uL, array<double, 4> &uR, double dx, double dt, string method, string EoS){ //u_left and u_right
    array<double, 4> flux, u_mid;
    if (method == "Godunov"){
        double vL, pL, vR, pR, rho_mid, v_mid, p_mid;
        array<double, 4> uL_pri = indi_conservedToPrimitive(uL, EoS);
        array<double, 4> uR_pri = indi_conservedToPrimitive(uR, EoS);
        double p_star = compute_p_star(uL_pri, uR_pri);
        array<double, 4> u_mid;
        double gamma;
        u_mid = calc_central_u(uL_pri, uR_pri, p_star);
        u_mid = indi_primitiveToConserved(u_mid, EoS);
        flux = f(u_mid);
    }

    return flux;
}


mat split(mat &u, mat &flux, double dt, double &dx, int &nxCells, string method, string EoS){
    mat uL, uR;
    if (method == "SLIC"){
        // tie(uL, uR) = dataReconstruct(u);
        // halfTimeUpdate(uL, uR, gamma, dx, dt);
        // //getXFLux
        // for (int i = 0; i != nxCells+1; i++){
        //     flux[i] = getXFlux(uR[i], uL[i+1], gammaL, gammaR, dx, dt, "FORCE");
        // }
        // for (int i = 1; i != nxCells+1; i++){
        //     for (int k = 0; k != 3; k++){
        //         u[i][k] = u[i][k] - dt/dx*(flux[i][k] - flux[i-1][k]);
        //     }
        // }
    }
    else {
        //getXFLux
        for (int i = 0; i != nxCells+1; i++){
            flux[i] = getXFlux(u[i], u[i+1], dx, dt, method, EoS);
        }
        for (int i = 1; i != nxCells+1; i++){
            for (int k = 0; k != 4; k++){
                u[i][k] = u[i][k] - dt/dx*(flux[i][k] - flux[i-1][k]);
            }
        }
    }

    //Need ensure boundary condition here
    transmissiveBC(u, nxCells);
    return u;
}



void solver(mat &u0, double T, double x0, double x1, int nxCells = 100, string method = "FORCE", string EoS="ideal_gas"){
    double dx = (x1-x0)/nxCells;
    double t = 0;
    double dt = 0;
    mat flux;
    mat u = u0;
    flux.resize(nxCells+1);
    u = primitiveToConserved(u);

    do{
        dt = calc_dt(u, dx);
        t += dt;
        //Transmissive BC
        transmissiveBC(u, nxCells);
        mat u_temp = u;

        split(u, flux, dt, dx, nxCells, method, EoS);

    }while(t<T);

    u = conservedToPrimitive(u);
    ofstream output1("rho.dat");
    ofstream output2("v.dat");
    ofstream output3("p.dat");
    ofstream output4("gamma.dat");

    for (int i = 1; i != nxCells + 1; i++){
        double x = x0 + (i-0.5)*dx;
        output1 << x << " " << u[i][0] << endl;
        output2 << x << " " << u[i][1] << endl;
        output3 << x << " " << u[i][2] << endl;
        output4 << x << " " << u[i][3] << endl;
        }

    output1.close();
    output2.close();
    output3.close();
    output4.close();

}


int main()
{
    //Test 1
    double x0 = 0, x1 = 1, T = 0.005, gammaL = 1.4, gammaR = 1.67, nxCells = 100;
    double dx = (x1 - x0)/nxCells;
    double p_infty = 3*1e8;
    mat u0;
    u0.resize(nxCells+2);
    for (int i = 0; i != nxCells+2; i++){
        double x = x0 + (i-0.5)*dx;
        u0[i][0] = ((x < 0.5) ? 1 : 0.5);
        u0[i][1] = ((x < 0.5) ? 1 : 1);
        u0[i][2] = ((x < 0.5) ? 1 : 1);
        u0[i][3] = ((x < 0.5) ? gammaL : gammaR);
    }




    solver(u0, T, x0, x1, nxCells, "Godunov", "ideal_gas");
    return 0;

}