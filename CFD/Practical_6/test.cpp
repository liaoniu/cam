#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
using namespace std;
using mat = vector<vector<double>>;
using myTuple = tuple<vector<double>, vector<double>, vector<double>>;

tuple<double, double> fK_and_dfK(double p_star, tuple<double, double, double> uK, double gamma){
    double rhoK, vK, pK;
    tie(rhoK, vK, pK) = uK;
    double AK = 2/((gamma + 1)*rhoK);
    double BK = (gamma - 1)/(gamma + 1)*pK;
    double fK, dfK;
    if (p_star > pK){
        fK = (p_star - pK)*sqrt(AK/(p_star + BK));
        dfK = pow(AK/(BK + p_star), 1/2)*(1 - (p_star - pK)/2/(BK + p_star));
    }
    else{
        double cs = sqrt(gamma*pK/rhoK);
        fK = 2*cs/(gamma - 1)*(pow(p_star/pK, (gamma - 1)/2/gamma) - 1);
        dfK = 1/rhoK/cs*pow(p_star/pK, -(gamma+1)/2/gamma);
    }
    return tie(fK, dfK);
}

tuple<double, double> f_and_df(double p_star, tuple<double, double, double> uL, tuple<double, double, double> uR, double gamma){
    double f, df;
    double vL, vR, dv;
    tie(ignore, vL, ignore) = uL;
    tie(ignore, vR, ignore) = uR;
    dv = vR - vL;
    double fL, fR, dfL, dfR;
    tie(fL, dfL) = fK_and_dfK(p_star, uL, gamma);
    tie(fR, dfR) = fK_and_dfK(p_star, uR, gamma);
    f = fL + fR + dv;
    df = dfL + dfR;
    return tie(f, df);
}


tuple<double, double, double, double, double, bool, bool> compute_stars_shock(tuple<double, double, double> uL, tuple<double, double, double> uR, double gamma, double p_star){
    double v_star, rhoL_star, rhoR_star, SL = 0, SR = 0;
    double rhoL, vL, pL, rhoR, vR, pR;
    bool left_shock = false, right_shock = false;
    tie(rhoL, vL, pL) = uL;
    tie(rhoR, vR, pR) = uR;
    double csL, csR;
    csL = sqrt(gamma*pL/rhoL);
    csR = sqrt(gamma*pR/rhoR);
    if (p_star > pL){
        SL = vL - csL*sqrt((gamma+1)/2/gamma*p_star/pL + (gamma-1)/2/gamma);
        rhoL_star = rhoL*((p_star/pL + (gamma-1)/(gamma+1))/((gamma - 1)/(gamma+1)*(p_star/pL)+1));
        v_star = (1 - rhoL/rhoL_star)*SL + vL*rhoL/rhoL_star;
        left_shock = true;
    }
    else{
        rhoL_star = rhoL*pow(p_star/pL, 1/gamma);
        v_star = vL - 2*csL/(gamma - 1)*(pow(p_star/pL, (gamma - 1)/2/gamma) - 1);
    }
    if (p_star > pR){
        SR = vR + csR*sqrt((gamma+1)/2/gamma*p_star/pR + (gamma-1)/2/gamma);
        rhoR_star = rhoR*((p_star/pR + (gamma-1)/(gamma+1))/((gamma - 1)/(gamma+1)*(p_star/pR)+1));
        v_star = (1 - rhoR/rhoR_star)*SR + vR*rhoR/rhoR_star;
        right_shock = true;
    }
    else{
        rhoR_star = rhoR*pow(p_star/pR, 1/gamma);
        v_star = vR + 2*csR/(gamma - 1)*(pow(p_star/pR, (gamma - 1)/2/gamma) - 1);
    }

    return tie(v_star, rhoL_star, rhoR_star, SL, SR, left_shock, right_shock);
}

tuple<double, double, double> primitiveInRarefaction(double x, double t, double gamma, tuple<double, double, double> uK, bool left){
    double rhoK, vK, pK;
    double rhoK_raref, vK_raref, pK_raref;
    tie(rhoK, vK, pK) = uK;
    double csK = sqrt(gamma*pK/rhoK);
    if (left){
        rhoK_raref = rhoK*pow(2/(gamma+1) - (gamma-1)/(gamma+1)/csK*(vK - x/t), 2/(gamma - 1));
        vK_raref = 2/(gamma + 1)*(-csK + (gamma - 1)/2*vK + x/t);
        pK_raref = pK*pow(2/(gamma+1) - (gamma-1)/(gamma+1)/csK*(vK - x/t), 2*gamma/(gamma - 1));
    }
    else{
        rhoK_raref = rhoK*pow(2/(gamma+1) + (gamma-1)/(gamma+1)/csK*(vK - x/t), 2/(gamma - 1));
        vK_raref = 2/(gamma + 1)*(-csK + (gamma - 1)/2*vK + x/t);
        pK_raref = pK*pow(2/(gamma+1) + (gamma-1)/(gamma+1)/csK*(vK - x/t), 2*gamma/(gamma - 1));
    }
    return tie(rhoK_raref, vK_raref, pK_raref);
}

tuple<double, double, double> calc_central_u(tuple<double, double, double> uL, tuple<double, double, double> uR, double gamma, double p_star){
    double rhoL, vL, pL, rhoR, vR, pR, v_star, rhoL_star, rhoR_star, SL, SR;
    bool left_shock, right_shock;
    tie(rhoL, vL, pL) = uL;
    tie(rhoR, vR, pR) = uR;
    double csL = sqrt(gamma*pL/rhoL);
    double csR = sqrt(gamma*pR/rhoR);
    tie(v_star, rhoL_star, rhoR_star, SL, SR, left_shock, right_shock) = compute_stars_shock(uL, uR, gamma, p_star);
    double csL_star = sqrt(gamma*p_star/rhoL_star);
    double csR_star = sqrt(gamma*p_star/rhoR_star);
    if (left_shock){
        if (SL > 0){
            return tie(rhoL, vL, pL);
        }
    }
    else{
        if (vL - csL > 0){
            return tie(rhoL, vL, pL);
        }
    }
    if (right_shock){
        if (SR < 0){
            return tie(rhoR, vR, pR);
        }
    }
    else{
        if (vR + csR < 0){
            return tie(rhoR, vR, pR);
        }
    }
    if (!left_shock && (v_star - csL_star > 0)){
        return primitiveInRarefaction(0, 0, gamma, uL, true);
    }
    if (!right_shock && (v_star + csR_star < 0)){
        return primitiveInRarefaction(0, 0, gamma, uR, false);
    }
    else{
        if (v_star < 0){
            return tie(rhoR_star, v_star, p_star);
        }
        else{
            return tie(rhoL_star, v_star, p_star);
        }
    }
    
}



double compute_p_star(tuple<double, double, double> uL, tuple<double, double, double> uR, double gamma, double epsilon = 1e-10){
    double rhoL, vL, pL, rhoR, vR, pR;
    tie(rhoL, vL, pL) = uL;
    tie(rhoR, vR, pR) = uR;
    double csL = sqrt(gamma*pL/rhoL);
    double csR = sqrt(gamma*pR/rhoR);
    double p_init = pow((csL+ csR - (gamma - 1)*(vR - vL)/2)/(csL/pow(pL, (gamma - 1)/2/gamma) + csR/pow(pR, (gamma - 1)/2/gamma)), 2*gamma/(gamma - 1));
    double p_new = p_init;
    double p_old;
    double f, df;
    do{
        p_old = p_new;
        tie(f, df) = f_and_df(p_new, uL, uR, gamma);
        p_new = p_new - f/df;
    }while(fabs(p_new - p_old)/2/fabs(p_new + p_old) > epsilon);
    return p_new;
}

int main()
{
    double rhoL, vL, pL, rhoR, vR, pR;
    rhoL = 1;
    vL = -2;
    pL = 0.4;
    rhoR = 1;
    vR = 2;
    pR = 0.4;
    tuple<double, double, double> uL = tie(rhoL, vL, pL);
    tuple<double, double, double> uR = tie(rhoR, vR, pR);
    double gamma = 1.4;
    double p_star = compute_p_star(uL, uR, gamma);
    cout << p_star << endl;

    double v_star, rhoL_star, rhoR_star, SL, SR;
    bool left_shock, right_shock;
    tie(v_star, rhoL_star, rhoR_star, SL, SR, left_shock, right_shock) = compute_stars_shock(uL, uR, gamma, p_star);

    double fL, fR, dfL, dfR;
    tie(fL, dfL) = fK_and_dfK(p_star, uL, gamma);
    tie(fR, dfR) = fK_and_dfK(p_star, uR, gamma);
    double vv = (vL+vR)/2 + (fR - fL)/2;

    double rho_mid, v_mid, p_mid;

    tie(rho_mid, v_mid, p_mid) = calc_central_u(uL, uR, gamma, p_star);


    return 0;



}