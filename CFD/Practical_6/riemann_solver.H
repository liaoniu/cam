#include "D:\me\cam\CFD\eulerSolver1D_v2.H"

// Exact Riemann Solver Routine



tuple<double, double> fK_and_dfK(double p_star, double rhoK, double vK, double pK, double csK)
{
    double fK, dfK;
    double AK = 2/(Gamma+1)/rhoK;
    double BK = (Gamma-1)/(Gamma+1)*pK;
    if (p_star > pK)
    {
        fK = (p_star - pK)*sqrt(AK/(p_star + BK));
        dfK = sqrt(AK/(BK + p_star))*(1 - (p_star - pK)/2/(BK + p_star));
    }
    else{
        fK = 2*csK/(Gamma - 1)*(pow(p_star/pK, (Gamma - 1)/2/Gamma) - 1);
        dfK = 1/rhoK/csK*pow(p_star/pK, -(Gamma+1)/2/Gamma);
    }
    return tie(fK, dfK);
}


double compute_p_star(const arr &uL, const arr &uR, double TOL=1e-3)
{
    // Note! Inputs takes primitive form
    double rhoL, vL, pL, csL;
    double rhoR, vR, pR, csR;
    tie(rhoL, vL, pL) = tie(uL[0], uL[1], uL[2]);
    tie(rhoR, vR, pR) = tie(uR[0], uR[1], uR[2]);
    double dv = vR - vL;
    csL = computeSoundSpeedFromEoS(rhoL, pL);
    csR = computeSoundSpeedFromEoS(rhoR, pR);
    double p_TR = pow((csL + csR - 0.5*(Gamma - 1)*(vR - vL)) / (csL/pow(pL, (Gamma - 1)/2/Gamma) + csR/pow(pR, (Gamma - 1)/2/Gamma)), 2*Gamma/(Gamma - 1));
    
    double p_old, p_new;
    p_old = p_TR;
    p_new = p_old;
    double f, df, err;
    do
    {
        p_old = p_new;
        double fL, fR, dfL, dfR;
        tie(fL, dfL) = fK_and_dfK(p_old, rhoL, vL, pL, csL);
        tie(fR, dfR) = fK_and_dfK(p_old, rhoR, vR, pR, csR);
        f = fL + fR + dv;
        df = dfL + dfR;
        p_new = p_old - f/df;
        err = fabs(p_new - p_old)/2/fabs(p_new + p_old);
    } while (err > TOL);

    return p_new;
}

tuple<arr, double> compute_inter_shock(double p_star, const arr &uK, bool left)
{
    double rhoK, vK, pK, SK;
    tie(rhoK, vK, pK) = tie(uK[0], uK[1], uK[2]);
    double csK = computeSoundSpeedFromEoS(rhoK, pK);
    double rhoK_star = rhoK*(p_star/pK + (Gamma-1)/(Gamma+1))/((Gamma-1)/(Gamma+1)*p_star/pK + 1);
    if (left)
        SK = vK - csK*sqrt((Gamma+1)/2/Gamma*p_star/pK + (Gamma-1)/2/Gamma);
    else
        SK = vK + csK*sqrt((Gamma+1)/2/Gamma*p_star/pK + (Gamma-1)/2/Gamma);
    double vK_star = (1 - rhoK/rhoK_star)*SK + vK*rhoK/rhoK_star;
    arr uK_star = {rhoK_star, vK_star, p_star};

    return tie(uK_star, SK);
}

arr compute_inter_raref(double p_star, const arr &uK, bool left)
{
    double rhoK, vK, pK, SK;
    tie(rhoK, vK, pK) = tie(uK[0], uK[1], uK[2]);
    double csK = computeSoundSpeedFromEoS(rhoK, pK);
    double rhoK_star = rhoK*pow(p_star/pK, 1/Gamma);
    double vK_star;
    if (left)
        vK_star = vK - 2*csK/(Gamma-1)*(pow(p_star/pK, (Gamma-1)/2/Gamma) - 1);
    else
        vK_star = vK + 2*csK/(Gamma-1)*(pow(p_star/pK, (Gamma-1)/2/Gamma) - 1);
    arr uK_star = {rhoK_star, vK_star, p_star};

    return uK_star;
}

arr primitiveInRarefaction(const arr &uK, bool left){
    double rhoK, vK, pK;
    double rhoK_raref, vK_raref, pK_raref;
    tie(rhoK, vK, pK) = tie(uK[0], uK[1], uK[2]);
    double csK = computeSoundSpeedFromEoS(rhoK, pK);
    if (left){
        rhoK_raref = rhoK*pow(2/(Gamma+1) + (Gamma-1)/(Gamma+1)/csK*(vK), 2/(Gamma - 1));
        vK_raref = 2/(Gamma + 1)*(csK + (Gamma - 1)/2*vK);
        pK_raref = pK*pow(2/(Gamma+1) + (Gamma-1)/(Gamma+1)/csK*(vK), 2*Gamma/(Gamma - 1));
    }
    else{
        rhoK_raref = rhoK*pow(2/(Gamma+1) - (Gamma-1)/(Gamma+1)/csK*(vK), 2/(Gamma - 1));
        vK_raref = 2/(Gamma + 1)*(-csK + (Gamma - 1)/2*vK);
        pK_raref = pK*pow(2/(Gamma+1) - (Gamma-1)/(Gamma+1)/csK*(vK), 2*Gamma/(Gamma - 1));
    }
    arr uK_raref = {rhoK_raref, vK_raref, pK_raref};
    return uK_raref;
}




arr exactFindUmid(const arr &uL, const arr &uR)
{
    // Note! Inputs takes primitive form
    // Return Also Primitive Form
    arr u_mid;
    arr uL_star, uR_star;
    double SL, SR;
    double p_star = compute_p_star(uL, uR, 1e-3);
    double rhoL, vL, pL, csL, csL_star;
    double rhoR, vR, pR, csR, csR_star;
    tie(rhoL, vL, pL) = tie(uL[0], uL[1], uL[2]);
    tie(rhoR, vR, pR) = tie(uR[0], uR[1], uR[2]);
    csL = computeSoundSpeedFromEoS(rhoL, pL);
    csR = computeSoundSpeedFromEoS(rhoR, pR);
    bool leftShock = false, rightShock = false;
    if (p_star > pL)
        leftShock = true;
    if (p_star > pR)
        rightShock = true;
    
    if (leftShock && rightShock)
    {
        tie(uL_star, SL) = compute_inter_shock(p_star, uL, true);
        tie(uR_star, SR) = compute_inter_shock(p_star, uR, false);
        if (SL > 0)
            return uL;
        if (SR < 0)
            return uR;
        if (uL_star[1] > 0)
            return uL_star;
        else
            return uR_star;
    }

    if (leftShock && !rightShock)
    {
        tie(uL_star, SL) = compute_inter_shock(p_star, uL, true);
        uR_star = compute_inter_raref(p_star, uR, false);
        if (SL > 0)
            return uL;
        if (vR + csR < 0)
            return uR;
        if (uL_star[1] > 0)
            return uL_star;
        csR_star = computeSoundSpeedFromEoS(uR_star[0], uR_star[2]);
        if (uR_star[1] + csR_star > 0)
            return uR_star;
        else
            return primitiveInRarefaction(uR, false);
    }

    if (!leftShock && rightShock)
    {
        uL_star = compute_inter_raref(p_star, uL, true);
        tie(uR_star, SR) = compute_inter_shock(p_star, uR, false);
        if (SR < 0)
            return uR;
        if (vL - csL > 0)
            return uL;
        if (uR_star[1] < 0)
            return uR_star;
        csL_star = computeSoundSpeedFromEoS(uL_star[0], uL_star[2]);
        if (uL_star[1] - csL_star < 0)
            return uL_star;
        else
            return primitiveInRarefaction(uL, true);
    }

    else
    {
        uL_star = compute_inter_raref(p_star, uL, true);
        uR_star = compute_inter_raref(p_star, uR, false);
        csL_star = computeSoundSpeedFromEoS(uL_star[0], uL_star[2]);
        csR_star = computeSoundSpeedFromEoS(uR_star[0], uR_star[2]);
        if (vL - csL > 0)
            return uL;
        if (vR + csR < 0)
            return uR;
        if (uL_star[1] - csL_star > 0)
            return primitiveInRarefaction(uL, true);
        if (uR_star[1] + csR_star < 0)
            return primitiveInRarefaction(uR, false);
        if (uL_star[1] > 0)
            return uL_star;
        else
            return uR_star;
    }
}




