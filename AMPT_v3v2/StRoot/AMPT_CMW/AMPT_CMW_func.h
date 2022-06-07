#ifndef AMPT_CMW_FUNC_H
#define AMPT_CMW_FUNC_H

Double_t eff_function(Double_t* x,Double_t* par)
{
    Double_t pt,y;
    Double_t A,B,C;A=par[0];B=par[1];C=par[2];
    pt=x[0];
    y=A*(exp(-pow(B/pt,C)));
    return y;
}

#endif
