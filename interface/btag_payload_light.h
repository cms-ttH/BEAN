

float mistag_CSVL(float eta,float x, float scale) {
  float mean=1.,min=1.,max=1.;
  double absEta = fabs(eta);

  if( absEta < 0.5 ){
    min = ((0.994425+(-8.66392e-05*x))+(-3.03813e-08*(x*x)))+(-3.52151e-10*(x*(x*x)));
    max = ((1.15628+(0.000437668*x))+(-1.69625e-06*(x*x)))+(1.00718e-09*(x*(x*x)));
    mean = ((1.07536+(0.000175506*x))+(-8.63317e-07*(x*x)))+(3.27516e-10*(x*(x*x)));
  }
  else if( absEta < 1.0 && absEta > 0.5 ){
    min = ((0.998088+(6.94916e-05*x))+(-4.82731e-07*(x*x)))+(1.63506e-10*(x*(x*x)));
    max = ((1.15882+(0.000579711*x))+(-2.12243e-06*(x*x)))+(1.53771e-09*(x*(x*x)));
    mean = ((1.07846+(0.00032458*x))+(-1.30258e-06*(x*x)))+(8.50608e-10*(x*(x*x)));
  }
  else if( absEta < 1.5 && absEta > 1.0 ){
    min = ((1.00294+(0.000289844*x))+(-7.9845e-07*(x*x)))+(5.38525e-10*(x*(x*x)));
    max = ((1.16292+(0.000659848*x))+(-2.07868e-06*(x*x)))+(1.72763e-09*(x*(x*x)));
    mean = ((1.08294+(0.000474818*x))+(-1.43857e-06*(x*x)))+(1.13308e-09*(x*(x*x)));
  }
  else if( absEta < 2.4 && absEta > 1.5 ){
    min = ((0.979816+(0.000138797*x))+(-3.14503e-07*(x*x)))+(2.38124e-10*(x*(x*x)));
    max = ((1.14357+(0.00020854*x))+(-7.43519e-07*(x*x)))+(8.73742e-10*(x*(x*x)));
    mean = ((1.0617+(0.000173654*x))+(-5.29009e-07*(x*x)))+(5.55931e-10*(x*(x*x)));
  }

  if( x>670 ){
    min = 0.925136;
    max = 1.09972;
    mean = 1.0124;
  }

  float SFl = mean;
  if( scale>0 ) SFl = mean + scale*(max-mean);
  else if( scale<0 ) SFl = mean + scale*(mean-min);

  return SFl*(0.979396 + 0.000205898*x + 2.49868e-07*x*x);
}

float mistag_CSVM(float eta,float x, float scale) {
  float mean=1.,min=1.,max=1.;
  double absEta = fabs(eta);

  if( absEta < 0.8 ){
    min = ((0.972455+(7.51396e-06*x))+(4.91857e-07*(x*x)))+(-1.47661e-09*(x*(x*x)));
    max = ((1.15116+(0.00122657*x))+(-3.63826e-06*(x*x)))+(2.08242e-09*(x*(x*x)));
    mean = ((1.06182+(0.000617034*x))+(-1.5732e-06*(x*x)))+(3.02909e-10*(x*(x*x)));
  }
  else if( absEta < 1.6 && absEta > 0.8 ){
    min = ((1.02055+(-0.000378856*x))+(1.49029e-06*(x*x)))+(-1.74966e-09*(x*(x*x)));
    max = ((1.20146+(0.000359543*x))+(-1.12866e-06*(x*x)))+(6.59918e-10*(x*(x*x)));
    mean = ((1.111+(-9.64191e-06*x))+(1.80811e-07*(x*x)))+(-5.44868e-10*(x*(x*x)));
  }
  else if( absEta < 2.4 && absEta > 1.6 ){
    min = ((0.983476+(-0.000607242*x))+(3.17997e-06*(x*x)))+(-4.01242e-09*(x*(x*x)));
    max = ((1.18654+(-0.000795808*x))+(3.69226e-06*(x*x)))+(-4.22347e-09*(x*(x*x)));
    mean = ((1.08498+(-0.000701422*x))+(3.43612e-06*(x*x)))+(-4.11794e-09*(x*(x*x)));
  }

  if( x>670 ){
    min = 0.844346;
    max = 1.05012;
    mean = 0.947232;
  }

  float SFl = mean;
  if( scale>0 ) SFl = mean + scale*(max-mean);
  else if( scale<0 ) SFl = mean + scale*(mean-min);

  return SFl*(1.10422 + -0.000523856*x + 1.14251e-06*x*x);
}

float mistag_CSVT(float eta,float x, float scale) {
  float mean=1.,min=1.,max=1.;
  double absEta = fabs(eta);

  if( absEta < 2.4 ){
    min = ((0.899715+(0.00102278*x))+(-2.46335e-06*(x*x)))+(9.71143e-10*(x*(x*x)));
    max = ((0.997077+(0.00473953*x))+(-1.34985e-05*(x*x)))+(1.0032e-08*(x*(x*x)));
    mean = ((0.948463+(0.00288102*x))+(-7.98091e-06*(x*x)))+(5.50157e-09*(x*(x*x)));
  }

  if( x>670 ){
    min = 0.771264;
    max = 1.13034;
    mean = 0.950785;
  }

  float SFl = mean;
  if( scale>0 ) SFl = mean + 1.5*scale*(max-mean);
  else if( scale<0 ) SFl = mean + 1.5*scale*(mean-min);

  return SFl*(1.19275 + -0.00191042*x + 2.92205e-06*x*x);
}
