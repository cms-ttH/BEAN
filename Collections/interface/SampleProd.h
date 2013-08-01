#ifndef ProductArea_SampleProd_h
#define ProductArea_SampleProd_h

#include <vector>


// a simple class
struct SampleProd
{
  explicit SampleProd(float v):et_(v),pt_(v),px_(v),py_(v),pz_(v),phi_(v),eta_(v),theta_(v) { }
  SampleProd():et_(-99),pt_(-99),px_(-99),py_(-99),pz_(-99),phi_(-99),eta_(-99),theta_(-99) { }
  float et_,pt_,px_,py_,pz_,phi_,eta_,theta_;
};

// this is our new product, it is simply a 
// collection of SampleProd held in an std::vector
typedef std::vector<SampleProd> SampleCollection;

#endif
