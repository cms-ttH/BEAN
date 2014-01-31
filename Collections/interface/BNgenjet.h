#ifndef ProductArea_BNgenjet_h
#define ProductArea_BNgenjet_h

#include <vector>
#include <string>

#include "BNPhysicsObject.h"

namespace bean {
    struct GenJet : bean::PhysicsObject
    {
        explicit GenJet(double v,  int c) :
            pt(v), eta(v), phi(v),
            px(v), py(v), pz(v),
            et(v), energy(v), mass(v),
            emEnergy(v), hadEnergy(v),
            invisibleEnergy(v), auxiliaryEnergy(v),
            charge(c) {};

        GenJet() :
            pt(-99), eta(-99), phi(-99),
            px(-99), py(-99), pz(-99),
            et(-99), energy(-99), mass(-99),
            emEnergy(-99), hadEnergy(-99),
            invisibleEnergy(-99), auxiliaryEnergy(-99),
            charge(-99) { }

        double pt, eta, phi;
        double px, py, pz;
        double et, energy, mass;
        double emEnergy, hadEnergy;
        double invisibleEnergy, auxiliaryEnergy;
        int charge;

    };
}

typedef bean::GenJet BNgenjet;
typedef std::vector<BNgenjet> BNgenjetCollection;

#endif
