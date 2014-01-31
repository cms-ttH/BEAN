#ifndef ProductArea_BNmcparticle_h
#define ProductArea_BNmcparticle_h

#include <vector>

#include "BNPhysicsObject.h"

namespace bean {
    struct MCParticle : PhysicsObject
    {
        explicit MCParticle(double v,  int c) :
            energy(v), et(v), pt(v), px(v), py(v), pz(v), phi(v), eta(v), theta(v),
            mass(v), vx(v), vy(v), vz(v), motherET(v), motherPT(v), motherPhi(v),
            motherEta(v), mother0ET(v), mother0PT(v), mother0Phi(v), mother0Eta(v),
            mother1ET(v), mother1PT(v), mother1Phi(v), mother1Eta(v),
            daughter0ET(v), daughter0PT(v), daughter0Phi(v), daughter0Eta(v),
            daughter1ET(v), daughter1PT(v), daughter1Phi(v), daughter1Eta(v),
            grandMotherET(v), grandMotherPT(v), grandMotherPhi(v),
            grandMotherEta(v), grandMother00ET(v), grandMother00PT(v),
            grandMother00Phi(v), grandMother00Eta(v), grandMother01ET(v),
            grandMother01PT(v), grandMother01Phi(v), grandMother01Eta(v),
            grandMother10ET(v), grandMother10PT(v), grandMother10Phi(v),
            grandMother10Eta(v), grandMother11ET(v), grandMother11PT(v),
            grandMother11Phi(v), grandMother11Eta(v), charge(c), id(c), status(c),
            motherId(c), motherCharge(c), mother0Id(c), mother0Status(c),
            mother0Charge(c), mother1Id(c), mother1Status(c), mother1Charge(c),
            daughter0Id(c), daughter0Status(c), daughter0Charge(c), daughter1Id(c),
            daughter1Status(c), daughter1Charge(c), grandMotherId(c),
            grandMotherCharge(c), grandMother00Id(c), grandMother00Status(c),
            grandMother00Charge(c), grandMother01Id(c), grandMother01Status(c),
            grandMother01Charge(c), grandMother10Id(c), grandMother10Charge(c),
            grandMother11Id(c), grandMother11Charge(c) {};

        MCParticle() :
            energy(-99), et(-99), pt(-99), px(-99), py(-99), pz(-99), phi(-99),
            eta(-99), theta(-99), mass(-99), vx(-99), vy(-99), vz(-99),
            motherET(-99), motherPT(-99), motherPhi(-99), motherEta(-99),
            mother0ET(-99), mother0PT(-99), mother0Phi(-99), mother0Eta(-99),
            mother1ET(-99), mother1PT(-99), mother1Phi(-99), mother1Eta(-99),
            daughter0ET(-99), daughter0PT(-99), daughter0Phi(-99),
            daughter0Eta(-99), daughter1ET(-99), daughter1PT(-99),
            daughter1Phi(-99), daughter1Eta(-99), grandMotherET(-99),
            grandMotherPT(-99), grandMotherPhi(-99), grandMotherEta(-99),
            grandMother00ET(-99), grandMother00PT(-99), grandMother00Phi(-99),
            grandMother00Eta(-99), grandMother01ET(-99), grandMother01PT(-99),
            grandMother01Phi(-99), grandMother01Eta(-99), grandMother10ET(-99),
            grandMother10PT(-99), grandMother10Phi(-99), grandMother10Eta(-99),
            grandMother11ET(-99), grandMother11PT(-99), grandMother11Phi(-99),
            grandMother11Eta(-99), charge(-99), id(-99), status(-99),
            motherId(-99), motherCharge(-99), mother0Id(-99), mother0Status(-99),
            mother0Charge(-99), mother1Id(-99), mother1Status(-99),
            mother1Charge(-99), daughter0Id(-99), daughter0Status(-99),
            daughter0Charge(-99), daughter1Id(-99), daughter1Status(-99),
            daughter1Charge(-99), grandMotherId(-99), grandMotherCharge(-99),
            grandMother00Id(-99), grandMother00Status(-99),
            grandMother00Charge(-99), grandMother01Id(-99),
            grandMother01Status(-99), grandMother01Charge(-99),
            grandMother10Id(-99), grandMother10Status(-99),
            grandMother10Charge(-99), grandMother11Id(-99),
            grandMother11Status(-99), grandMother11Charge(-99) {};

        double energy, et, pt, px, py, pz;
        double phi, eta, theta, mass;
        double vx, vy, vz;
        double motherET, motherPT, motherPhi, motherEta;
        double mother0ET, mother0PT, mother0Phi, mother0Eta;
        double mother1ET, mother1PT, mother1Phi, mother1Eta;
        double daughter0ET, daughter0PT, daughter0Phi, daughter0Eta;
        double daughter1ET, daughter1PT, daughter1Phi, daughter1Eta;
        double grandMotherET, grandMotherPT, grandMotherPhi, grandMotherEta;
        double grandMother00ET, grandMother00PT, grandMother00Phi, grandMother00Eta;
        double grandMother01ET, grandMother01PT, grandMother01Phi, grandMother01Eta;
        double grandMother10ET, grandMother10PT, grandMother10Phi, grandMother10Eta;
        double grandMother11ET, grandMother11PT, grandMother11Phi, grandMother11Eta;
        int charge, id, status;
        int motherId, motherCharge;
        int mother0Id, mother0Status, mother0Charge;
        int mother1Id, mother1Status, mother1Charge;
        int daughter0Id, daughter0Status, daughter0Charge;
        int daughter1Id, daughter1Status, daughter1Charge;
        int grandMotherId, grandMotherCharge;
        int grandMother00Id, grandMother00Status, grandMother00Charge;
        int grandMother01Id, grandMother01Status, grandMother01Charge;
        int grandMother10Id, grandMother10Status, grandMother10Charge;
        int grandMother11Id, grandMother11Status, grandMother11Charge;
    };
}

typedef bean::MCParticle BNmcparticle;
typedef std::vector<BNmcparticle> BNmcparticleCollection;

#endif
