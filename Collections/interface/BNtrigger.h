#ifndef ProductArea_BNtrigger_h
#define ProductArea_BNtrigger_h

#include <vector>
#include <string>

struct BNtrigger
{
  explicit BNtrigger(int c, std::string n) :
     pass(c), prescale(c), name(n) {};
  BNtrigger() :
     pass(-99), prescale(-99), name("empty") {};

  int pass, prescale;
  std::string name;
};

typedef std::vector<BNtrigger> BNtriggerCollection;

#endif
