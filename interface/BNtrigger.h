#ifndef ProductArea_BNtrigger_h
#define ProductArea_BNtrigger_h

#include <vector>
#include <string>

// a simple class
struct BNtrigger
{
  explicit BNtrigger(int c, std::string n):pass(c),prescale(c),name(n) { }
  BNtrigger():pass(-99),prescale(-99),name("empty") { }
  int    pass,prescale;
  std::string name;
};

// this is our new product, it is simply a 
// collection of BNtrigger held in an std::vector
typedef std::vector<BNtrigger> BNtriggerCollection;

#endif
