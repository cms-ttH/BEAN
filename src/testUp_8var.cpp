#include "testUp_8var.hpp"
#include <cmath>

double testUp_8var::Value(int index,double in0,double in1,double in2,double in3,double in4,double in5,double in6,double in7) {
   input0 = (in0 - -0.00165308)/1.61685;
   input1 = (in1 - 0.00171885)/1.73727;
   input2 = (in2 - 138.351)/132.525;
   input3 = (in3 - 2.7945)/1.11845;
   input4 = (in4 - 100.997)/115.907;
   input5 = (in5 - 2.75303)/1.36459;
   input6 = (in6 - 80.953)/100.095;
   input7 = (in7 - 2.65524)/1.4474;
   switch(index) {
     case 0:
         return neuron0xca141e0();
     default:
         return 0.;
   }
}

double testUp_8var::Value(int index, double* input) {
   input0 = (input[0] - -0.00165308)/1.61685;
   input1 = (input[1] - 0.00171885)/1.73727;
   input2 = (input[2] - 138.351)/132.525;
   input3 = (input[3] - 2.7945)/1.11845;
   input4 = (input[4] - 100.997)/115.907;
   input5 = (input[5] - 2.75303)/1.36459;
   input6 = (input[6] - 80.953)/100.095;
   input7 = (input[7] - 2.65524)/1.4474;
   switch(index) {
     case 0:
         return neuron0xca141e0();
     default:
         return 0.;
   }
}

double testUp_8var::neuron0xca20410() {
   return input0;
}

double testUp_8var::neuron0xca20720() {
   return input1;
}

double testUp_8var::neuron0xca20a30() {
   return input2;
}

double testUp_8var::neuron0xca20d40() {
   return input3;
}

double testUp_8var::neuron0xca21050() {
   return input4;
}

double testUp_8var::neuron0xc991190() {
   return input5;
}

double testUp_8var::neuron0xc9914a0() {
   return input6;
}

double testUp_8var::neuron0xc9917e0() {
   return input7;
}

double testUp_8var::input0xc991b20() {
   double input = -0.223505;
   input += synapse0xc9a46a0();
   input += synapse0xc991e30();
   input += synapse0xc991e70();
   input += synapse0xc991eb0();
   input += synapse0xc991ef0();
   input += synapse0xc991f30();
   input += synapse0xc991f70();
   input += synapse0xc991fb0();
   return input;
}

double testUp_8var::neuron0xc991b20() {
   double input = input0xc991b20();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double testUp_8var::input0xca1d510() {
   double input = 0.190823;
   input += synapse0xca1d820();
   input += synapse0xca1d860();
   input += synapse0xca1d8a0();
   input += synapse0xca1d8e0();
   input += synapse0xca1d920();
   input += synapse0xca1d960();
   input += synapse0xca1d9a0();
   input += synapse0xca1d9e0();
   return input;
}

double testUp_8var::neuron0xca1d510() {
   double input = input0xca1d510();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double testUp_8var::input0xca1da20() {
   double input = 0.29544;
   input += synapse0xca1dd30();
   input += synapse0xca29030();
   input += synapse0x10106170();
   input += synapse0x101061b0();
   input += synapse0xc991ff0();
   input += synapse0xc992030();
   input += synapse0xc992070();
   input += synapse0xc9920b0();
   return input;
}

double testUp_8var::neuron0xca1da20() {
   double input = input0xca1da20();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double testUp_8var::input0xca1de80() {
   double input = -0.269691;
   input += synapse0xca1e190();
   input += synapse0xca1e1d0();
   input += synapse0xca1e210();
   input += synapse0xca1e250();
   input += synapse0xca1e290();
   input += synapse0xca1e2d0();
   input += synapse0xca1e310();
   input += synapse0xca1e350();
   return input;
}

double testUp_8var::neuron0xca1de80() {
   double input = input0xca1de80();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double testUp_8var::input0xca1e390() {
   double input = 0.404981;
   input += synapse0xca1e6a0();
   input += synapse0xca390f0();
   input += synapse0xc9a4520();
   input += synapse0xca1dd70();
   input += synapse0xca1ddb0();
   input += synapse0xca1ddf0();
   input += synapse0xca1de30();
   input += synapse0xc9a4830();
   return input;
}

double testUp_8var::neuron0xca1e390() {
   double input = input0xca1e390();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double testUp_8var::input0xca1e8f0() {
   double input = 0.0166029;
   input += synapse0xca1ec00();
   input += synapse0xca1ec40();
   input += synapse0xca1ec80();
   input += synapse0xca1ecc0();
   input += synapse0xca1ed00();
   input += synapse0xca1ed40();
   input += synapse0xca1ed80();
   input += synapse0xca1edc0();
   return input;
}

double testUp_8var::neuron0xca1e8f0() {
   double input = input0xca1e8f0();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double testUp_8var::input0xca1ee00() {
   double input = -0.431461;
   input += synapse0xca1f110();
   input += synapse0xca1f150();
   input += synapse0xca1f190();
   input += synapse0xca1f1d0();
   input += synapse0xca1f210();
   input += synapse0xca1f250();
   input += synapse0xca1f290();
   input += synapse0xca1f2d0();
   return input;
}

double testUp_8var::neuron0xca1ee00() {
   double input = input0xca1ee00();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double testUp_8var::input0xca11740() {
   double input = -0.252423;
   input += synapse0xca11a50();
   input += synapse0xca11a90();
   input += synapse0xca11ad0();
   input += synapse0xca11b10();
   input += synapse0xca11b50();
   input += synapse0xca11b90();
   input += synapse0xca11bd0();
   input += synapse0xca11c10();
   return input;
}

double testUp_8var::neuron0xca11740() {
   double input = input0xca11740();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double testUp_8var::input0xca11c50() {
   double input = -0.0721609;
   input += synapse0x10105fc0();
   input += synapse0x10106000();
   input += synapse0xc9a4870();
   input += synapse0xc9a48b0();
   input += synapse0xc9a48f0();
   input += synapse0xc9a4930();
   input += synapse0xc9a4970();
   input += synapse0xca1f310();
   return input;
}

double testUp_8var::neuron0xca11c50() {
   double input = input0xca11c50();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double testUp_8var::input0xca1e6e0() {
   double input = -0.386897;
   input += synapse0xca125a0();
   input += synapse0xca125e0();
   input += synapse0xca12620();
   input += synapse0xca12660();
   input += synapse0xca126a0();
   input += synapse0xca126e0();
   input += synapse0xca12720();
   input += synapse0xca12760();
   return input;
}

double testUp_8var::neuron0xca1e6e0() {
   double input = input0xca1e6e0();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double testUp_8var::input0xca127a0() {
   double input = 0.32746;
   input += synapse0xca12ae0();
   input += synapse0xca12b20();
   input += synapse0xca12b60();
   input += synapse0xca12ba0();
   input += synapse0xca12be0();
   input += synapse0xca12c20();
   input += synapse0xca12c60();
   input += synapse0xca12ca0();
   return input;
}

double testUp_8var::neuron0xca127a0() {
   double input = input0xca127a0();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double testUp_8var::input0xca12ce0() {
   double input = 0.417442;
   input += synapse0xca13020();
   input += synapse0xca13060();
   input += synapse0xca130a0();
   input += synapse0xca130e0();
   input += synapse0xca13120();
   input += synapse0xca13160();
   input += synapse0xca131a0();
   input += synapse0xca131e0();
   return input;
}

double testUp_8var::neuron0xca12ce0() {
   double input = input0xca12ce0();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double testUp_8var::input0xca13220() {
   double input = -0.458196;
   input += synapse0xca13560();
   input += synapse0xca135a0();
   input += synapse0xca135e0();
   input += synapse0xca13620();
   input += synapse0xca13660();
   input += synapse0xca136a0();
   input += synapse0xca136e0();
   input += synapse0xca13720();
   return input;
}

double testUp_8var::neuron0xca13220() {
   double input = input0xca13220();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double testUp_8var::input0xca13760() {
   double input = 0.316729;
   input += synapse0xca13aa0();
   input += synapse0xca13ae0();
   input += synapse0xca13b20();
   input += synapse0xca13b60();
   input += synapse0xca13ba0();
   input += synapse0xca13be0();
   input += synapse0xca13c20();
   input += synapse0xca13c60();
   return input;
}

double testUp_8var::neuron0xca13760() {
   double input = input0xca13760();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double testUp_8var::input0xca13ca0() {
   double input = -0.172563;
   input += synapse0xca13fe0();
   input += synapse0xca14020();
   input += synapse0xca14060();
   input += synapse0xca140a0();
   input += synapse0xca140e0();
   input += synapse0xca14120();
   input += synapse0xca14160();
   input += synapse0xca141a0();
   return input;
}

double testUp_8var::neuron0xca13ca0() {
   double input = input0xca13ca0();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double testUp_8var::input0xca141e0() {
   double input = 0.23802;
   input += synapse0xca14520();
   input += synapse0xca14560();
   input += synapse0xca145a0();
   input += synapse0xca145e0();
   input += synapse0xca14620();
   input += synapse0xca14660();
   input += synapse0xca146a0();
   input += synapse0xca146e0();
   input += synapse0xca14720();
   input += synapse0xca1e870();
   input += synapse0xca1e8b0();
   input += synapse0xca1f3e0();
   input += synapse0xca1f420();
   input += synapse0xca29160();
   input += synapse0xca291a0();
   return input;
}

double testUp_8var::neuron0xca141e0() {
   double input = input0xca141e0();
   return (input * 1)+0;
}

double testUp_8var::synapse0xc9a46a0() {
   return (neuron0xca20410()*0.224733);
}

double testUp_8var::synapse0xc991e30() {
   return (neuron0xca20720()*0.0693316);
}

double testUp_8var::synapse0xc991e70() {
   return (neuron0xca20a30()*0.124538);
}

double testUp_8var::synapse0xc991eb0() {
   return (neuron0xca20d40()*-0.329444);
}

double testUp_8var::synapse0xc991ef0() {
   return (neuron0xca21050()*-0.265944);
}

double testUp_8var::synapse0xc991f30() {
   return (neuron0xc991190()*-0.504027);
}

double testUp_8var::synapse0xc991f70() {
   return (neuron0xc9914a0()*0.126659);
}

double testUp_8var::synapse0xc991fb0() {
   return (neuron0xc9917e0()*0.414205);
}

double testUp_8var::synapse0xca1d820() {
   return (neuron0xca20410()*0.0875955);
}

double testUp_8var::synapse0xca1d860() {
   return (neuron0xca20720()*0.290021);
}

double testUp_8var::synapse0xca1d8a0() {
   return (neuron0xca20a30()*-0.0514226);
}

double testUp_8var::synapse0xca1d8e0() {
   return (neuron0xca20d40()*-0.244169);
}

double testUp_8var::synapse0xca1d920() {
   return (neuron0xca21050()*-0.336007);
}

double testUp_8var::synapse0xca1d960() {
   return (neuron0xc991190()*0.158296);
}

double testUp_8var::synapse0xca1d9a0() {
   return (neuron0xc9914a0()*0.0388821);
}

double testUp_8var::synapse0xca1d9e0() {
   return (neuron0xc9917e0()*0.376108);
}

double testUp_8var::synapse0xca1dd30() {
   return (neuron0xca20410()*-0.230289);
}

double testUp_8var::synapse0xca29030() {
   return (neuron0xca20720()*0.0938741);
}

double testUp_8var::synapse0x10106170() {
   return (neuron0xca20a30()*0.332573);
}

double testUp_8var::synapse0x101061b0() {
   return (neuron0xca20d40()*0.302418);
}

double testUp_8var::synapse0xc991ff0() {
   return (neuron0xca21050()*0.186268);
}

double testUp_8var::synapse0xc992030() {
   return (neuron0xc991190()*-0.0521467);
}

double testUp_8var::synapse0xc992070() {
   return (neuron0xc9914a0()*-0.0999771);
}

double testUp_8var::synapse0xc9920b0() {
   return (neuron0xc9917e0()*-0.191849);
}

double testUp_8var::synapse0xca1e190() {
   return (neuron0xca20410()*-0.00236579);
}

double testUp_8var::synapse0xca1e1d0() {
   return (neuron0xca20720()*-0.364589);
}

double testUp_8var::synapse0xca1e210() {
   return (neuron0xca20a30()*-0.337518);
}

double testUp_8var::synapse0xca1e250() {
   return (neuron0xca20d40()*0.521269);
}

double testUp_8var::synapse0xca1e290() {
   return (neuron0xca21050()*-0.32081);
}

double testUp_8var::synapse0xca1e2d0() {
   return (neuron0xc991190()*0.281799);
}

double testUp_8var::synapse0xca1e310() {
   return (neuron0xc9914a0()*0.054788);
}

double testUp_8var::synapse0xca1e350() {
   return (neuron0xc9917e0()*0.411501);
}

double testUp_8var::synapse0xca1e6a0() {
   return (neuron0xca20410()*-0.365537);
}

double testUp_8var::synapse0xca390f0() {
   return (neuron0xca20720()*-0.183734);
}

double testUp_8var::synapse0xc9a4520() {
   return (neuron0xca20a30()*0.35888);
}

double testUp_8var::synapse0xca1dd70() {
   return (neuron0xca20d40()*0.388761);
}

double testUp_8var::synapse0xca1ddb0() {
   return (neuron0xca21050()*0.228885);
}

double testUp_8var::synapse0xca1ddf0() {
   return (neuron0xc991190()*-0.0481367);
}

double testUp_8var::synapse0xca1de30() {
   return (neuron0xc9914a0()*-0.0111447);
}

double testUp_8var::synapse0xc9a4830() {
   return (neuron0xc9917e0()*0.264687);
}

double testUp_8var::synapse0xca1ec00() {
   return (neuron0xca20410()*-0.294037);
}

double testUp_8var::synapse0xca1ec40() {
   return (neuron0xca20720()*0.41738);
}

double testUp_8var::synapse0xca1ec80() {
   return (neuron0xca20a30()*0.262689);
}

double testUp_8var::synapse0xca1ecc0() {
   return (neuron0xca20d40()*-0.265291);
}

double testUp_8var::synapse0xca1ed00() {
   return (neuron0xca21050()*-0.368972);
}

double testUp_8var::synapse0xca1ed40() {
   return (neuron0xc991190()*-0.359903);
}

double testUp_8var::synapse0xca1ed80() {
   return (neuron0xc9914a0()*-0.305203);
}

double testUp_8var::synapse0xca1edc0() {
   return (neuron0xc9917e0()*0.440525);
}

double testUp_8var::synapse0xca1f110() {
   return (neuron0xca20410()*0.0139094);
}

double testUp_8var::synapse0xca1f150() {
   return (neuron0xca20720()*0.156326);
}

double testUp_8var::synapse0xca1f190() {
   return (neuron0xca20a30()*0.449083);
}

double testUp_8var::synapse0xca1f1d0() {
   return (neuron0xca20d40()*0.26697);
}

double testUp_8var::synapse0xca1f210() {
   return (neuron0xca21050()*0.317291);
}

double testUp_8var::synapse0xca1f250() {
   return (neuron0xc991190()*0.275406);
}

double testUp_8var::synapse0xca1f290() {
   return (neuron0xc9914a0()*0.262371);
}

double testUp_8var::synapse0xca1f2d0() {
   return (neuron0xc9917e0()*0.177175);
}

double testUp_8var::synapse0xca11a50() {
   return (neuron0xca20410()*0.0346303);
}

double testUp_8var::synapse0xca11a90() {
   return (neuron0xca20720()*0.151766);
}

double testUp_8var::synapse0xca11ad0() {
   return (neuron0xca20a30()*-0.369438);
}

double testUp_8var::synapse0xca11b10() {
   return (neuron0xca20d40()*0.318608);
}

double testUp_8var::synapse0xca11b50() {
   return (neuron0xca21050()*0.374113);
}

double testUp_8var::synapse0xca11b90() {
   return (neuron0xc991190()*0.410684);
}

double testUp_8var::synapse0xca11bd0() {
   return (neuron0xc9914a0()*-0.21595);
}

double testUp_8var::synapse0xca11c10() {
   return (neuron0xc9917e0()*0.356553);
}

double testUp_8var::synapse0x10105fc0() {
   return (neuron0xca20410()*0.23501);
}

double testUp_8var::synapse0x10106000() {
   return (neuron0xca20720()*-0.44031);
}

double testUp_8var::synapse0xc9a4870() {
   return (neuron0xca20a30()*0.000366896);
}

double testUp_8var::synapse0xc9a48b0() {
   return (neuron0xca20d40()*0.512966);
}

double testUp_8var::synapse0xc9a48f0() {
   return (neuron0xca21050()*-0.0117369);
}

double testUp_8var::synapse0xc9a4930() {
   return (neuron0xc991190()*0.494143);
}

double testUp_8var::synapse0xc9a4970() {
   return (neuron0xc9914a0()*0.398463);
}

double testUp_8var::synapse0xca1f310() {
   return (neuron0xc9917e0()*-0.0599023);
}

double testUp_8var::synapse0xca125a0() {
   return (neuron0xca20410()*-0.163892);
}

double testUp_8var::synapse0xca125e0() {
   return (neuron0xca20720()*0.326437);
}

double testUp_8var::synapse0xca12620() {
   return (neuron0xca20a30()*-0.015027);
}

double testUp_8var::synapse0xca12660() {
   return (neuron0xca20d40()*-0.48688);
}

double testUp_8var::synapse0xca126a0() {
   return (neuron0xca21050()*0.0981572);
}

double testUp_8var::synapse0xca126e0() {
   return (neuron0xc991190()*0.00448099);
}

double testUp_8var::synapse0xca12720() {
   return (neuron0xc9914a0()*0.0950754);
}

double testUp_8var::synapse0xca12760() {
   return (neuron0xc9917e0()*-0.14896);
}

double testUp_8var::synapse0xca12ae0() {
   return (neuron0xca20410()*0.013189);
}

double testUp_8var::synapse0xca12b20() {
   return (neuron0xca20720()*-0.103336);
}

double testUp_8var::synapse0xca12b60() {
   return (neuron0xca20a30()*0.230339);
}

double testUp_8var::synapse0xca12ba0() {
   return (neuron0xca20d40()*-0.0319027);
}

double testUp_8var::synapse0xca12be0() {
   return (neuron0xca21050()*-0.20446);
}

double testUp_8var::synapse0xca12c20() {
   return (neuron0xc991190()*-0.343628);
}

double testUp_8var::synapse0xca12c60() {
   return (neuron0xc9914a0()*0.155758);
}

double testUp_8var::synapse0xca12ca0() {
   return (neuron0xc9917e0()*-0.236613);
}

double testUp_8var::synapse0xca13020() {
   return (neuron0xca20410()*0.115706);
}

double testUp_8var::synapse0xca13060() {
   return (neuron0xca20720()*0.229905);
}

double testUp_8var::synapse0xca130a0() {
   return (neuron0xca20a30()*0.0384846);
}

double testUp_8var::synapse0xca130e0() {
   return (neuron0xca20d40()*0.161893);
}

double testUp_8var::synapse0xca13120() {
   return (neuron0xca21050()*0.210048);
}

double testUp_8var::synapse0xca13160() {
   return (neuron0xc991190()*0.101615);
}

double testUp_8var::synapse0xca131a0() {
   return (neuron0xc9914a0()*-0.29161);
}

double testUp_8var::synapse0xca131e0() {
   return (neuron0xc9917e0()*0.15847);
}

double testUp_8var::synapse0xca13560() {
   return (neuron0xca20410()*-0.225062);
}

double testUp_8var::synapse0xca135a0() {
   return (neuron0xca20720()*-0.167485);
}

double testUp_8var::synapse0xca135e0() {
   return (neuron0xca20a30()*-0.179065);
}

double testUp_8var::synapse0xca13620() {
   return (neuron0xca20d40()*-0.243777);
}

double testUp_8var::synapse0xca13660() {
   return (neuron0xca21050()*0.280087);
}

double testUp_8var::synapse0xca136a0() {
   return (neuron0xc991190()*-0.495266);
}

double testUp_8var::synapse0xca136e0() {
   return (neuron0xc9914a0()*-0.0650865);
}

double testUp_8var::synapse0xca13720() {
   return (neuron0xc9917e0()*0.288873);
}

double testUp_8var::synapse0xca13aa0() {
   return (neuron0xca20410()*0.368761);
}

double testUp_8var::synapse0xca13ae0() {
   return (neuron0xca20720()*0.15478);
}

double testUp_8var::synapse0xca13b20() {
   return (neuron0xca20a30()*0.0250133);
}

double testUp_8var::synapse0xca13b60() {
   return (neuron0xca20d40()*-0.337637);
}

double testUp_8var::synapse0xca13ba0() {
   return (neuron0xca21050()*-0.427704);
}

double testUp_8var::synapse0xca13be0() {
   return (neuron0xc991190()*0.235051);
}

double testUp_8var::synapse0xca13c20() {
   return (neuron0xc9914a0()*0.224364);
}

double testUp_8var::synapse0xca13c60() {
   return (neuron0xc9917e0()*-0.348137);
}

double testUp_8var::synapse0xca13fe0() {
   return (neuron0xca20410()*0.137724);
}

double testUp_8var::synapse0xca14020() {
   return (neuron0xca20720()*-0.343124);
}

double testUp_8var::synapse0xca14060() {
   return (neuron0xca20a30()*0.35907);
}

double testUp_8var::synapse0xca140a0() {
   return (neuron0xca20d40()*-0.0686533);
}

double testUp_8var::synapse0xca140e0() {
   return (neuron0xca21050()*0.280467);
}

double testUp_8var::synapse0xca14120() {
   return (neuron0xc991190()*0.358314);
}

double testUp_8var::synapse0xca14160() {
   return (neuron0xc9914a0()*0.0974945);
}

double testUp_8var::synapse0xca141a0() {
   return (neuron0xc9917e0()*0.0866062);
}

double testUp_8var::synapse0xca14520() {
   return (neuron0xc991b20()*0.0184669);
}

double testUp_8var::synapse0xca14560() {
   return (neuron0xca1d510()*-0.00262066);
}

double testUp_8var::synapse0xca145a0() {
   return (neuron0xca1da20()*0.124156);
}

double testUp_8var::synapse0xca145e0() {
   return (neuron0xca1de80()*0.055817);
}

double testUp_8var::synapse0xca14620() {
   return (neuron0xca1e390()*-0.217577);
}

double testUp_8var::synapse0xca14660() {
   return (neuron0xca1e8f0()*0.019634);
}

double testUp_8var::synapse0xca146a0() {
   return (neuron0xca1ee00()*0.027855);
}

double testUp_8var::synapse0xca146e0() {
   return (neuron0xca11740()*0.2653);
}

double testUp_8var::synapse0xca14720() {
   return (neuron0xca11c50()*0.0185044);
}

double testUp_8var::synapse0xca1e870() {
   return (neuron0xca1e6e0()*0.0526754);
}

double testUp_8var::synapse0xca1e8b0() {
   return (neuron0xca127a0()*0.310598);
}

double testUp_8var::synapse0xca1f3e0() {
   return (neuron0xca12ce0()*-0.0813102);
}

double testUp_8var::synapse0xca1f420() {
   return (neuron0xca13220()*0.0272276);
}

double testUp_8var::synapse0xca29160() {
   return (neuron0xca13760()*-0.0980805);
}

double testUp_8var::synapse0xca291a0() {
   return (neuron0xca13ca0()*0.0629244);
}

