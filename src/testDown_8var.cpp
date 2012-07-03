#include "testDown_8var.hpp"
#include <cmath>

double testDown_8var::Value(int index,double in0,double in1,double in2,double in3,double in4,double in5,double in6,double in7) {
   input0 = (in0 - -0.00216753)/1.63676;
   input1 = (in1 - 3.2558e-05)/1.76583;
   input2 = (in2 - 152.626)/146.98;
   input3 = (in3 - 2.87781)/1.15782;
   input4 = (in4 - 110.379)/125.195;
   input5 = (in5 - 2.79505)/1.36887;
   input6 = (in6 - 88.6684)/109.583;
   input7 = (in7 - 2.69394)/1.45485;
   switch(index) {
     case 0:
         return neuron0x62677b0();
     default:
         return 0.;
   }
}

double testDown_8var::Value(int index, double* input) {
   input0 = (input[0] - -0.00216753)/1.63676;
   input1 = (input[1] - 3.2558e-05)/1.76583;
   input2 = (input[2] - 152.626)/146.98;
   input3 = (input[3] - 2.87781)/1.15782;
   input4 = (input[4] - 110.379)/125.195;
   input5 = (input[5] - 2.79505)/1.36887;
   input6 = (input[6] - 88.6684)/109.583;
   input7 = (input[7] - 2.69394)/1.45485;
   switch(index) {
     case 0:
         return neuron0x62677b0();
     default:
         return 0.;
   }
}

double testDown_8var::neuron0x99bfac0() {
   return input0;
}

double testDown_8var::neuron0x99bfdd0() {
   return input1;
}

double testDown_8var::neuron0x627bf20() {
   return input2;
}

double testDown_8var::neuron0x627c230() {
   return input3;
}

double testDown_8var::neuron0x627c540() {
   return input4;
}

double testDown_8var::neuron0x627c880() {
   return input5;
}

double testDown_8var::neuron0x627cbc0() {
   return input6;
}

double testDown_8var::neuron0x97ea930() {
   return input7;
}

double testDown_8var::input0x97ead80() {
   double input = 0.169665;
   input += synapse0x1449cbe0();
   input += synapse0x9441a20();
   input += synapse0x97eb000();
   input += synapse0x97eb040();
   input += synapse0x97eb080();
   input += synapse0x97eb0c0();
   input += synapse0x97eb100();
   input += synapse0x97eb140();
   return input;
}

double testDown_8var::neuron0x97ead80() {
   double input = input0x97ead80();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double testDown_8var::input0x97eb180() {
   double input = -0.10427;
   input += synapse0x97eb490();
   input += synapse0x97eb4d0();
   input += synapse0x97eb510();
   input += synapse0x97eb550();
   input += synapse0x97eb590();
   input += synapse0x97eb5d0();
   input += synapse0x97eb610();
   input += synapse0x97eb650();
   return input;
}

double testDown_8var::neuron0x97eb180() {
   double input = input0x97eb180();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double testDown_8var::input0x97eb690() {
   double input = -0.570568;
   input += synapse0x62791d0();
   input += synapse0x99c00e0();
   input += synapse0x62f3080();
   input += synapse0x62f30c0();
   input += synapse0x6279320();
   input += synapse0x6279360();
   input += synapse0x62793a0();
   input += synapse0x62793e0();
   return input;
}

double testDown_8var::neuron0x97eb690() {
   double input = input0x97eb690();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double testDown_8var::input0x6279420() {
   double input = 0.215718;
   input += synapse0x6279730();
   input += synapse0x6279770();
   input += synapse0x62797b0();
   input += synapse0x62797f0();
   input += synapse0x6279830();
   input += synapse0x6279870();
   input += synapse0x62798b0();
   input += synapse0x62798f0();
   return input;
}

double testDown_8var::neuron0x6279420() {
   double input = input0x6279420();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double testDown_8var::input0x6279930() {
   double input = 0.340174;
   input += synapse0x6279c40();
   input += synapse0x1449ca60();
   input += synapse0x127f8340();
   input += synapse0x6279210();
   input += synapse0x6279250();
   input += synapse0x6279290();
   input += synapse0x62792d0();
   input += synapse0x1449cd70();
   return input;
}

double testDown_8var::neuron0x6279930() {
   double input = input0x6279930();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double testDown_8var::input0x6279e90() {
   double input = -0.0782271;
   input += synapse0x627a1a0();
   input += synapse0x627a1e0();
   input += synapse0x627a220();
   input += synapse0x627a260();
   input += synapse0x627a2a0();
   input += synapse0x627a2e0();
   input += synapse0x627a320();
   input += synapse0x627a360();
   return input;
}

double testDown_8var::neuron0x6279e90() {
   double input = input0x6279e90();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double testDown_8var::input0x627a3a0() {
   double input = -0.416538;
   input += synapse0x627a6b0();
   input += synapse0x627a6f0();
   input += synapse0x627a730();
   input += synapse0x627a770();
   input += synapse0x627a7b0();
   input += synapse0x627a7f0();
   input += synapse0x627a830();
   input += synapse0x627a870();
   return input;
}

double testDown_8var::neuron0x627a3a0() {
   double input = input0x627a3a0();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double testDown_8var::input0x627a8b0() {
   double input = -0.359877;
   input += synapse0x627abc0();
   input += synapse0x627ac00();
   input += synapse0x627ac40();
   input += synapse0x627ac80();
   input += synapse0x627acc0();
   input += synapse0x627ad00();
   input += synapse0x627ad40();
   input += synapse0x627ad80();
   return input;
}

double testDown_8var::neuron0x627a8b0() {
   double input = input0x627a8b0();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double testDown_8var::input0x627adc0() {
   double input = 0.357678;
   input += synapse0x62f2ed0();
   input += synapse0x62f2f10();
   input += synapse0x6279c80();
   input += synapse0x6279cc0();
   input += synapse0x6279d00();
   input += synapse0x6279d40();
   input += synapse0x6279d80();
   input += synapse0x6279dc0();
   return input;
}

double testDown_8var::neuron0x627adc0() {
   double input = input0x627adc0();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double testDown_8var::input0x62658c0() {
   double input = -0.252647;
   input += synapse0x6265b70();
   input += synapse0x6265bb0();
   input += synapse0x6265bf0();
   input += synapse0x6265c30();
   input += synapse0x6265c70();
   input += synapse0x6265cb0();
   input += synapse0x6265cf0();
   input += synapse0x6265d30();
   return input;
}

double testDown_8var::neuron0x62658c0() {
   double input = input0x62658c0();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double testDown_8var::input0x6265d70() {
   double input = -0.308998;
   input += synapse0x62660b0();
   input += synapse0x62660f0();
   input += synapse0x6266130();
   input += synapse0x6266170();
   input += synapse0x62661b0();
   input += synapse0x62661f0();
   input += synapse0x6266230();
   input += synapse0x6266270();
   return input;
}

double testDown_8var::neuron0x6265d70() {
   double input = input0x6265d70();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double testDown_8var::input0x62662b0() {
   double input = -0.332606;
   input += synapse0x62665f0();
   input += synapse0x6266630();
   input += synapse0x6266670();
   input += synapse0x62666b0();
   input += synapse0x62666f0();
   input += synapse0x6266730();
   input += synapse0x6266770();
   input += synapse0x62667b0();
   return input;
}

double testDown_8var::neuron0x62662b0() {
   double input = input0x62662b0();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double testDown_8var::input0x62667f0() {
   double input = 0.0127758;
   input += synapse0x6266b30();
   input += synapse0x6266b70();
   input += synapse0x6266bb0();
   input += synapse0x6266bf0();
   input += synapse0x6266c30();
   input += synapse0x6266c70();
   input += synapse0x6266cb0();
   input += synapse0x6266cf0();
   return input;
}

double testDown_8var::neuron0x62667f0() {
   double input = input0x62667f0();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double testDown_8var::input0x6266d30() {
   double input = 0.439498;
   input += synapse0x6267070();
   input += synapse0x62670b0();
   input += synapse0x62670f0();
   input += synapse0x6267130();
   input += synapse0x6267170();
   input += synapse0x62671b0();
   input += synapse0x62671f0();
   input += synapse0x6267230();
   return input;
}

double testDown_8var::neuron0x6266d30() {
   double input = input0x6266d30();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double testDown_8var::input0x6267270() {
   double input = 0.385622;
   input += synapse0x62675b0();
   input += synapse0x62675f0();
   input += synapse0x6267630();
   input += synapse0x6267670();
   input += synapse0x62676b0();
   input += synapse0x62676f0();
   input += synapse0x6267730();
   input += synapse0x6267770();
   return input;
}

double testDown_8var::neuron0x6267270() {
   double input = input0x6267270();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double testDown_8var::input0x62677b0() {
   double input = 0.341397;
   input += synapse0x6267af0();
   input += synapse0x6267b30();
   input += synapse0x6267b70();
   input += synapse0x6267bb0();
   input += synapse0x6267bf0();
   input += synapse0x6267c30();
   input += synapse0x6267c70();
   input += synapse0x6267cb0();
   input += synapse0x6267cf0();
   input += synapse0x97eac40();
   input += synapse0x97eac80();
   input += synapse0x97eacc0();
   input += synapse0x97ead00();
   input += synapse0x97ead40();
   input += synapse0x62654b0();
   return input;
}

double testDown_8var::neuron0x62677b0() {
   double input = input0x62677b0();
   return (input * 1)+0;
}

double testDown_8var::synapse0x1449cbe0() {
   return (neuron0x99bfac0()*0.389373);
}

double testDown_8var::synapse0x9441a20() {
   return (neuron0x99bfdd0()*-0.48957);
}

double testDown_8var::synapse0x97eb000() {
   return (neuron0x627bf20()*-0.347346);
}

double testDown_8var::synapse0x97eb040() {
   return (neuron0x627c230()*-0.253495);
}

double testDown_8var::synapse0x97eb080() {
   return (neuron0x627c540()*0.386313);
}

double testDown_8var::synapse0x97eb0c0() {
   return (neuron0x627c880()*0.20958);
}

double testDown_8var::synapse0x97eb100() {
   return (neuron0x627cbc0()*0.446689);
}

double testDown_8var::synapse0x97eb140() {
   return (neuron0x97ea930()*0.217763);
}

double testDown_8var::synapse0x97eb490() {
   return (neuron0x99bfac0()*-0.213524);
}

double testDown_8var::synapse0x97eb4d0() {
   return (neuron0x99bfdd0()*-0.0276012);
}

double testDown_8var::synapse0x97eb510() {
   return (neuron0x627bf20()*-0.123119);
}

double testDown_8var::synapse0x97eb550() {
   return (neuron0x627c230()*0.493958);
}

double testDown_8var::synapse0x97eb590() {
   return (neuron0x627c540()*0.0578306);
}

double testDown_8var::synapse0x97eb5d0() {
   return (neuron0x627c880()*-0.216421);
}

double testDown_8var::synapse0x97eb610() {
   return (neuron0x627cbc0()*-0.345164);
}

double testDown_8var::synapse0x97eb650() {
   return (neuron0x97ea930()*-0.339928);
}

double testDown_8var::synapse0x62791d0() {
   return (neuron0x99bfac0()*-0.33531);
}

double testDown_8var::synapse0x99c00e0() {
   return (neuron0x99bfdd0()*0.0753224);
}

double testDown_8var::synapse0x62f3080() {
   return (neuron0x627bf20()*0.00461538);
}

double testDown_8var::synapse0x62f30c0() {
   return (neuron0x627c230()*0.398061);
}

double testDown_8var::synapse0x6279320() {
   return (neuron0x627c540()*0.35078);
}

double testDown_8var::synapse0x6279360() {
   return (neuron0x627c880()*0.19307);
}

double testDown_8var::synapse0x62793a0() {
   return (neuron0x627cbc0()*-0.376767);
}

double testDown_8var::synapse0x62793e0() {
   return (neuron0x97ea930()*0.269392);
}

double testDown_8var::synapse0x6279730() {
   return (neuron0x99bfac0()*-0.14597);
}

double testDown_8var::synapse0x6279770() {
   return (neuron0x99bfdd0()*0.149361);
}

double testDown_8var::synapse0x62797b0() {
   return (neuron0x627bf20()*0.500791);
}

double testDown_8var::synapse0x62797f0() {
   return (neuron0x627c230()*-0.0353816);
}

double testDown_8var::synapse0x6279830() {
   return (neuron0x627c540()*0.622072);
}

double testDown_8var::synapse0x6279870() {
   return (neuron0x627c880()*-0.107594);
}

double testDown_8var::synapse0x62798b0() {
   return (neuron0x627cbc0()*0.0013542);
}

double testDown_8var::synapse0x62798f0() {
   return (neuron0x97ea930()*0.000698204);
}

double testDown_8var::synapse0x6279c40() {
   return (neuron0x99bfac0()*-0.342214);
}

double testDown_8var::synapse0x1449ca60() {
   return (neuron0x99bfdd0()*-0.00185553);
}

double testDown_8var::synapse0x127f8340() {
   return (neuron0x627bf20()*-0.03626);
}

double testDown_8var::synapse0x6279210() {
   return (neuron0x627c230()*-0.219365);
}

double testDown_8var::synapse0x6279250() {
   return (neuron0x627c540()*0.348899);
}

double testDown_8var::synapse0x6279290() {
   return (neuron0x627c880()*-0.0959341);
}

double testDown_8var::synapse0x62792d0() {
   return (neuron0x627cbc0()*-0.0356587);
}

double testDown_8var::synapse0x1449cd70() {
   return (neuron0x97ea930()*-0.0658842);
}

double testDown_8var::synapse0x627a1a0() {
   return (neuron0x99bfac0()*0.359046);
}

double testDown_8var::synapse0x627a1e0() {
   return (neuron0x99bfdd0()*0.331923);
}

double testDown_8var::synapse0x627a220() {
   return (neuron0x627bf20()*0.339166);
}

double testDown_8var::synapse0x627a260() {
   return (neuron0x627c230()*0.238351);
}

double testDown_8var::synapse0x627a2a0() {
   return (neuron0x627c540()*-0.351606);
}

double testDown_8var::synapse0x627a2e0() {
   return (neuron0x627c880()*0.118595);
}

double testDown_8var::synapse0x627a320() {
   return (neuron0x627cbc0()*0.0162918);
}

double testDown_8var::synapse0x627a360() {
   return (neuron0x97ea930()*0.346751);
}

double testDown_8var::synapse0x627a6b0() {
   return (neuron0x99bfac0()*0.158194);
}

double testDown_8var::synapse0x627a6f0() {
   return (neuron0x99bfdd0()*0.365581);
}

double testDown_8var::synapse0x627a730() {
   return (neuron0x627bf20()*0.363621);
}

double testDown_8var::synapse0x627a770() {
   return (neuron0x627c230()*0.00775549);
}

double testDown_8var::synapse0x627a7b0() {
   return (neuron0x627c540()*0.214232);
}

double testDown_8var::synapse0x627a7f0() {
   return (neuron0x627c880()*-0.281382);
}

double testDown_8var::synapse0x627a830() {
   return (neuron0x627cbc0()*0.0867076);
}

double testDown_8var::synapse0x627a870() {
   return (neuron0x97ea930()*0.405421);
}

double testDown_8var::synapse0x627abc0() {
   return (neuron0x99bfac0()*-0.0762182);
}

double testDown_8var::synapse0x627ac00() {
   return (neuron0x99bfdd0()*0.149686);
}

double testDown_8var::synapse0x627ac40() {
   return (neuron0x627bf20()*-0.332379);
}

double testDown_8var::synapse0x627ac80() {
   return (neuron0x627c230()*0.110111);
}

double testDown_8var::synapse0x627acc0() {
   return (neuron0x627c540()*-0.0811025);
}

double testDown_8var::synapse0x627ad00() {
   return (neuron0x627c880()*-0.087762);
}

double testDown_8var::synapse0x627ad40() {
   return (neuron0x627cbc0()*0.119123);
}

double testDown_8var::synapse0x627ad80() {
   return (neuron0x97ea930()*-0.0480286);
}

double testDown_8var::synapse0x62f2ed0() {
   return (neuron0x99bfac0()*-0.154619);
}

double testDown_8var::synapse0x62f2f10() {
   return (neuron0x99bfdd0()*-0.0536008);
}

double testDown_8var::synapse0x6279c80() {
   return (neuron0x627bf20()*0.135606);
}

double testDown_8var::synapse0x6279cc0() {
   return (neuron0x627c230()*0.459735);
}

double testDown_8var::synapse0x6279d00() {
   return (neuron0x627c540()*-0.26862);
}

double testDown_8var::synapse0x6279d40() {
   return (neuron0x627c880()*-0.214857);
}

double testDown_8var::synapse0x6279d80() {
   return (neuron0x627cbc0()*0.34867);
}

double testDown_8var::synapse0x6279dc0() {
   return (neuron0x97ea930()*-0.431706);
}

double testDown_8var::synapse0x6265b70() {
   return (neuron0x99bfac0()*0.393696);
}

double testDown_8var::synapse0x6265bb0() {
   return (neuron0x99bfdd0()*-0.204502);
}

double testDown_8var::synapse0x6265bf0() {
   return (neuron0x627bf20()*0.116519);
}

double testDown_8var::synapse0x6265c30() {
   return (neuron0x627c230()*0.0163439);
}

double testDown_8var::synapse0x6265c70() {
   return (neuron0x627c540()*-0.130644);
}

double testDown_8var::synapse0x6265cb0() {
   return (neuron0x627c880()*-0.191666);
}

double testDown_8var::synapse0x6265cf0() {
   return (neuron0x627cbc0()*0.285123);
}

double testDown_8var::synapse0x6265d30() {
   return (neuron0x97ea930()*-0.0618732);
}

double testDown_8var::synapse0x62660b0() {
   return (neuron0x99bfac0()*-0.136011);
}

double testDown_8var::synapse0x62660f0() {
   return (neuron0x99bfdd0()*0.0208816);
}

double testDown_8var::synapse0x6266130() {
   return (neuron0x627bf20()*0.28422);
}

double testDown_8var::synapse0x6266170() {
   return (neuron0x627c230()*0.192218);
}

double testDown_8var::synapse0x62661b0() {
   return (neuron0x627c540()*-0.0724817);
}

double testDown_8var::synapse0x62661f0() {
   return (neuron0x627c880()*-0.303366);
}

double testDown_8var::synapse0x6266230() {
   return (neuron0x627cbc0()*0.349301);
}

double testDown_8var::synapse0x6266270() {
   return (neuron0x97ea930()*0.116522);
}

double testDown_8var::synapse0x62665f0() {
   return (neuron0x99bfac0()*-0.403199);
}

double testDown_8var::synapse0x6266630() {
   return (neuron0x99bfdd0()*0.136147);
}

double testDown_8var::synapse0x6266670() {
   return (neuron0x627bf20()*-0.160282);
}

double testDown_8var::synapse0x62666b0() {
   return (neuron0x627c230()*0.0646685);
}

double testDown_8var::synapse0x62666f0() {
   return (neuron0x627c540()*-0.490997);
}

double testDown_8var::synapse0x6266730() {
   return (neuron0x627c880()*0.0632854);
}

double testDown_8var::synapse0x6266770() {
   return (neuron0x627cbc0()*0.41705);
}

double testDown_8var::synapse0x62667b0() {
   return (neuron0x97ea930()*-0.366285);
}

double testDown_8var::synapse0x6266b30() {
   return (neuron0x99bfac0()*0.0131497);
}

double testDown_8var::synapse0x6266b70() {
   return (neuron0x99bfdd0()*-0.215072);
}

double testDown_8var::synapse0x6266bb0() {
   return (neuron0x627bf20()*0.136329);
}

double testDown_8var::synapse0x6266bf0() {
   return (neuron0x627c230()*-0.299338);
}

double testDown_8var::synapse0x6266c30() {
   return (neuron0x627c540()*0.0422075);
}

double testDown_8var::synapse0x6266c70() {
   return (neuron0x627c880()*-0.0937092);
}

double testDown_8var::synapse0x6266cb0() {
   return (neuron0x627cbc0()*-0.041287);
}

double testDown_8var::synapse0x6266cf0() {
   return (neuron0x97ea930()*-0.161816);
}

double testDown_8var::synapse0x6267070() {
   return (neuron0x99bfac0()*-0.108866);
}

double testDown_8var::synapse0x62670b0() {
   return (neuron0x99bfdd0()*0.0819662);
}

double testDown_8var::synapse0x62670f0() {
   return (neuron0x627bf20()*0.0930976);
}

double testDown_8var::synapse0x6267130() {
   return (neuron0x627c230()*0.454644);
}

double testDown_8var::synapse0x6267170() {
   return (neuron0x627c540()*0.328714);
}

double testDown_8var::synapse0x62671b0() {
   return (neuron0x627c880()*0.147577);
}

double testDown_8var::synapse0x62671f0() {
   return (neuron0x627cbc0()*-0.382583);
}

double testDown_8var::synapse0x6267230() {
   return (neuron0x97ea930()*-0.0862269);
}

double testDown_8var::synapse0x62675b0() {
   return (neuron0x99bfac0()*-0.543976);
}

double testDown_8var::synapse0x62675f0() {
   return (neuron0x99bfdd0()*0.364287);
}

double testDown_8var::synapse0x6267630() {
   return (neuron0x627bf20()*0.475494);
}

double testDown_8var::synapse0x6267670() {
   return (neuron0x627c230()*0.276835);
}

double testDown_8var::synapse0x62676b0() {
   return (neuron0x627c540()*0.443653);
}

double testDown_8var::synapse0x62676f0() {
   return (neuron0x627c880()*0.378547);
}

double testDown_8var::synapse0x6267730() {
   return (neuron0x627cbc0()*0.218047);
}

double testDown_8var::synapse0x6267770() {
   return (neuron0x97ea930()*-0.272174);
}

double testDown_8var::synapse0x6267af0() {
   return (neuron0x97ead80()*-0.0176205);
}

double testDown_8var::synapse0x6267b30() {
   return (neuron0x97eb180()*-0.17456);
}

double testDown_8var::synapse0x6267b70() {
   return (neuron0x97eb690()*-0.172125);
}

double testDown_8var::synapse0x6267bb0() {
   return (neuron0x6279420()*-0.148405);
}

double testDown_8var::synapse0x6267bf0() {
   return (neuron0x6279930()*0.193311);
}

double testDown_8var::synapse0x6267c30() {
   return (neuron0x6279e90()*-0.0386847);
}

double testDown_8var::synapse0x6267c70() {
   return (neuron0x627a3a0()*0.110551);
}

double testDown_8var::synapse0x6267cb0() {
   return (neuron0x627a8b0()*0.0292452);
}

double testDown_8var::synapse0x6267cf0() {
   return (neuron0x627adc0()*0.0913651);
}

double testDown_8var::synapse0x97eac40() {
   return (neuron0x62658c0()*0.0022672);
}

double testDown_8var::synapse0x97eac80() {
   return (neuron0x6265d70()*-0.0256277);
}

double testDown_8var::synapse0x97eacc0() {
   return (neuron0x62662b0()*0.0173415);
}

double testDown_8var::synapse0x97ead00() {
   return (neuron0x62667f0()*0.0609495);
}

double testDown_8var::synapse0x97ead40() {
   return (neuron0x6266d30()*0.292384);
}

double testDown_8var::synapse0x62654b0() {
   return (neuron0x6267270()*-0.0383951);
}

