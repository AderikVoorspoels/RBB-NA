/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2020 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "Colvar.h"
#include "ActionRegister.h"
#include "tools/Pbc.h"
#include "tools/Vector.h"
#include "tools/Keywords.h"

#include <string>
#include <cmath>
//#include <unordered_map>
#include <vector>

using namespace std;

namespace PLMD {
namespace colvar {

//+PLUMEDOC Tsukuba
/*
This file provides a template for if you want to introduce a new CV.

<!-----You should add a description of your CV here---->
This CV will calculate the Tsukuba convention inter bp rotations(tilt, roll, twist)

\par Examples

<!---You should put an example of how to use your CV here--->error: expected unqualified-id before ‘switch’

\plumedfile
# This should be a sample input.
rotations: TSUKUBAINTERROTFITTED BASE1=1,2,3 BASE2=4,5,6 BASE3=7,8,9 BASE4=10,11,12
PRINT ARG=rotations.roll,rotations.tilt,rotations.twist STRIDE=100 FILE=COLVAR
#here if the order of atoms per base is N1(Y)/N9(R) C'1 C2(Y)/C4(R) where Y are the Pyrimidenes (T and C), R are the Purines (A and G)
\endplumedfile
<!---You should reference here the other actions used in this example--->
(see also \ref PRINT)

*/
//+ENDPLUMEDOC

class cirtoNA : public Colvar {
  bool components;
  bool pbc;
  bool RNA;
  bool FITTED;
  enum mode {pair=0, step=1, sequence=2, empty=3};
  mode INPUTMODE;
  enum BaseType {A = 'A', C='C', G='G', T='T', U='U'};
  vector<char> BASES;
  
public:
  explicit cirtoNA(const ActionOptions&);
// active methods:
  void calculate() override;
  static void registerKeywords(Keywords& keys);
//  static std::unordered_map<std::string,BaseType> const BaseTable;

private:
  void completeBase(char type, AtomNumber start, vector<AtomNumber> *atoms, AtomNumber *next);
  
  void fit(char type, vector<VectorGeneric<3>> atoms, VectorGeneric<3> *fit1, VectorGeneric<3> *fit2, VectorGeneric<3> *fit3);
  double Determinant3x3(vector< VectorGeneric<3> > Matrix);
  VectorGeneric<3> EigenVals(vector<VectorGeneric<3>> Matrix);
  void EigenValsAndVectorsSymmetric(vector<VectorGeneric<3>> Matrix, VectorGeneric<3>* EigenValues, vector<VectorGeneric<3>>* EigenVectors);
  
  double Norm(VectorGeneric<3> vec);
  VectorGeneric<3> Normalize(VectorGeneric<3> vec);
  VectorGeneric<3> AxisRotation(VectorGeneric<3> v, VectorGeneric<3> k, double theta);

  vector<VectorGeneric<3>> Skew(VectorGeneric<3> phi);

  vector<VectorGeneric<3>> RotationMatrix(vector<VectorGeneric<3>> frame1, vector<VectorGeneric<3>> frame2);
  VectorGeneric<3> RotationAxis(vector<VectorGeneric<3>> Q);
  double RotationAngle(vector<VectorGeneric<3>> Q);
  VectorGeneric<3> RotationVector(vector<VectorGeneric<3>> Q);
  vector<VectorGeneric<3>> AxisAngleToMatrix(VectorGeneric<3> RotationVector);

  vector<VectorGeneric<3>> midFrame(vector<VectorGeneric<3>> frame1, VectorGeneric<3> RotationVector);

  vector<VectorGeneric<3>> BaseFrame(VectorGeneric<3> glycosidic, VectorGeneric<3> internal);
  VectorGeneric<3> RefPoint(VectorGeneric<3> glycosidic, VectorGeneric<3> normal, VectorGeneric<3> nitrogen);
  
  vector<VectorGeneric<3>> DNormalizeDvec(VectorGeneric<3> normedVec, VectorGeneric<3> vec);
  VectorGeneric<3> DNormDvec(double norm, VectorGeneric<3> vec);

  vector<VectorGeneric<3>> DRotationAngleDQ(double theta, vector<VectorGeneric<3>> Q);
  vector<vector<VectorGeneric<3>>> DRotationAxisDQ(VectorGeneric<3> Axis, vector<VectorGeneric<3>> Q);
  vector<vector<VectorGeneric<3>>> DVectorDTransformation(VectorGeneric<3> RotationVector, vector<VectorGeneric<3>> Q);
  vector<vector<VectorGeneric<3>>> DTransformationDRot(vector<VectorGeneric<3>> Q, vector<VectorGeneric<3>> frame, double direction);
  vector<vector<VectorGeneric<3>>> DMidDRot(vector<VectorGeneric<3>>  midFrame, vector<VectorGeneric<3>> frame);
};

PLUMED_REGISTER_ACTION(cirtoNA,"CIRTONA")

void cirtoNA::registerKeywords(Keywords& keys) {
  Colvar::registerKeywords( keys );
  keys.addFlag("FITTED",false,"fit the bases with an ideal base before calulation of relative deformation");
  keys.addFlag("RNA",false,"if checking RNA set to true complete A in sequence with coplementary U");
  
  keys.add("atoms","BASE1","the atoms setting the refference frame of the first base");
  keys.add("optional","TYPE1", "the type of the first base (A, C, G, T or U)");
  keys.add("atoms","BASE2","the atoms setting the refference frame of the second base (opposite from first base)");
  keys.add("optional","TYPE2", "the type of the first base (A, C, G, T or U)");
  keys.add("atoms","BASE3","the atoms setting the refference frame of the third base (same strand as first base)");
  keys.add("optional","TYPE3", "the type of the first base (A, C, G, T or U)");
  keys.add("atoms","BASE4","the atoms setting the refference frame of the fourth base (opposite from third base)");
  keys.add("optional","TYPE4", "the type of the first base (A, C, G, T or U)");
  keys.add("atoms","ENDS","the four atoms c'1 delimiting the strand of DNA/RNA that is to be analyzed (first the 5'-3' then 3'-5')");
  keys.add("optional","SEQUENCE","the squence of the strand to be computed on");

 // keys.addFlag("COMPONENTSSTEPROT",true,"calculate the tilt, roll and twist store as label.tilt, label.roll and label.twist, will be set to true always"); //Output the rotational deformations in the bpstep
 // keys.addFlag("COMPONENTSSTEPTRANS",true,"calculate the tilt, roll and twist store as label.tilt, label.roll and label.twist, will be set to true always"); //Output the translational deformations in the bpstep
 // keys.addFlag("COMPONENTSPAIRROT",true,"calculate the tilt, roll and twist store as label.tilt, label.roll and label.twist, will be set to true always"); //Output the rotational deformations in the bp
 // keys.addFlag("COMPONENTSPAIRTRANS",true,"calculate the tilt, roll and twist store as label.tilt, label.roll and label.twist, will be set to true always"); //Output the translational deformations in the bp
  keys.addOutputComponent("tilt", "COMPONENTSSTEPROT", "the tilt component of the rotation");
  keys.addOutputComponent("roll", "COMPONENTSSTEPROT", "the roll component of the rotation");
  keys.addOutputComponent("twist", "COMPONENTSSTEPROT", "the twist component of the rotation");
  
  keys.addOutputComponent("shift", "COMPONENTSSTEPTRANS", "the tilt component of the rotation");
  keys.addOutputComponent("slide", "COMPONENTSSTEPTRANS", "the roll component of the rotation");
  keys.addOutputComponent("rise", "COMPONENTSSTEPTRANS", "the twist component of the rotation");
  
  keys.addOutputComponent("buckle", "COMPONENTSPAIRROT", "the tilt component of the rotation");
  keys.addOutputComponent("propeller", "COMPONENTSPAIRROT", "the roll component of the rotation");
  keys.addOutputComponent("opening", "COMPONENTSPAIRROT", "the twist component of the rotation");
  
  keys.addOutputComponent("shear", "COMPONENTSPAIRTRANS", "the tilt component of the rotation");
  keys.addOutputComponent("stretch", "COMPONENTSPAIRTRANS", "the roll component of the rotation");
  keys.addOutputComponent("stagger", "COMPONENTSPAIRTRANS", "the twist component of the rotation");
  
  for(int i=0; i<50; i++){
    keys.addOutputComponent("tilt"+to_string(i), "SEQCOMPONENTSSTEPROT"+to_string(i), "the tilt component of the rotation of step "+to_string(i));
    keys.addOutputComponent("roll"+to_string(i), "SEQCOMPONENTSSTEPROT"+to_string(i), "the roll component of the rotation of step "+to_string(i));
    keys.addOutputComponent("twist"+to_string(i), "SEQCOMPONENTSSTEPROT"+to_string(i), "the twist component of the rotation of step "+to_string(i));
    keys.addOutputComponent("shift"+to_string(i), "SEQCOMPONENTSSTEPTRANS"+to_string(i), "the tilt component of the rotation of step "+to_string(i));
    keys.addOutputComponent("slide"+to_string(i), "SEQCOMPONENTSSTEPTRANS"+to_string(i), "the roll component of the rotation of step "+to_string(i));
    keys.addOutputComponent("rise"+to_string(i), "SEQCOMPONENTSSTEPTRANS"+to_string(i), "the twist component of the rotation of step "+to_string(i));
  }

  for(int i=0; i<51; i++){
    keys.addOutputComponent("buckle"+to_string(i), "SEQCOMPONENTSPAIRROT"+to_string(i), "the tilt component of the rotation of pair "+to_string(i));
    keys.addOutputComponent("propeller"+to_string(i), "SEQCOMPONENTSPAIRROT"+to_string(i), "the roll component of the rotation of pair "+to_string(i));
    keys.addOutputComponent("opening"+to_string(i), "SEQCOMPONENTSPAIRROT"+to_string(i), "the twist component of the rotation of pair "+to_string(i));
    keys.addOutputComponent("shear"+to_string(i), "SEQCOMPONENTSPAIRTRANS"+to_string(i), "the tilt component of the rotation of pair "+to_string(i));
    keys.addOutputComponent("stretch"+to_string(i), "SEQCOMPONENTSPAIRTRANS"+to_string(i), "the roll component of the rotation of pair "+to_string(i));
    keys.addOutputComponent("stagger"+to_string(i), "SEQCOMPONENTSPAIRTRANS"+to_string(i), "the twist component of the rotation of pair "+to_string(i));
  }
  
}


void cirtoNA::completeBase(char type, AtomNumber c1, vector<AtomNumber> *atoms, AtomNumber *next){
  switch(type){
    case('A'):{
      (*atoms).resize(10);
      int start=c1.index();
      AtomNumber a2 = PLMD::AtomNumber::index(start + 2); //N9   
      AtomNumber a3 = PLMD::AtomNumber::index(start + 15);  //C4
      AtomNumber a4 = PLMD::AtomNumber::index(start + 3);  //C8
      AtomNumber a5 = PLMD::AtomNumber::index(start + 11);  //N1
      AtomNumber a6 = PLMD::AtomNumber::index(start + 12);  //C2
      AtomNumber a7 = PLMD::AtomNumber::index(start + 14);  //N3
      AtomNumber a8 = PLMD::AtomNumber::index(start + 6);  //C5
      AtomNumber a9 = PLMD::AtomNumber::index(start + 7);  //C6
      AtomNumber a10 = PLMD::AtomNumber::index(start + 5);  //N7
      vector<AtomNumber> Atoms{c1, a2, a3, a4, a5, a6, a7, a8, a9, a10};
      *atoms=Atoms;
      *next=PLMD::AtomNumber::index(start+32);
      break;
    }
    case('C'):{
      (*atoms).resize(7);
      int start=c1.index();
      AtomNumber a2 = PLMD::AtomNumber::index(start + 2);  //N1
      AtomNumber a3 = PLMD::AtomNumber::index(start + 12);  //C2
      AtomNumber a4 = PLMD::AtomNumber::index(start + 3);  //C6
      AtomNumber a5 = PLMD::AtomNumber::index(start + 11);  //N3
      AtomNumber a6 = PLMD::AtomNumber::index(start + 7);  //C4
      AtomNumber a7 = PLMD::AtomNumber::index(start + 5);  //C5
      vector<AtomNumber> Atoms{c1, a2, a3, a4, a5, a6, a7};
      *atoms=Atoms;
      *next=PLMD::AtomNumber::index(start+30);
      break;
    }
    case('G'):{
      (*atoms).resize(10);
      int start=c1.index();
      AtomNumber a2 = PLMD::AtomNumber::index(start + 2);  //N9
      AtomNumber a3 = PLMD::AtomNumber::index(start + 16);  //C4
      AtomNumber a4 = PLMD::AtomNumber::index(start + 3);  //C8
      AtomNumber a5 = PLMD::AtomNumber::index(start + 9);  //N1
      AtomNumber a6 = PLMD::AtomNumber::index(start + 11);  //C2
      AtomNumber a7 = PLMD::AtomNumber::index(start + 15);  //N3
      AtomNumber a8 = PLMD::AtomNumber::index(start + 6);  //C5
      AtomNumber a9 = PLMD::AtomNumber::index(start + 7);  //C6
      AtomNumber a10 = PLMD::AtomNumber::index(start + 5);  //N7
      vector<AtomNumber> Atoms{c1, a2, a3, a4, a5, a6, a7, a8, a9, a10};
      *atoms=Atoms;
      *next=PLMD::AtomNumber::index(start+33);
      break;
    }
    case('T'):{
      (*atoms).resize(7);
      int start=c1.index();
      AtomNumber a2 = PLMD::AtomNumber::index(start + 2);  //N1
      AtomNumber a3 = PLMD::AtomNumber::index(start + 14);  //C2
      AtomNumber a4 = PLMD::AtomNumber::index(start + 3);  //C6
      AtomNumber a5 = PLMD::AtomNumber::index(start + 12);  //N3
      AtomNumber a6 = PLMD::AtomNumber::index(start + 10);  //C4
      AtomNumber a7 = PLMD::AtomNumber::index(start + 5);  //C5
      vector<AtomNumber> Atoms{c1, a2, a3, a4, a5, a6, a7};
      *atoms=Atoms;
      *next=PLMD::AtomNumber::index(start+32);
      break;
    }
    case('U'):{
      (*atoms).resize(7);
      int start=c1.index();
      AtomNumber a2 = PLMD::AtomNumber::index(start + 2);  //N1
      AtomNumber a3 = PLMD::AtomNumber::index(start + 11);  //C2
      AtomNumber a4 = PLMD::AtomNumber::index(start + 3);  //C6
      AtomNumber a5 = PLMD::AtomNumber::index(start + 9);  //N3
      AtomNumber a6 = PLMD::AtomNumber::index(start + 7);  //C4
      AtomNumber a7 = PLMD::AtomNumber::index(start + 5);  //C5
      vector<AtomNumber> Atoms{c1, a2, a3, a4, a5, a6, a7};
      *atoms=Atoms;
      *next=PLMD::AtomNumber::index(start+29);
      break;
    }
      return;
  }
}


cirtoNA::cirtoNA(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  pbc(true),
  components(true),
  RNA(false),
  FITTED(false)
{
  vector<AtomNumber> atoms,base1,base2,base3,base4,ends;
  AtomNumber i;
  string dummy;
  atoms.resize(0);

  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

  bool fitted;
  parseFlag("FITTED",fitted);
  parse("FITTED",dummy);
  FITTED = fitted;
  
  bool ribo;
  parseFlag("RNA", ribo);
  parse("RNA",dummy);
  RNA = ribo;

  // take the inputted atoms
  parseAtomList("BASE1",base1);
  parseAtomList("BASE2",base2);
  parseAtomList("BASE3",base3);
  parseAtomList("BASE4",base4);
  parseAtomList("ENDS",ends);
  // use context to see the mode of operation used
  
  INPUTMODE = mode::empty;
  if(ends.empty()){
    if(base1.empty() || base2.empty()){error("please specify base1 and base2 if not using sequence mode");} 
    if( (!(base3.empty() ))||(! (base4.empty())) ){
      INPUTMODE = mode::step;
      if(base3.empty() || base4.empty()){error("please specify both base3 and base4 if using step mode");}
    }else{
    INPUTMODE = mode::pair;
    }
  }else if((base1.empty() && base2.empty() && base3.empty() && base4.empty())){
    INPUTMODE = mode::sequence;
  }else{error("base and ends keywords are incompatible");}    



  //do input parsing based on the mode of operation
  switch(INPUTMODE){
    case(pair):{// request the required atoms as needed and create output handlinmatplotlib close figure
      BASES.resize(2);
      vector<AtomNumber> atomsBase1 = base1;
      vector<AtomNumber> atomsBase2 = base2;
      if(FITTED){
        if(!(base1.size()==1)){error("number of atoms in base1 should be 1 (the carbon involved in the glycosidic bond C'1) if you want to specify it for fitting.");}
        if(!(base2.size()==1)){error("number of atoms in base2 should be 1 (the carbon involved in the glycosidic bond C'1) if you want to specify it for fitting.");}
        string TypeIn;
        parse("TYPE1",TypeIn);
        BASES[0] = TypeIn[0];
        completeBase(BASES[0], base1[0], &atomsBase1, &i);
        
        parse("TYPE2",TypeIn);
        BASES[1] = TypeIn[0];
        completeBase(BASES[1], base2[0], &atomsBase2, &i);
      }else{
        if(!(base1.size()==3)){error("number of atoms in base1 should be 3 (C'1, N1/N9 ,C2/C4) if you want to specify it without fitting.");}
        if(!(base2.size()==3)){error("number of atoms in base2 should be 3 (C'1, N1/N9 ,C2/C4) if you want to specify it without fitting.");}
      }

      atoms.resize(atomsBase1.size()+atomsBase2.size());
      for(int j = 0; j < atomsBase1.size(); j++){
        atoms[j] = atomsBase1[j];
      }
      for(int j = 0; j < atomsBase2.size(); j++){
        atoms[j+atomsBase1.size()] = atomsBase2[j];
      }

      addComponentWithDerivatives("buckle"); componentIsNotPeriodic("buckle");
      addComponentWithDerivatives("propeller"); componentIsNotPeriodic("propeller");
      addComponentWithDerivatives("opening"); componentIsNotPeriodic("opening");
      
      addComponentWithDerivatives("shear"); componentIsNotPeriodic("shear");
      addComponentWithDerivatives("stretch"); componentIsNotPeriodic("stretch");
      addComponentWithDerivatives("stagger"); componentIsNotPeriodic("stagger");
      
      log.printf("Tsukuba Basepair between base1 %d and base2 %d \n", atomsBase1[0].serial(),atomsBase2[0].serial());
      break;
    }    
    case(step):{// request the required atoms as needed and create output handling
      BASES.resize(4);
      int b, j, N;
      vector<vector<PLMD::AtomNumber>> atomsBases = {base1, base2, base3, base4};
      string TypeIn;
      AtomNumber c1;

      if(FITTED) {
          for (b = 0; b < 4; b++) {
              if (!(atomsBases[b].size() == 1)) {
                  error("number of atoms in base" + to_string(b + 1) +
                        "should be 1 (the carbon involved in the glycosidic bond C'1) if you want to specify for fitting.");
              }

              parse("TYPE" + to_string(b+1) ,TypeIn);
              BASES[b] = TypeIn[0];
              c1 = atomsBases[b][0];
              completeBase(BASES[b], c1, &atomsBases[b], &i);

          }
      }else{
          for (b = 0; b < 4; b++) {
              if (!(atomsBases[b].size() == 1)) {
                  error("number of atoms in base" + to_string(b + 1) +
                        " should be 3 (C'1, N1/N9 ,C2/C4) if you want to specify without fitting.");
              }
          }
      }

      N = atomsBases[0].size()+atomsBases[1].size()+atomsBases[2].size()+atomsBases[3].size();
      atoms.resize(N);
      printf("N: %d = %d + %d + %d + %d \n", N, atomsBases[0].size(), atomsBases[1].size(), atomsBases[2].size(), atomsBases[3].size());
      N=0;

      for(b = 0; b < 4; b++) {
          for (int j = 0; j < atomsBases[b].size(); j++) {
              atoms[N] = atomsBases[b][j];
              N++;
              printf("%d ", atomsBases[b][j]);
          }
          printf("\n");
      }
      
      addComponentWithDerivatives("tilt"); componentIsNotPeriodic("tilt");
      addComponentWithDerivatives("roll"); componentIsNotPeriodic("roll");
      addComponentWithDerivatives("twist"); componentIsNotPeriodic("twist");
      
      addComponentWithDerivatives("shift"); componentIsNotPeriodic("shift");
      addComponentWithDerivatives("slide"); componentIsNotPeriodic("slide");
      addComponentWithDerivatives("rise"); componentIsNotPeriodic("rise");

      log.printf("Tsukuba Basepairstep between base1 %d, base2 %d, base3 %d and base4 %d\n", atomsBases[0][0].serial(), atomsBases[1][0].serial(), atomsBases[2][0].serial(), atomsBases[3][0].serial());
      break;

    }
    case(sequence):{// do a check on the sequence, request the required atoms as needed, and create output handling
      if(!(ends.size()==4)){
        error("please specify the four ends of the double stranded nucleic acid");
      }
      string seqIn;
      parse("SEQUENCE",seqIn);
      int j, N, n;
      int length = seqIn.size();
      BASES.resize(2*length);
      vector<AtomNumber> Base53;
      vector<vector<PLMD::AtomNumber>> Bases53;
      Bases53.resize(length);
      vector<AtomNumber> Base35;
      vector<vector<PLMD::AtomNumber>> Bases35;
      Bases35.resize(length);
      PLMD::AtomNumber current53 = ends[0];
      PLMD::AtomNumber current35 = ends[2];
      PLMD::AtomNumber next53;
      PLMD::AtomNumber next35;
      char type53;
      char type35;
      
      N = 0;
      for(j=0; j<length; j++){// loop to compose the complementary sequence
        type53 = seqIn[j];
        
        switch(type53){
          case('A'): {
              if (RNA) { type35 = 'U'; } else { type35 = 'T'; }
              break;
          }
          case('C'): {
              type35 = 'G';
              break;
          }
          case('G'): {
              type35 = 'C';
              break;
          }
          case('T'): {
              type35 = 'A';
              break;
          }
          case('U'): {
              type35 = 'A';
              break;
          }
        }
       
        BASES[2*j] = type53;
        BASES[2*j+1] = type35;

      }

      for(j=0; j<length; j++){// loop to get the atoms

        completeBase(BASES[2*j], current53, &Base53, &next53);
        if(RNA){next53 = PLMD::AtomNumber::index(next53.index() + 1);}
        current53 = next53;
        Bases53[j] = Base53;

        completeBase(BASES[2*(length-j-1)+1], current35, &Base35, &next35);// in the complementary strand you start at the end of the sequence
        if(RNA){next35 = PLMD::AtomNumber::index(next35.index() + 1);}
        current35 = next35;
        Bases35[length-j-1] = Base35;

        atoms.resize(atoms.size()+Bases53[j].size()+Bases35[length-j-1].size());
      }

      N=0; // set the starting index in the atom list
      for(j=0; j<length; j++){// loop putting the atoms in the order conventional in the two above modes (alternating bases between strands to make complementary bases near each other)
        
        for(n=0; n<Bases53[j].size(); n++){
          atoms[N] = Bases53[j][n]; // add the atoms of the base in the 5->3 strand to the bottom of the atom list
          N++;
        }

        for(n=0; n<Bases35[j].size(); n++){
          atoms[N] = Bases35[j][n]; // add the atoms of the base in the 3->5 strand to the bottom of the atom list
          N++;
        }

        // add length bp parameters as components
        addComponentWithDerivatives("buckle"+to_string(j)); componentIsNotPeriodic("buckle"+to_string(j));
        addComponentWithDerivatives("propeller"+to_string(j)); componentIsNotPeriodic("propeller"+to_string(j));
        addComponentWithDerivatives("opening"+to_string(j)); componentIsNotPeriodic("opening"+to_string(j));
         
        addComponentWithDerivatives("shear"+to_string(j)); componentIsNotPeriodic("shear"+to_string(j));
        addComponentWithDerivatives("stretch"+to_string(j)); componentIsNotPeriodic("stretch"+to_string(j));
        addComponentWithDerivatives("stagger"+to_string(j)); componentIsNotPeriodic("stagger"+to_string(j));
      
        if(!(j==0)){  // add length-1 bp step parameters as components
          addComponentWithDerivatives("tilt"+to_string(j-1)); componentIsNotPeriodic("tilt"+to_string(j-1));
          addComponentWithDerivatives("roll"+to_string(j-1)); componentIsNotPeriodic("roll"+to_string(j-1));
          addComponentWithDerivatives("twist"+to_string(j-1)); componentIsNotPeriodic("twist"+to_string(j-1));
      
          addComponentWithDerivatives("shift"+to_string(j-1)); componentIsNotPeriodic("shift"+to_string(j-1));
          addComponentWithDerivatives("slide"+to_string(j-1)); componentIsNotPeriodic("slide"+to_string(j-1));
          addComponentWithDerivatives("rise"+to_string(j-1)); componentIsNotPeriodic("rise"+to_string(j-1));
        }
      }     
      log.printf("cirtoNA between ends %d, %d, %d and %d\n", ends[0].serial(),ends[1].serial(),ends[2].serial(),ends[3].serial());
      break;
    }
  }
  checkRead();

  // some log output
  if(pbc){ log.printf("  using periodic boundary conditions\n");
  }else{   log.printf("  without periodic boundary conditions\n");};

  
  if(RNA){ log.printf("  Nucleic acid is RNA\n");
  }else{   log.printf("  Nucleic acid is DNA\n");};
  
  
  if(FITTED){ log.printf("  using fitting to optimize parameters\n");
  }else{   log.printf("  without using fitting to optimize parameters\n");};


  log.printf("  %d Atoms will be used\n", atoms.size());

  log.printf("  version is 4.may.22\n");

  //get the atoms, then calculation can start

  requestAtoms(atoms);
}


//fitting
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void cirtoNA::fit(char type, vector<VectorGeneric<3>> atoms, VectorGeneric<3> *fit1, VectorGeneric<3> *fit2, VectorGeneric<3> *fit3) {
  vector<VectorGeneric<3> > modelBaseCoords, T, eigenVectorsH, eigenVectorsK, Tt, TxTt, TtxT, R;
  VectorGeneric<3> eigenValues, eigenValuesDeSquared, fitAtom1, fitAtom2, fitAtom3;
  VectorGeneric<3> MC = Vector(0, 0, 0);
  int natoms, l, n, i, j;
  double DetT;

  T.resize(3);
  Tt.resize(3);
  TtxT.resize(3);
  TxTt.resize(3);
  R.resize(3);
  eigenVectorsH.resize(3);
  eigenVectorsK.resize(3);
  if(type == 'A' or type == 'G'){
    natoms = 10;
  }else{natoms = 7;}


  for (n = 0; n < natoms; n++) {
    MC += atoms[n];
  }
  MC = MC/natoms;
  for (n = 0; n < natoms; n++) {
    atoms[n] = atoms[n] - MC;
  }


     // define the model that the fit is to
  switch(type) {
    case ('A'): {
      modelBaseCoords = {Vector(1.57340, -2.41044, -0.12190), Vector(0.37786, -1.52945, -0.00531),
                        Vector(0.38110, -0.16006, -0.04011), Vector(-0.94428, -1.87750, 0.15888),
                        Vector(-0.11123, 2.43959, -0.04964), Vector(1.10298, 1.90641, -0.17768),
                        Vector(1.45731, 0.64079, -0.18694), Vector(-0.91754, 0.23706, 0.10110),
                        Vector(-1.16848, 1.61858, 0.09523), Vector(-1.75109, -0.86494, 0.22638)};
          // modelBaseNames = 
        break;
    }
    case('C'): {
      modelBaseCoords = {Vector(1.94866, -1.45161, -0.19029), Vector(0.74385, -0.58352, -0.07229),
                        Vector(0.93020, 0.79520, -0.12929), Vector(-0.49500, -1.12092, 0.08671),
                        Vector(-0.15547, 1.60399, -0.02329), Vector(-1.37969, 1.08626, 0.13371),
                        Vector(-1.59254, -0.32937, 0.19471)};
          // modelBaseNames =
      break;
    }
    case('G'): {
      modelBaseCoords = {Vector(1.58195, -2.39594, -0.12320), Vector(0.37714, -1.52785, -0.00520),
                        Vector(0.37551, -0.14975, -0.04020), Vector(-0.94239, -1.88689, 0.15880),
                        Vector(-0.10165, 2.41153, -0.05020), Vector(1.18367, 1.92662, -0.18720),
                        Vector(1.47916, 0.62851, -0.18920), Vector(-0.93480, 0.24071, 0.10280),
                        Vector(-1.25149, 1.62390, 0.10480), Vector(-1.76708, -0.87080, 0.22880)};
          // modelBaseNames =
      break;
    }
    case('T'): {
      modelBaseCoords = {Vector(1.95363, -1.45659, -0.19073), Vector(0.75787, -0.57587, -0.07419),
                        Vector(0.97215, 0.78019, -0.13405), Vector(-0.49400, -1.10811, 0.08633),
                        Vector(-0.15841, 1.56515, -0.02137), Vector(-1.45350, 1.11665, 0.14102),
                        Vector(-1.57773, -0.32139, 0.19302)};
          // modelBaseNames =
      break;
    }
    case('U'): {
      modelBaseCoords = {Vector(1.95363, -1.45659, -0.19073), Vector(0.75787, -0.57587, -0.07419),   
                        Vector(0.97215, 0.78019, -0.13405), Vector(-0.49400, -1.10811, 0.08633),
                        Vector(-0.15841, 1.56515, -0.02137), Vector(-1.45350, 1.11665, 0.14102),
                        Vector(-1.57773, -0.32139, 0.19302)};

      break;
    }
  }
  // Construct the Matrices U and Ut
  //U = {Vector(0, 0, 0), Vector(0, 0, 0), Vector(0, 0, 0)};
  //Ut = {Vector(0, 0, 0), Vector(0, 0, 0), Vector(0, 0, 0)};
  double weight = 1.00/((double)natoms);

  for(i=0; i<3; i++){
    for(j=0; j<3; j++){
      for(n=0; n<natoms; n++){
        T[i][j] += weight * atoms[n][j] * modelBaseCoords[n][i];
      }
      Tt[j][i] = T[i][j];
    }
  }

  DetT = Determinant3x3( T );

  // Construct the Matrices UxUt and UtxU
  for(i=0; i<3; i++){
    for(j=0; j<3; j++){
      TxTt[i][j] = T[0][i]*Tt[j][0] + T[1][i]*Tt[j][1] + T[2][i]*Tt[j][2];
      TtxT[i][j] = Tt[0][i]*T[j][0] + Tt[1][i]*T[j][1] + Tt[2][i]*T[j][2];
    }
  }

  // calculate the eigenVectors and Values
  EigenValsAndVectorsSymmetric(TxTt, &eigenValues, &eigenVectorsH);
  eigenValuesDeSquared = { sqrt(eigenValues[0]) ,sqrt(eigenValues[1]) ,sqrt(eigenValues[2]) };

  //EigenValsAndVectorsSymmetric(TtxT, &eigenValues, &eigenVectorsK);
  for ( i = 0; i < 2; i++) {
    for ( j = 0; j < 3; j++) {
      eigenVectorsK[i][j] = Tt[0][j]*eigenVectorsH[i][0] + Tt[2][j]*eigenVectorsH[i][2] + Tt[1][j]*eigenVectorsH[i][1];
      }
    eigenVectorsK[i] = eigenVectorsK[i]/eigenVectorsK[i].modulo();
  }
  // construct rotation Matrix

  if(DetT >= 0){
    eigenVectorsK[2] = crossProduct(eigenVectorsK[0], eigenVectorsK[1]);
    eigenVectorsK[2] = eigenVectorsK[2]/eigenVectorsK[2].modulo();
    for(i=0; i<3; i++){
      for(j=0; j<3; j++){
        R[i][j] = eigenVectorsK[0][i]*eigenVectorsH[0][j] + eigenVectorsK[1][i]*eigenVectorsH[1][j] + eigenVectorsK[2][i]*eigenVectorsH[2][j];
      }
    }
  }else{
    eigenVectorsK[2] = -1*crossProduct(eigenVectorsK[0], eigenVectorsK[1]);
    eigenVectorsK[2] = eigenVectorsK[2]/eigenVectorsK[2].modulo();
    for(i=0; i<3; i++){
      for(j=0; j<3; j++){
        R[i][j] = eigenVectorsK[0][i]*eigenVectorsH[0][j] + eigenVectorsK[1][i]*eigenVectorsH[1][j] - eigenVectorsK[2][i]*eigenVectorsH[2][j];
      }
    }
  }
  // Calculate best fit of the Three needed Atoms 
  for(i=0; i<3; i++){
    fitAtom1[i] = R[0][i]*modelBaseCoords[0][0] + R[1][i]*modelBaseCoords[0][1] + R[2][i]*modelBaseCoords[0][2]; // rotate the model N9 or N1 (=0) base to be as close as possible to the simulation atom
    fitAtom2[i] = R[0][i]*modelBaseCoords[1][0] + R[1][i]*modelBaseCoords[1][1] + R[2][i]*modelBaseCoords[1][2]; // rotate the model C*1 (=1) base to be as close as possible to the simulation atom
    fitAtom3[i] = R[0][i]*modelBaseCoords[2][0] + R[1][i]*modelBaseCoords[2][1] + R[2][i]*modelBaseCoords[2][2]; // rotate the model C4 or C2 (=2) base to be as close as possible to the simulation atom
  }

  *fit1=(fitAtom1+MC);
  *fit2=(fitAtom2+MC);
  *fit3=(fitAtom3+MC);

  return;
}


double cirtoNA::Determinant3x3(vector< VectorGeneric<3> > Matrix){
  
  double D = Matrix[0][0]*Matrix[1][1]*Matrix[2][2] + Matrix[1][0]*Matrix[2][1]*Matrix[0][2] + Matrix[2][0]*Matrix[0][1]*Matrix[1][2];
  D -= Matrix[0][0]*Matrix[2][1]*Matrix[1][2] + Matrix[1][1]*Matrix[0][2]*Matrix[2][0] + Matrix[0][1]*Matrix[1][0]*Matrix[2][2];
  return D;
}


void cirtoNA::EigenValsAndVectorsSymmetric(vector<VectorGeneric<3>> Matrix, VectorGeneric<3>* EigenValues, vector<VectorGeneric<3>>* EigenVectors){
  VectorGeneric<3> eigVals;
  double temp, p1, q, p2, p, r, phi;
  int i, j, k;
  vector<VectorGeneric<3>> eye = {Vector(1, 0, 0), Vector(0, 1, 0), Vector(0, 0, 1)};
  p1 = Matrix[1][0]*Matrix[0][1] + Matrix[2][0]*Matrix[0][2] + Matrix[2][1]*Matrix[1][2];
  q = (Matrix[0][0]+Matrix[1][1]+Matrix[2][2])/3;               // trace(A) is the sum of all diagonal values
  p2 = (Matrix[0][0] - q)*(Matrix[0][0] - q) + (Matrix[1][1] - q)*(Matrix[1][1] - q) + (Matrix[2][2] - q)*(Matrix[2][2] - q) + 2*p1;
  p = sqrt(p2 / 6);
  vector<VectorGeneric<3>> B;
  B.resize(3);
  //B = (1 / p) * (A - q * I)    I is the identity matrix
  for(j=0; j<3; j++){
    for(i=0; i<3; i++){
      B[i][j] = (Matrix[i][j]-q*eye[i][j])/p;
    }
  }
    //B = (1 / p) * (A - q * I)    I is the identity matrix
  r = Determinant3x3( B ) / 2;

  if (r <= -1){
    phi = pi / 3;
  }else if(r >= 1){
    phi = 0;
  }else {
    phi = acos(r) / 3;
  }
  // the eigenvalues satisfy eig3 <= eig2 <= eig1

  eigVals[0] = p * 2 * cos(phi ) + q ;
  eigVals[2] = p * 2 * cos(phi + 2*pi/3 ) + q ;
  eigVals[1] = 3 * q - eigVals[0] - eigVals[2];

    //now eigenVectors
  VectorGeneric<3> EigenVector1, EigenVector2, EigenVector3, nonZeroColl, nonParrallelVec;
  VectorGeneric<3> checkVector;
  
  vector<VectorGeneric<3>> M = Matrix;

  for(i=0; i<3;  i++){
    M[i][i] -= eigVals[0]; // Matrix -I*lambda
  }

  for(i=0; i<3; i++){
    EigenVector1 = crossProduct(M[i%3], M[(i+1)%3]);
    if(EigenVector1.modulo()>=0.000001){ //two collumns were independent, and thus calculated eigen vector is in null space of M and eigen space of matrix
      EigenVector1 = EigenVector1/EigenVector1.modulo();
      break;
    }
  }

  if(EigenVector1.modulo()<0.000001) { //no two independent collumns
    if(!(eigVals[0]==eigVals[1]) ){
      eigVals[2] = eigVals[1];
      eigVals[1] = eigVals[0];
    }
    for (int i = 0; i < 3; ++i) {
      nonZeroColl = M[i];
      if(nonZeroColl.modulo()>=0.000001){
        if( (nonZeroColl[1]==0) && (nonZeroColl[2]==0) ) {
          nonParrallelVec = {0, 1, 0};
        }else{
          nonParrallelVec = {1, 0, 0};
        }  
        EigenVector1 = crossProduct(nonZeroColl, nonParrallelVec);
        EigenVector2 = crossProduct(EigenVector1, nonZeroColl);
        break;
      }
    }
    EigenVector3 = crossProduct(EigenVector1, EigenVector2);

    EigenVector1 = EigenVector1/EigenVector1.modulo();
    EigenVector2 = EigenVector2/EigenVector2.modulo();
    EigenVector3 = EigenVector3/EigenVector3.modulo();

    *EigenValues = eigVals;
    *EigenVectors = {EigenVector1, EigenVector2, EigenVector3};
    return;
  }


  M = Matrix;
  for(i=0; i<3;  i++){
    M[i][i] -= eigVals[1];
  }

  for(i=0; i<3; i++){
    EigenVector2 = crossProduct(M[i%3], M[(i+1)%3]);
    if(EigenVector2.modulo()>=0.0001){ //two collumns were independent, and thus calculated eigen vector is in null space of M and eigen space of matrix
      EigenVector2 = EigenVector2/EigenVector2.modulo();
      break;
    }
  }

  if(EigenVector2.modulo()<0.0001) { //no two independent collumns
    eigVals[2] = eigVals[1];
    for (int i = 0; i < 3; ++i) {
      nonZeroColl = M[i];
      if(nonZeroColl.modulo()>=0.0001){
        if( (nonZeroColl[1]==0) && (nonZeroColl[2]==0) ){
          nonParrallelVec = {0, 1, 0};
        }else{
          nonParrallelVec = {1, 0, 0};
        }
        EigenVector2 = crossProduct(nonZeroColl, nonParrallelVec);
        break;
      }
    }
    EigenVector3 = crossProduct(EigenVector1, EigenVector2);

    EigenVector1 = EigenVector1/EigenVector1.modulo();
    EigenVector2 = EigenVector2/EigenVector2.modulo();
    EigenVector3 = EigenVector3/EigenVector3.modulo();

    *EigenValues = eigVals;
    *EigenVectors = {EigenVector1, EigenVector2, EigenVector3};

    return;
  }

  EigenVector3 = crossProduct(EigenVector1, EigenVector2);

  EigenVector1 = EigenVector1/EigenVector1.modulo();
  EigenVector2 = EigenVector2/EigenVector2.modulo();
  EigenVector3 = EigenVector3/EigenVector3.modulo();

  *EigenValues = eigVals;
  *EigenVectors = {EigenVector1, EigenVector2, EigenVector3};

  return;
}


// calculator
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//helperfunctions
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double cirtoNA::Norm(VectorGeneric<3> vec){
  double norm = vec.modulo();
  return norm;
}

VectorGeneric<3> cirtoNA::Normalize(VectorGeneric<3> vec){
  double norm = vec.modulo();
  VectorGeneric<3> normalized = vec/norm;
  return normalized;
}

VectorGeneric<3> cirtoNA::AxisRotation(VectorGeneric<3> v, VectorGeneric<3> k, double theta){
  VectorGeneric<3> vrot = v*cos(theta) + crossProduct(k,v)*sin(theta) + k*dotProduct(k,v)*(1-cos(theta));
  return vrot;
}

vector<VectorGeneric<3>> cirtoNA::Skew(VectorGeneric<3> phi){
  vector<VectorGeneric<3>> S(3, Vector(0,0,0));
  S = {Vector(0, phi[2], -1*phi[1]), Vector(-1*phi[2], 0, phi[0]), Vector(phi[1], -1*phi[0], 0)};

  return S;
}

vector<VectorGeneric<3>> cirtoNA::RotationMatrix(vector<VectorGeneric<3>> frame1, vector<VectorGeneric<3>> frame2){

  vector<VectorGeneric<3>> Q(3, Vector(0,0,0));
  for(int i = 0; i<3; i++){
    for(int j = 0; j<3; j++){
      Q[i][j] = 0;
      for(int k = 0; k<3; k++){
        Q[i][j] += frame2[k][j]*frame1[k][i];
        //dotProduct(frame1[i],frame2[j]); //altered 13.14.05 matrix product indices verwisseld
      }
    }
  }

  return Q;
}

VectorGeneric<3> cirtoNA::RotationAxis(vector<VectorGeneric<3>> Q){

  VectorGeneric<3> axis;
  axis[0] = Q[1][2]-Q[2][1];
  axis[1] = Q[2][0]-Q[0][2];
  axis[2] = Q[0][1]-Q[1][0];

  return axis;
}

double cirtoNA::RotationAngle(vector<VectorGeneric<3>> Q){

  double cost=(Q[0][0]+Q[1][1]+Q[2][2]-1)/2;
  if(-1>cost){cost = -0.99999999999999;}
  if(1<cost){cost = 0.99999999999999;}
  double theta=acos(cost);

  return theta;
}

VectorGeneric<3> cirtoNA::RotationVector(vector<VectorGeneric<3>> Q){
  VectorGeneric<3> Axis, axis;
  double theta;

  Axis = RotationAxis(Q);
  axis = Normalize(Axis);
  theta = RotationAngle(Q);

  return theta*axis;
}

vector<VectorGeneric<3>> cirtoNA::AxisAngleToMatrix(VectorGeneric<3> RotationVector){
  vector<VectorGeneric<3>> Q(3, Vector(0,0,0));
  double theta = RotationVector.modulo();
  VectorGeneric<3> axis = RotationVector/theta;
  double cost=cos(theta);
  double sint=sin(theta);

  for(int i=0; i<3; i++){
    Q[i][i] = cost + axis[i]*axis[i]*(1-cost);
    Q[i][(i+1)%3] = axis[i]*axis[(i+1)%3]*(1-cost) + axis[(i+2)%3]*sint;
    Q[i][(i+2)%3] = axis[i]*axis[(i+2)%3]*(1-cost) - axis[(i+1)%3]*sint;
  }

 return Q;
}

vector<VectorGeneric<3>> cirtoNA::midFrame(vector<VectorGeneric<3>> frame1, VectorGeneric<3> RotationVector){
  vector<VectorGeneric<3>> mid(3, Vector(0,0,0));
  double theta = RotationVector.modulo();
  VectorGeneric<3> axis = RotationVector/theta;

  for(int i=0; i<3; i++){
    mid[i]=AxisRotation(frame1[i], axis, theta/2); //altered last 13.15 replaced -theta by theta
    mid[i]=mid[i]/mid[i].modulo();
  }

  return mid;
}

vector<VectorGeneric<3>> cirtoNA::BaseFrame(VectorGeneric<3> glycosidic, VectorGeneric<3> internal){
  vector<VectorGeneric<3>> frame(3, Vector(0,0,0));
  VectorGeneric<3> normalVector, longVector, shortVector, Gly, Int;
  double normnormal, normlong, normshort, normgly, normint;

  normgly=glycosidic.modulo();
  normint=internal.modulo();

  Gly=glycosidic/normgly;
  Int=internal/normint;

  normalVector=crossProduct(Gly, Int);
  normnormal=normalVector.modulo();
  frame[2] = normalVector/normnormal;
  
  longVector = AxisRotation(Gly, frame[2], -0.9496336460101146); // altered in 14.20 version (angle accuracy increased)
  normlong=longVector.modulo();
  frame[1] = longVector/normlong;

  shortVector = crossProduct(frame[1], frame[2]); // altered in 14.20 version (vectors swapped)
  normshort=shortVector.modulo();
  frame[0] = shortVector/normshort;

  return frame;

}

VectorGeneric<3> cirtoNA::RefPoint(VectorGeneric<3> glycosidic, VectorGeneric<3> normal, VectorGeneric<3> nitrogen){
  VectorGeneric<3> displacement, position;

  displacement = AxisRotation( glycosidic/glycosidic.modulo(), normal, 2.4691172928);
  position = nitrogen + 0.4702*(displacement/displacement.modulo());
  return position;

}

//Let's redo the derivatives
// Derivative function convention: deriv function does not take derivs as in put
// name of function is D[Func_Name]D[Input] for derivative of [Func_Name] to [Input]
// Output is given as an array of shape([Input])xshape([Func_Name])
// Input is given as the Output of [Func_Name] followed by its Inputs (in the same order) 
// all output is done directly (no refferencing)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// vector basics
vector<VectorGeneric<3>> cirtoNA::DNormalizeDvec(VectorGeneric<3> normedVec, VectorGeneric<3> vec){
  double norm = vec.modulo() ;
  VectorGeneric<3> derNormedx = (norm*Vector(1,0,0)-normedVec*vec[0])/(norm*norm);
  VectorGeneric<3> derNormedy = (norm*Vector(0,1,0)-normedVec*vec[1])/(norm*norm);
  VectorGeneric<3> derNormedz = (norm*Vector(0,0,1)-normedVec*vec[2])/(norm*norm);
  vector<VectorGeneric<3>> derNormed = {derNormedx, derNormedy, derNormedz} ;

  return derNormed;
}

VectorGeneric<3> cirtoNA::DNormDvec(double norm, VectorGeneric<3> vec){
  VectorGeneric<3> Der = Normalize(vec); 
  return Der;  // D(sqrt(x**2+y**2+z**2))Dx = 2x/2sqrt(x**2+y**2+z**2) = x/sqrt(x**2+y**2+z**2)
}

// construction of axis angle

vector<VectorGeneric<3>> cirtoNA::DRotationAngleDQ(double theta, vector<VectorGeneric<3>> Q){
  //calculation done from rotation matrix, angle axis correpondence
  vector<VectorGeneric<3>> Ders(3, Vector(0,0,0));

  for(int i =0; i<3; i++){
    Ders[i][i] = -0.5/sin(theta);
  }

  return Ders;
}

vector<vector<VectorGeneric<3>>> cirtoNA::DRotationAxisDQ(VectorGeneric<3> Axis, vector<VectorGeneric<3>> Q){
  vector<VectorGeneric<3>> dummy33(3, Vector(0,0,0));
  vector<vector<VectorGeneric<3>>> Ders(3, dummy33);

  for(int i =0; i<3; i++){
    Ders[(i+1)%3][(i+2)%3][i] = 1;
    Ders[(i+2)%3][(i+1)%3][i] = -1;
  }

  return Ders;
}

vector<vector<VectorGeneric<3>>> cirtoNA::DVectorDTransformation(VectorGeneric<3> RotationVector, vector<VectorGeneric<3>> Q){
  int k, l, i, j;
  double theta, sint, cost, cott;
  VectorGeneric<3> Axis, axis;
  vector<VectorGeneric<3>> Zero3x3(3, Vector(0,0,0)), daxisdAxis(3, Vector(0,0,0)), dthetadQ(3, Vector(0,0,0));

  vector<vector<VectorGeneric<3>>> dAxisdQ(3, Zero3x3), dVectordQ(3, Zero3x3);
  
  Axis = RotationAxis(Q);
  theta = RotationAngle(Q);
  sint = sin(theta);
  cost = cos(theta);
  cott = cost/sint;

  // Calculate derivatives
  dAxisdQ = DRotationAxisDQ(Axis, Q);
  dthetadQ = DRotationAngleDQ(theta, Q);

  //resolve chian rule between for normalizing the axis

  for( k = 0; k < 3; k++){
    for( l = 0; l < 3; l++){
      for( i = 0; i < 3; i++){
        dVectordQ[k][l][i] += theta*dAxisdQ[k][l][i];
      }
    }
  }

  //resolve product rule for multiplying by theta

  for( k = 0; k < 3; k++){
    l=k;  // theta only has non zero derivatives with respect to the diagonal
    for( i = 0; i < 3; i++){
      dVectordQ[k][l][i] += (1-theta*cott)*dthetadQ[k][l]*Axis[i];
    }
  }


  for( k = 0; k < 3; k++){
    for( l = 0; l < 3; l++){
      for( i = 0; i < 3; i++){
          dVectordQ[k][l][i] = dVectordQ[k][l][i]/(2*sint);
      }
    }
  }

  return dVectordQ;
}

//Rotation matrices that are functions of frames derived with respect to rotations of those frames

vector<vector<VectorGeneric<3>>> cirtoNA::DTransformationDRot(vector<VectorGeneric<3>> frame1, vector<VectorGeneric<3>> frame2, double direction){
  vector<VectorGeneric<3>> Sx(3, Vector(0,0,0)), Sy(3, Vector(0,0,0)), Sz(3, Vector(0,0,0)), dummy33(3, Vector(0,0,0)), DtransformDrotx(3, Vector(0,0,0)), DtransformDroty(3, Vector(0,0,0)), DtransformDrotz(3, Vector(0,0,0));
  vector<VectorGeneric<3>> Sxframe1t(3, Vector(0,0,0)), Syframe1t(3, Vector(0,0,0)), Szframe1t(3, Vector(0,0,0));
  vector<vector<VectorGeneric<3>>> DtransformDrotation(3, dummy33);

  Sx = Skew(Vector(1, 0, 0));
  Sy = Skew(Vector(0, 1, 0));
  Sz = Skew(Vector(0, 0, 1));
  
  if (direction==-1){
    for(int k = 0; k<3; k++){
      for(int l = 0; l<3; l++){
        for(int j =0; j<3; j++){
          Sxframe1t[l][k] += Sx[j][k]*frame1[j][l];
          Syframe1t[l][k] += Sy[j][k]*frame1[j][l];
          Szframe1t[l][k] += Sz[j][k]*frame1[j][l];
        }
      }
    }
  }else if(direction==1){
    for(int k = 0; k<3; k++){
      for(int l = 0; l<3; l++){
        for(int j = 0; j<3; j++){
          Sxframe1t[l][k] += Sx[k][j]*frame1[j][l];
          Syframe1t[l][k] += Sy[k][j]*frame1[j][l];
          Szframe1t[l][k] += Sz[k][j]*frame1[j][l];
        }
      }
    }
  }


  for(int k = 0; k<3; k++){
    for(int l = 0; l<3; l++){
      for(int j = 0; j<3; j++){
        DtransformDrotx[k][l] +=  frame2[j][l]*Sxframe1t[k][j];
        DtransformDroty[k][l] +=  frame2[j][l]*Syframe1t[k][j];
        DtransformDrotz[k][l] +=  frame2[j][l]*Szframe1t[k][j];
      }
    }
  }

  DtransformDrotation = {DtransformDrotx, DtransformDroty, DtransformDrotz};

  return DtransformDrotation;
}

vector<vector<VectorGeneric<3>>> cirtoNA::DMidDRot(vector<VectorGeneric<3>>  midFrame, vector<VectorGeneric<3>> frame){
  vector<VectorGeneric<3>> Sx(3, Vector(0,0,0)), Sy(3, Vector(0,0,0)), Sz(3, Vector(0,0,0)), DmidDrotx(3, Vector(0,0,0)), DmidDroty(3, Vector(0,0,0)), DmidDrotz(3, Vector(0,0,0));
  vector<VectorGeneric<3>> dummy33(3, Vector(0,0,0)), SxM(3, Vector(0,0,0)), SyM(3, Vector(0,0,0)), SzM(3, Vector(0,0,0)), MtSxM(3, Vector(0,0,0)), MtSyM(3, Vector(0,0,0)), MtSzM(3, Vector(0,0,0));
  vector<vector<VectorGeneric<3>>> DmidDrotation(3, dummy33);
  vector<VectorGeneric<3>> I = {Vector(1,0,0), Vector(0,1,0), Vector(0,0,1)};
  VectorGeneric<3> phi;

  Sx = Skew(Vector(1, 0, 0));
  Sy = Skew(Vector(0, 1, 0));
  Sz = Skew(Vector(0, 0, 1));
  
  for(int k = 0; k<3; k++){
    for(int l = 0; l<3; l++){
      for(int j = 0; j<3; j++){
        DmidDrotx[k][l] += 0.5*midFrame[j][l]*Sx[k][j];
        DmidDroty[k][l] += 0.5*midFrame[j][l]*Sy[k][j];
        DmidDrotz[k][l] += 0.5*midFrame[j][l]*Sz[k][j];
      }
    }
  }
  DmidDrotation = {DmidDrotx, DmidDroty, DmidDrotz};

  return DmidDrotation;
}

//actual calculation
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void cirtoNA::calculate(){
  int i, j, l, k, n, N, M, x;

  switch(cirtoNA::INPUTMODE){

 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 ///////////////////////////////////////////Pair//////////////////////////////////////////////////////////////////////////////////////// 
    case(pair):{  // To be converted back from new sequence setup
      //initialize
      break;
      }
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 ///////////////////////////////////////////Step////////////////////////////////////////////////////////////////////////////////////////      
    case(step):{  // to be redone (simplified) from Sequence mode
      //initialize
      break;
    }
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 /////////////////////////////////////////Sequence//////////////////////////////////////////////////////////////////////////////////////

    case(sequence):{
     //initialize
      //define needed variables

      int npairs, nbases, nsteps, natoms;  // global numbers
      int pair, base, step;  // counters

     //

     // get size of chain for initialization
      nbases = BASES.size();
      npairs = nbases/2;
      nsteps = npairs-1;

     //

     //set size of vector variables and initialize correctly

      // dummys to set size of other varibles
      VectorGeneric<3> Zero3 = Vector(0, 0, 0);
      vector<double> Zero10(10, 0.0);
      vector<VectorGeneric<3>> Zero2x3(2, Vector(0, 0, 0)), Zero3x3(3, Vector(0, 0, 0)), Zero10x3(10, Vector(0, 0, 0));
      vector<vector<VectorGeneric<3>>> Zero2x2x3(2, Zero2x3), Zero3x3x3(3, Zero3x3), Zero2x3x3(2, Zero3x3);
      vector<vector<vector<VectorGeneric<3>>>> Zero3x3x3x3(3, Zero3x3x3), Zero2x3x3x3(2, Zero3x3x3), Zero2x2x3x3(2, Zero2x3x3);
      vector<vector<vector<vector<VectorGeneric<3>>>>> Zero2x2x3x3x3(2, Zero2x3x3x3);

      vector<double> tilt(nsteps, 0), roll(nsteps, 0), twist(nsteps, 0), shift(nsteps, 0), slide(nsteps, 0), rise(nsteps, 0);  // output Steps
      vector<double> buckle(npairs, 0), propeller(npairs, 0), opening(npairs, 0), shear(npairs, 0), stretch(npairs, 0), stagger(npairs, 0); // output Pairs

      VectorGeneric<3> glycosidic, internal;  // important intermediate variables to get to the base frames

      // information related to bases
      vector<int> AtomsInBase(nbases, 0);
      vector<vector<double>> AtomsMasses(nbases, Zero10);
      vector<VectorGeneric<3>> BaseCenter(nbases, Vector(0, 0, 0)), BaseRefs(nbases, Vector(0, 0, 0));
      vector<vector<VectorGeneric<3>>> AtomsBases(nbases, Zero10x3), refAtoms(nbases, Zero3x3), BaseRedefs(nbases, Zero3x3), Baseframes(nbases, Zero3x3);

      // information on Pairs and Steps
      vector<vector<VectorGeneric<3>>> Pairframes(npairs, Zero3x3), Stepframes(nsteps, Zero3x3), TransformationMatrixPair(npairs, Zero3x3), TransformationMatrixStep(nsteps, Zero3x3);
      vector<VectorGeneric<3>> RotationvectorPair(npairs, Vector(0, 0, 0)), pairTranslations(npairs, Vector(0, 0, 0)), RotationvectorStep(nsteps, Vector(0, 0, 0)), stepTranslations(nsteps, Vector(0, 0, 0));

      // derivatives of Pair and step mids and rotations to underlying baseframe rotations
      vector<vector<vector<VectorGeneric<3>>>> DRotationVectorDPairTransformation(npairs, Zero3x3x3), DRotationVectorDStepTransformation(nsteps, Zero3x3x3);
     
      vector<vector<vector<VectorGeneric<3>>>> DRotationVectorDPairRotation(npairs, Zero2x3x3);
      vector<vector<vector<vector<VectorGeneric<3>>>>> DRotationVectorDStepRotation(nsteps, Zero2x2x3x3);

      vector<vector<vector<vector<VectorGeneric<3>>>>> DPairTransformationDRotation(npairs,Zero2x3x3x3), DPairMidDRotation(npairs,Zero2x3x3x3);
      vector<vector<vector<vector<vector<VectorGeneric<3>>>>>> DStepTransformationDRotation(nsteps,Zero2x2x3x3x3), DStepMidDRotation(nsteps,Zero2x2x3x3x3);

      // derivatives giving torques on the bqse frame
      vector<VectorGeneric<3>> DtiltDBaseRot(4, Vector(0, 0, 0)), DrollDBaseRot(4, Vector(0, 0, 0)), DtwistDBaseRot(4, Vector(0, 0, 0));
      vector<VectorGeneric<3>> DbuckleDBaseRot(2, Vector(0, 0, 0)), DpropellerDBaseRot(2, Vector(0, 0, 0)), DopeningDBaseRot(2, Vector(0, 0, 0));

      // derivatives used to transfer torques on the base frame in to forces
      vector<vector<VectorGeneric<3>>> DtiltDbaseframe(4, Zero3x3), DrollDbaseframe(4, Zero3x3), DtwistDbaseframe(4, Zero3x3);
      vector<vector<VectorGeneric<3>>> DbuckleDbaseframe(2, Zero3x3), DpropellerDbaseframe(2, Zero3x3), DopeningDbaseframe(2, Zero3x3);

      // output derivatives
      VectorGeneric<3>  DtiltDatom, DrollDatom, DtwistDatom, DbuckleDatom, DpropellerDatom, DopeningDatom;

      vector<VectorGeneric<3>> HalfRotationMatrix(3, Vector(0, 0, 0));

     //
      
     // initialize pointers for output
      vector<Value*> valuestilt, valuesroll, valuestwist;
      vector<Value*> valuesshift, valuesslide, valuesrise;
      valuestilt.resize(nsteps); valuesroll.resize(nsteps); valuestwist.resize(nsteps); 
      valuesshift.resize(nsteps); valuesslide.resize(nsteps); valuesrise.resize(nsteps);

      vector<Value*> valuesbuckle, valuespropeller, valuesopening;
      vector<Value*> valuesshear, valuesstretch, valuesstagger;
      valuesbuckle.resize(npairs); valuespropeller.resize(npairs); valuesopening.resize(npairs); 
      valuesshear.resize(npairs); valuesstretch.resize(npairs); valuesstagger.resize(npairs);

     //

     // pbc handling

      if(pbc){makeWhole();}

     //

     // fitting and retrieving atoms
      N=0;
      for(base=0;base<(nbases-1);base=base+2){ // for all pairs in the DNA chain check the first base, second base should be complementary
        switch(BASES[base]){  //statements to load the required atoms depinding on the sequence provided

          case('A'): {
            AtomsInBase[base]=10;
            AtomsBases[base].resize(10);
            AtomsMasses[base].resize(10);

            for (n = 0; n < 10; n++) {
              AtomsBases[base][n] = getPosition(N);
              AtomsMasses[base][n] = getMass(N);
              N++;
            }

            AtomsInBase[base+1]=7;
            AtomsBases[base+1].resize(7);
            AtomsMasses[base+1].resize(7);

            for (n = 0; n < 7; n++) {
              AtomsBases[base + 1][n] = getPosition(N);
              AtomsMasses[base + 1][n] = getMass(N);
              N++;
            }

            break;
          }

          case('C'): {
            AtomsInBase[base]=7;
            AtomsBases[base].resize(7);
            AtomsMasses[base].resize(7);

            for (n = 0; n < 7; n++) {
              AtomsBases[base][n] = getPosition(N);
              AtomsMasses[base][n] = getMass(N);
              N++;
            }

            AtomsInBase[base+1]=10;
            AtomsBases[base+1].resize(10);
            AtomsMasses[base+1].resize(10);

            for (n = 0; n < 10; n++) {
              AtomsBases[base + 1][n] = getPosition(N);
              AtomsMasses[base + 1][n] = getMass(N);
              N++;
            }

            break;
          }

          case('G'): {
            AtomsInBase[base]=10;
            AtomsBases[base].resize(10);
            AtomsMasses[base].resize(10);

            for (n = 0; n < 10; n++) {
              AtomsBases[base][n] = getPosition(N);
              AtomsMasses[base][n] = getMass(N);
              N++;
            }

            AtomsInBase[base+1]=7;
            AtomsBases[base+1].resize(7);
            AtomsMasses[base+1].resize(7);

            for (n = 0; n < 7; n++) {
              AtomsBases[base + 1][n] = getPosition(N);
              AtomsMasses[base + 1][n] = getMass(N);
              N++;
            }

            break;
          }

          case('T'): {
            AtomsInBase[base]=7;
            AtomsBases[base].resize(7);
            AtomsMasses[base].resize(7);

            for (n = 0; n < 7; n++) {
              AtomsBases[base][n] = getPosition(N);
              AtomsMasses[base][n] = getMass(N);
              N++;
            }

            AtomsInBase[base+1]=10;
            AtomsBases[base+1].resize(10);
            AtomsMasses[base+1].resize(10);

            for (n = 0; n < 10; n++) {
              AtomsBases[base + 1][n] = getPosition(N);
              AtomsMasses[base + 1][n] = getMass(N);
              N++;
            }

            break;
          }

          case('U'): {
            AtomsInBase[base]=7;
            AtomsBases[base].resize(7);
            AtomsMasses[base].resize(7);

            for (n = 0; n < 7; n++) {
              AtomsBases[base][n] = getPosition(N);
              AtomsMasses[base][n] = getMass(N);
              N++;
            }

            AtomsInBase[base+1]=10;
            AtomsBases[base+1].resize(10);
            AtomsMasses[base+1].resize(10);

            for (n = 0; n < 10; n++) {
              AtomsBases[base + 1][n] = getPosition(N);
              AtomsMasses[base + 1][n] = getMass(N);
              N++;
            }

            break;
          }
        }
        if(FITTED){  // if fitting was requested fit every base and return the three relavant positions for the construction of the base frames
          fit(BASES[base], AtomsBases[base], &refAtoms[base][0], &refAtoms[base][1], &refAtoms[base][2]);
          fit(BASES[base+1], AtomsBases[base+1], &refAtoms[base+1][0], &refAtoms[base+1][1], &refAtoms[base+1][2]);
          
        }else{  // if no fitting was requested the first three atoms give the relevant positions for constructing the frames
          refAtoms[base] = {AtomsBases[base][0],AtomsBases[base][1],AtomsBases[base][2]};
          refAtoms[base+1] = {AtomsBases[base+1][0],AtomsBases[base+1][1],AtomsBases[base+1][2]};
        }
      }
      natoms = N; //identify the total number of atoms used

     //

     // create base frames

      for(base=0;base<nbases;base++){  // for every base create a referrence frame and point
        glycosidic=delta(refAtoms[base][1],refAtoms[base][0]);  // define the glycosidic bond
        internal=delta(refAtoms[base][1],refAtoms[base][2]);  // get the internal bond used in frame construction
        
        Baseframes[base] = BaseFrame(glycosidic, internal);  // construct the base frame
        BaseRefs[base] = RefPoint(glycosidic, Baseframes[base][2], AtomsBases[base][1]); // the refference point is defined as a rotation of the glycsidic about the normal
        BaseCenter[base] = Vector(0,0,0);
        for(n=0;n<AtomsInBase[base];n++){
          BaseCenter[base] += AtomsBases[base][n];
        }
        BaseCenter[base] = BaseCenter[base]/AtomsInBase[base];  // get the geometric center of the base

        if((base%2) == 1){  // deal with the 3'->5', 5'->3' anti-symmetry making the tangents of opposing bases point in the same direction()
            Baseframes[base][2] = -1*Baseframes[base][2]; //making the tangents of opposing bases point in the same direction
            Baseframes[base][1] = -1*Baseframes[base][1]; //also making the e2 vectors (connecting backbones) not point against but "with" each other
        }
        //define atoms as function of frame
        
        BaseRedefs[base].resize(AtomsInBase[base]);
        for(n=0;n<AtomsInBase[base];n++){
          for(i=0;i<3;i++){
            BaseRedefs[base][n][i] = dotProduct(AtomsBases[base][n]-BaseRefs[base], Baseframes[base][i]);  // for future reference we calculate the position of the atoms in the base frame
          }
        }
      }    
         
     //

     // a diagram for refference:

      //  Left    High    Right
      //   ||              ||
      //   ||==2  pair  3==||
      //   ||              || 
      //   ||     Step     ||
      //   ||              ||
      //   ||==0  pair  1==|| 
      //   ||              ||
      //  Left    Low     Right
      //  frame1(pair) = Left, frame1(Step) = Low

     //


     // apply rodrigues to frames

      for(pair=0;pair<npairs;pair++){ //input of loop is two base frames, output are the rotation and translation vector between them and their mid frame, appropriate derivatives are also set
        TransformationMatrixPair[pair] = RotationMatrix(Baseframes[2*pair], Baseframes[2*pair+1]);  // the rotation matrix that transforms the left base to the right (given in the lab frame)
        RotationvectorPair[pair] = RotationVector(TransformationMatrixPair[pair]); // convert the rotation to axis angle representation
        Pairframes[pair] = midFrame(Baseframes[2*pair], RotationvectorPair[pair]); // find the frame halfway between the two base frames
        pairTranslations[pair] = BaseRefs[2*pair+1]-BaseRefs[2*pair];  // Calculate the translation vector between the bases

        //log.printf("RotationvectorPair of pair %d = %f, %f, %f\n", pair, RotationvectorPair[pair][0],RotationvectorPair[pair][1],RotationvectorPair[pair][2]);  
        //log.printf("Pairframes of pair %d =  %f, %f, %f\n", pair, Pairframes[pair][0][0], Pairframes[pair][0][1],Pairframes[pair][1][2]);  

        //bookkeepderivatives

        DRotationVectorDPairTransformation[pair] = DVectorDTransformation(RotationvectorPair[pair], TransformationMatrixPair[pair]); // calculate the derivative from transforming the matrix into axis angle

        for(base=0;base<2;base++){
          DPairTransformationDRotation[pair][base] = DTransformationDRot(Baseframes[2*pair], Baseframes[2*pair+1], pow(-1.0,base+1)); // calculate the derivative of the transformation to rotations of the base frames

          for(i=0;i<3;i++){
            for(x=0;x<3;x++){
              for(k=0;k<3;k++){
                for(l=0;l<3;l++){
                  DRotationVectorDPairRotation[pair][base][x][i] += DPairTransformationDRotation[pair][base][x][k][l]*DRotationVectorDPairTransformation[pair][k][l][i];  //this is the derivative of a vector (axis-angle) to rotations of the base frames through the chain rule
                }
              }
            }
          }

          DPairMidDRotation[pair][base] = DMidDRot(Pairframes[pair], Baseframes[2*pair+base]); // get the derivative of the midframe to baseframe rotations
        }
      }
   
  
      for(step=0;step<nsteps;step++){ //input of loop is two pair frames, output is the rotation and translation vector between them and their mid frame, appropriate derivatives are also set
        TransformationMatrixStep[step] = RotationMatrix(Pairframes[step], Pairframes[step+1]);
        RotationvectorStep[step] = RotationVector(TransformationMatrixStep[step]);
        Stepframes[step] = midFrame(Pairframes[step], RotationvectorStep[step]);
        stepTranslations[step] = 0.5*(BaseRefs[2*step+2]+BaseRefs[2*step+3] - BaseRefs[2*step]-BaseRefs[2*step+1]);

        //bookkeepderivatives
        HalfRotationMatrix = AxisAngleToMatrix(0.5*RotationvectorStep[step]);

        DRotationVectorDStepTransformation[step] = DVectorDTransformation(RotationvectorStep[step], TransformationMatrixStep[step]);
       
        for(base=0;base<2;base++){

          for(x=0;x<3;x++){
            for(i=0;i<3;i++){
              for(j=0;j<3;j++){
                for(k=0;k<3;k++){
                  DStepTransformationDRotation[step][0][base][x][i][j] += Pairframes[step+1][k][j]*DPairMidDRotation[step][base][x][k][i];
                  DStepTransformationDRotation[step][1][base][x][i][j] += DPairMidDRotation[step+1][base][x][k][j]*Pairframes[step][k][i];
                }
              }
            }
          }

          for(i=0;i<3;i++){
            for(x=0;x<3;x++){
              for(k=0;k<3;k++){
                for(l=0;l<3;l++){
                  DRotationVectorDStepRotation[step][0][base][x][i] += DStepTransformationDRotation[step][0][base][x][k][l]*DRotationVectorDStepTransformation[step][k][l][i];
                  DRotationVectorDStepRotation[step][1][base][x][i] += DStepTransformationDRotation[step][1][base][x][k][l]*DRotationVectorDStepTransformation[step][k][l][i];
                }
              }
            }
          }

          for(x=0;x<3;x++){
            for(i=0;i<3;i++){
              for(j=0;j<3;j++){
                for(k=0;k<3;k++){
                  DStepMidDRotation[step][1][base][x][i][j] += 0.5*HalfRotationMatrix[k][j]*DPairMidDRotation[step+1][base][x][i][k];
                  DStepMidDRotation[step][0][base][x][i][j] += 0.5*HalfRotationMatrix[j][k]*DPairMidDRotation[step][base][x][i][k];
                }
              }
            }
          }
        
        }
      }

     // calculate parameters
 
      for(pair=0;pair<npairs;pair++){
        buckle[pair] = dotProduct(RotationvectorPair[pair],Pairframes[pair][0]);
        propeller[pair] = dotProduct(RotationvectorPair[pair],Pairframes[pair][1]);
        opening[pair] = dotProduct(RotationvectorPair[pair],Pairframes[pair][2]);

        shear[pair] = dotProduct(pairTranslations[pair],Pairframes[pair][0]);
        stretch[pair] = dotProduct(pairTranslations[pair],Pairframes[pair][1]);
        stagger[pair] = dotProduct(pairTranslations[pair],Pairframes[pair][2]);
      }

      
      for(step=0;step<nsteps;step++){
        tilt[step] = dotProduct(RotationvectorStep[step],Stepframes[step][0]);
        roll[step] = dotProduct(RotationvectorStep[step],Stepframes[step][1]);
        twist[step] = dotProduct(RotationvectorStep[step],Stepframes[step][2]);

        shift[step] = dotProduct(stepTranslations[step],Stepframes[step][0]);
        slide[step] = dotProduct(stepTranslations[step],Stepframes[step][1]);
        rise[step] = dotProduct(stepTranslations[step],Stepframes[step][2]);
      }
      
     //
     // bookkeepderivatives and outputting them

      N=0;  //counter of the atoms within the step
      M=0;  //counter of atoms that are before the step

      for(step=0;step<nsteps;step++){      
      	valuestilt[step]=getPntrToComponent("tilt"+to_string(step));
      	valuesroll[step]=getPntrToComponent("roll"+to_string(step));
      	valuestwist[step]=getPntrToComponent("twist"+to_string(step));
        
      	valuesshift[step]=getPntrToComponent("shift"+to_string(step));
      	valuesslide[step]=getPntrToComponent("slide"+to_string(step));
      	valuesrise[step]=getPntrToComponent("rise"+to_string(step));

        valuestilt[step]->set(tilt[step]); 
        valuesroll[step]->set(roll[step]); 
        valuestwist[step]->set(twist[step]);
 
        valuesshift[step]->set(shift[step]);         
        valuesslide[step]->set(slide[step]); 
        valuesrise[step]->set(rise[step]);
      	
        N=M; //set atom counter equal to the first atom in this step
        
       //now the derivatives of the parameters of this step are defined for all atoms
        for(n=0;n<M;n++){//atoms before the step give derivative 0
          setAtomsDerivatives(valuestilt[step],N,Vector(0,0,0));  
          setAtomsDerivatives(valuesroll[step],N,Vector(0,0,0));
          setAtomsDerivatives(valuestwist[step],N,Vector(0,0,0));
        // derivatives of the translations are calculated much simpler to be in the direction of the appropriate midframe vector
          setAtomsDerivatives(valuesshift[step],N,Vector(0,0,0));
          setAtomsDerivatives(valuesslide[step],N,Vector(0,0,0));
          setAtomsDerivatives(valuesrise[step],N,Vector(0,0,0));       
        }
        
        for(base=0;base<4;base++){
          if(base<2){pair=0;}else{pair=1;}

          for(x=0;x<3;x++){
            DtiltDBaseRot[base][x] =  dotProduct( RotationvectorStep[step], DStepMidDRotation[step][pair][base%2][x][0]); 
            DtiltDBaseRot[base][x] += dotProduct( DRotationVectorDStepRotation[step][pair][base%2][x], Stepframes[step][0]);

            DrollDBaseRot[base][x] =  dotProduct( RotationvectorStep[step], DStepMidDRotation[step][pair][base%2][x][1]);
            DrollDBaseRot[base][x] += dotProduct( DRotationVectorDStepRotation[step][pair][base%2][x], Stepframes[step][1]);

            DtwistDBaseRot[base][x] =  dotProduct( RotationvectorStep[step], DStepMidDRotation[step][pair][base%2][x][2]);
            DtwistDBaseRot[base][x] += dotProduct( DRotationVectorDStepRotation[step][pair][base%2][x], Stepframes[step][2]);
          }

          //invert atom definition
          VectorGeneric<3> InertiaRedef = Vector(0.0,0.0,0.0);
          for(n=0;n<AtomsInBase[2*step+base];n++){
            for(x=0;x<3;x++){
              InertiaRedef[x] += AtomsMasses[2*step+base][n]*(pow(BaseRedefs[2*step+base][n][(x+1)%3],2) + pow(BaseRedefs[2*step+base][n][(x+2)%3],2));
            }
          }   

          for(x=0;x<3;x++){
            DtiltDbaseframe[base][x] = Baseframes[2*step+base][(x+2)%3]*DtiltDBaseRot[base][(x+1)%3]/InertiaRedef[(x+1)%3] - Baseframes[2*step+base][(x+1)%3]*DtiltDBaseRot[base][(x+2)%3]/InertiaRedef[(x+2)%3] ;
            DrollDbaseframe[base][x] = Baseframes[2*step+base][(x+2)%3]*DrollDBaseRot[base][(x+1)%3]/InertiaRedef[(x+1)%3] - Baseframes[2*step+base][(x+1)%3]*DrollDBaseRot[base][(x+2)%3]/InertiaRedef[(x+2)%3] ;
            DtwistDbaseframe[base][x] = Baseframes[2*step+base][(x+2)%3]*DtwistDBaseRot[base][(x+1)%3]/InertiaRedef[(x+1)%3] - Baseframes[2*step+base][(x+1)%3]*DtwistDBaseRot[base][(x+2)%3]/InertiaRedef[(x+2)%3] ;
          }


          for(n=0;n<AtomsInBase[2*step+base];n++){ // small n keeps track of the index within the base
            DtiltDatom = Vector(0,0,0);
            DrollDatom = Vector(0,0,0);
            DtwistDatom = Vector(0,0,0);
            if(FITTED or n<3){//IF fitting was not used the output will only depend on the first three atoms in each base

              for(x=0;x<3;x++){
                DtiltDatom += DtiltDbaseframe[base][x]*BaseRedefs[2*step+base][n][x]*AtomsMasses[2*step+base][n];
                DrollDatom += DrollDbaseframe[base][x]*BaseRedefs[2*step+base][n][x]*AtomsMasses[2*step+base][n];
                DtwistDatom += DtwistDbaseframe[base][x]*BaseRedefs[2*step+base][n][x]*AtomsMasses[2*step+base][n];
              }
              
              // now the derivatives of atom N can be outputted
              setAtomsDerivatives(valuestilt[step],N,DtiltDatom);
              setAtomsDerivatives(valuesroll[step],N,DrollDatom);
              setAtomsDerivatives(valuestwist[step],N,DtwistDatom);

            }else{
              //If no fitting was done derivatives past n==3 give 0
              setAtomsDerivatives(valuestilt[step],N,Vector(0,0,0));  
              setAtomsDerivatives(valuesroll[step],N,Vector(0,0,0));
              setAtomsDerivatives(valuestwist[step],N,Vector(0,0,0));
            }

            // derivatives of the translations are calculated much simpler to be in the direction of the appropriate midframe vector they point along the apropriate step frame vector and as trans = refTop-refBottom with a minus for base 1 and 2 and plus for 3 and 4
           
            if(base<2){
              setAtomsDerivatives(valuesshift[step],N,-1*Stepframes[step][0]);
              setAtomsDerivatives(valuesslide[step],N,-1*Stepframes[step][1]);
              setAtomsDerivatives(valuesrise[step],N,-1*Stepframes[step][2]);
              M++;  //increase the counter keeping track of the starting atom of the next step when going over the bottom two bases
            }else{
              setAtomsDerivatives(valuesshift[step],N,Stepframes[step][0]);
              setAtomsDerivatives(valuesslide[step],N,Stepframes[step][1]);
              setAtomsDerivatives(valuesrise[step],N,Stepframes[step][2]);
            }
            N++;  // large N is the overall index in the atom list
          }

        }


        //log.printf("DrollDatom of atom %d =  %f, %f, %f\n", N, DrollDatom[0], DrollDatom[1], DrollDatom[2]);  
       
        for(n=N;n<natoms;n++){//atoms after the step give derivative 0
          setAtomsDerivatives(valuestilt[step],n,Vector(0,0,0));  
          setAtomsDerivatives(valuesroll[step],n,Vector(0,0,0));
          setAtomsDerivatives(valuestwist[step],n,Vector(0,0,0));

        // derivatives of the translations are calculated much simpler to be in the direction of the appropriate midframe vector
          setAtomsDerivatives(valuesshift[step],n,Vector(0,0,0));
          setAtomsDerivatives(valuesslide[step],n,Vector(0,0,0));
          setAtomsDerivatives(valuesrise[step],n,Vector(0,0,0));       
        }
      
      }

      N=0;  //counter of the atoms within the pairs
      M=0;  //counter of the atoms that are before the pairs

      for(pair=0;pair<npairs;pair++){
        valuesbuckle[pair]=getPntrToComponent("buckle"+to_string(pair));
        valuespropeller[pair]=getPntrToComponent("propeller"+to_string(pair));
        valuesopening[pair]=getPntrToComponent("opening"+to_string(pair));

        valuesshear[pair]=getPntrToComponent("shear"+to_string(pair));
        valuesstretch[pair]=getPntrToComponent("stretch"+to_string(pair));
        valuesstagger[pair]=getPntrToComponent("stagger"+to_string(pair));


        valuesbuckle[pair]->set(buckle[pair]); 
        valuespropeller[pair]->set(propeller[pair]); 
        valuesopening[pair]->set(opening[pair]);
 
        valuesshear[pair]->set(shear[pair]); 
        valuesstretch[pair]->set(stretch[pair]); 
        valuesstagger[pair]->set(stagger[pair]);

        N=M;

        for(n=0;n<M;n++){//atoms before the step give derivative 0
          setAtomsDerivatives(valuesbuckle[pair],n,Vector(0,0,0));  
          setAtomsDerivatives(valuespropeller[pair],n,Vector(0,0,0));
          setAtomsDerivatives(valuesopening[pair],n,Vector(0,0,0));

          // derivatives of the translations are calculated much simpler to be in the direction of the appropriate midframe vector
          setAtomsDerivatives(valuesshear[pair],n,Vector(0,0,0));
          setAtomsDerivatives(valuesstretch[pair],n,Vector(0,0,0));
          setAtomsDerivatives(valuesstagger[pair],n,Vector(0,0,0));      
        }
       
        for(base=0;base<2;base++){

          for(x=0;x<3;x++){
            DbuckleDBaseRot[base][x] =  dotProduct( RotationvectorPair[pair], DPairMidDRotation[pair][base][x][0]); 
            DbuckleDBaseRot[base][x] += dotProduct( DRotationVectorDPairRotation[pair][base][x], Pairframes[pair][0]);

            DpropellerDBaseRot[base][x] =  dotProduct( RotationvectorPair[pair], DPairMidDRotation[pair][base][x][1]);
            DpropellerDBaseRot[base][x] += dotProduct( DRotationVectorDPairRotation[pair][base][x], Pairframes[pair][1]);

            DopeningDBaseRot[base][x] =  dotProduct( RotationvectorPair[pair], DPairMidDRotation[pair][base][x][2]);
            DopeningDBaseRot[base][x] += dotProduct( DRotationVectorDPairRotation[pair][base][x], Pairframes[pair][2]);
          }


          for(l=0;l<3;l++){ // derivative needs to be calculated to every vector of the frames
          // high depends on the highest bqse numbers(2step +2) and (2step+3), Right denotes the lower number in a pair (2step+0) and (2step+2)
            DbuckleDbaseframe[base][l] = crossProduct(DbuckleDBaseRot[base], Baseframes[2*pair + base][l]);
            DpropellerDbaseframe[base][l] = crossProduct(DpropellerDBaseRot[base], Baseframes[2*pair + base][l]);
            DopeningDbaseframe[base][l] = crossProduct(DopeningDBaseRot[base], Baseframes[2*pair + base][l]);
          }

          //invert atom definition
          for(n=0;n<AtomsInBase[2*pair+base];n++){ // small n keeps track of the index within the base
            if(FITTED or n<3){ //IF fitting was not used the output will omly depend on the first three atoms in each base
              for(k=0;k<3;k++){ // the outputted derivatives requested are gradients and thus vectors containing the derivative of the cv to every coordinate of the atoms position as a component.
                DbuckleDatom[k] = 0;  
                DpropellerDatom[k] = 0; 
                DopeningDatom[k] = 0;  
                for(l=0;l<3;l++){  // the atom positions have been defined in terms of the local frame and the ref point. The partial derivative of the frame vectors to the atom position can be found trough implicit derivation of that defining eqution as the coordinates of the atoms. applying the chain rule then gives the derivative of the CV to the atom position. as th cv depends on all three frame vectors this step involves summing over those vectors here represented by k.  
                  DbuckleDatom[k] += DbuckleDbaseframe[base][l][k]*BaseRedefs[2*pair+base][n][l] ;  
                  DpropellerDatom[k] += DpropellerDbaseframe[base][l][k]*BaseRedefs[2*pair+base][n][l] ;
                  DopeningDatom[k] += DopeningDbaseframe[base][l][k]*BaseRedefs[2*pair+base][n][l] ;
                }
              }
            
             // now the derivatives of atom N can be outputted
              setAtomsDerivatives(valuesbuckle[pair],N,DbuckleDatom);  
              setAtomsDerivatives(valuespropeller[pair],N,DpropellerDatom);
              setAtomsDerivatives(valuesopening[pair],N,DopeningDatom);
            }else{
              setAtomsDerivatives(valuesbuckle[pair],N,Vector(0,0,0));  
              setAtomsDerivatives(valuespropeller[pair],N,Vector(0,0,0));
              setAtomsDerivatives(valuesopening[pair],N,Vector(0,0,0));
            }

           // derivatives of the translations are calculated much simpler to be in the direction of the appropriate midframe vector
            setAtomsDerivatives(valuesshear[pair],N,pow(-1,base)*Pairframes[pair][0]);
            setAtomsDerivatives(valuesstretch[pair],N,pow(-1,base)*Pairframes[pair][1]);
            setAtomsDerivatives(valuesstagger[pair],N,pow(-1,base)*Pairframes[pair][2]);
           
            N++;  // large N is the overall index in the atom list
          }
        }
        
        for(n=N;n<natoms;n++){//atoms after the step give derivative 0
          setAtomsDerivatives(valuesbuckle[pair],n,Vector(0,0,0));  
          setAtomsDerivatives(valuespropeller[pair],n,Vector(0,0,0));
          setAtomsDerivatives(valuesopening[pair],n,Vector(0,0,0));

          // derivatives of the translations are calculated much simpler to be in the direction of the appropriate midframe vector
          setAtomsDerivatives(valuesshear[pair],n,Vector(0,0,0));
          setAtomsDerivatives(valuesstretch[pair],n,Vector(0,0,0));
          setAtomsDerivatives(valuesstagger[pair],n,Vector(0,0,0));       
        }
       
      }

     //
     //
     // end
      break;
    }
  }

}



} // colvar
} // PLMD



