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

//+PLUMEDOC rbbNAnoder
/*
This file provides a template for if you want to introduce a new CV.

<!-----You should add a description of your CV here---->
This CV will calculate the rigid base coordinates of a DNA or RNA strand.
they come in 12N-6 components with N the number of base pairs analysed.
they are named here as in the literature:
(tilt, roll, twist, shift slide rise, buckle, proppeller, opening, shear, stretch, and stagger)
the first three are base pair step rotations, the folowing three base pair step translations, 
the next are the base pair rotations and the final three describe the base pair translations
this verion is for analysis alone and doesn't calculate derivatives

\par Examples

<!---You should put an example of how to use your CV here--->error: expected unqualified-id before ‘switch’

\plumedfile
# This should be a sample input.
rotations: rbbNA BASE1=1,2,3 BASE2=4,5,6 BASE3=7,8,9 BASE4=10,11,12
PRINT ARG=rotations.roll,rotations.tilt,rotations.twist STRIDE=100 FILE=COLVAR
#here if the order of atoms per base is N1(Y)/N9(R) C'1 C2(Y)/C4(R) where Y are the Pyrimidenes (T and C), R are the Purines (A and G)
\endplumedfile
<!---You should reference here the other actions used in this example--->
(see also \ref PRINT)

*/
//+ENDPLUMEDOC

class rbbNAnoder : public Colvar {
  bool components;
  bool pbc;
  bool RNA;
  bool FITTED;
  enum mode {pair=0, step=1, sequence=2, empty=3};
  mode INPUTMODE;
  enum BaseType {A = 'A', C='C', G='G', T='T', U='U'};
  vector<char> BASES;
  
public:
  explicit rbbNAnoder(const ActionOptions&);
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

  vector<VectorGeneric<3>> RotationMatrix(vector<VectorGeneric<3>> frame1, vector<VectorGeneric<3>> frame2);
  VectorGeneric<3> RotationAxis(vector<VectorGeneric<3>> Q);
  double RotationAngle(vector<VectorGeneric<3>> Q);
  VectorGeneric<3> RotationVector(vector<VectorGeneric<3>> Q);

  vector<VectorGeneric<3>> midFrame(vector<VectorGeneric<3>> frame1, VectorGeneric<3> RotationVector);

  vector<VectorGeneric<3>> BaseFrame(VectorGeneric<3> glycosidic, VectorGeneric<3> internal);
  VectorGeneric<3> RefPoint(VectorGeneric<3> glycosidic, VectorGeneric<3> normal, VectorGeneric<3> nitrogen);
};

PLUMED_REGISTER_ACTION(rbbNAnoder,"RBBNANODER")

void rbbNAnoder::registerKeywords(Keywords& keys) {
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


void rbbNAnoder::completeBase(char type, AtomNumber c1, vector<AtomNumber> *atoms, AtomNumber *next){
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


rbbNAnoder::rbbNAnoder(const ActionOptions&ao):
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

      addComponent("buckle"); componentIsNotPeriodic("buckle");
      addComponent("propeller"); componentIsNotPeriodic("propeller");
      addComponent("opening"); componentIsNotPeriodic("opening");
      
      addComponent("shear"); componentIsNotPeriodic("shear");
      addComponent("stretch"); componentIsNotPeriodic("stretch");
      addComponent("stagger"); componentIsNotPeriodic("stagger");
      
      log.printf("rbb-NA noder Basepair between base1 %d and base2 %d \n", atomsBase1[0].serial(),atomsBase2[0].serial());
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
      
      addComponent("tilt"); componentIsNotPeriodic("tilt");
      addComponent("roll"); componentIsNotPeriodic("roll");
      addComponent("twist"); componentIsNotPeriodic("twist");
      
      addComponent("shift"); componentIsNotPeriodic("shift");
      addComponent("slide"); componentIsNotPeriodic("slide");
      addComponent("rise"); componentIsNotPeriodic("rise");

      log.printf("rbb-NA noder Basepairstep between base1 %d, base2 %d, base3 %d and base4 %d\n", atomsBases[0][0].serial(), atomsBases[1][0].serial(), atomsBases[2][0].serial(), atomsBases[3][0].serial());
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
        addComponent("buckle"+to_string(j)); componentIsNotPeriodic("buckle"+to_string(j));
        addComponent("propeller"+to_string(j)); componentIsNotPeriodic("propeller"+to_string(j));
        addComponent("opening"+to_string(j)); componentIsNotPeriodic("opening"+to_string(j));
         
        addComponent("shear"+to_string(j)); componentIsNotPeriodic("shear"+to_string(j));
        addComponent("stretch"+to_string(j)); componentIsNotPeriodic("stretch"+to_string(j));
        addComponent("stagger"+to_string(j)); componentIsNotPeriodic("stagger"+to_string(j));
      
        if(!(j==0)){  // add length-1 bp step parameters as components
          addComponent("tilt"+to_string(j-1)); componentIsNotPeriodic("tilt"+to_string(j-1));
          addComponent("roll"+to_string(j-1)); componentIsNotPeriodic("roll"+to_string(j-1));
          addComponent("twist"+to_string(j-1)); componentIsNotPeriodic("twist"+to_string(j-1));
      
          addComponent("shift"+to_string(j-1)); componentIsNotPeriodic("shift"+to_string(j-1));
          addComponent("slide"+to_string(j-1)); componentIsNotPeriodic("slide"+to_string(j-1));
          addComponent("rise"+to_string(j-1)); componentIsNotPeriodic("rise"+to_string(j-1));
        }
      }     
      log.printf("rbb-NA noder between ends %d, %d, %d and %d\n", ends[0].serial(),ends[1].serial(),ends[2].serial(),ends[3].serial());
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

  log.printf("  version is 19.July.22 [fixed name]\n");

  //get the atoms, then calculation can start

  requestAtoms(atoms);
}


//fitting
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void rbbNAnoder::fit(char type, vector<VectorGeneric<3>> atoms, VectorGeneric<3> *fit1, VectorGeneric<3> *fit2, VectorGeneric<3> *fit3) {
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
  // Construct the Matrices T and Tt
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

  // Construct the Matrices TxTt and TtxT
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


double rbbNAnoder::Determinant3x3(vector< VectorGeneric<3> > Matrix){
  
  double D = Matrix[0][0]*Matrix[1][1]*Matrix[2][2] + Matrix[1][0]*Matrix[2][1]*Matrix[0][2] + Matrix[2][0]*Matrix[0][1]*Matrix[1][2];
  D -= Matrix[0][0]*Matrix[2][1]*Matrix[1][2] + Matrix[1][1]*Matrix[0][2]*Matrix[2][0] + Matrix[0][1]*Matrix[1][0]*Matrix[2][2];
  return D;
}


void rbbNAnoder::EigenValsAndVectorsSymmetric(vector<VectorGeneric<3>> Matrix, VectorGeneric<3>* EigenValues, vector<VectorGeneric<3>>* EigenVectors){
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

double rbbNAnoder::Norm(VectorGeneric<3> vec){
  double norm = vec.modulo();
  return norm;
}

VectorGeneric<3> rbbNAnoder::Normalize(VectorGeneric<3> vec){
  double norm = vec.modulo();
  VectorGeneric<3> normalized = vec/norm;
  return normalized;
}

VectorGeneric<3> rbbNAnoder::AxisRotation(VectorGeneric<3> v, VectorGeneric<3> k, double theta){
  VectorGeneric<3> vrot = v*cos(theta) + crossProduct(k,v)*sin(theta) + k*dotProduct(k,v)*(1-cos(theta));
  return vrot;
}


vector<VectorGeneric<3>> rbbNAnoder::RotationMatrix(vector<VectorGeneric<3>> frame1, vector<VectorGeneric<3>> frame2){

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

VectorGeneric<3> rbbNAnoder::RotationAxis(vector<VectorGeneric<3>> Q){

  VectorGeneric<3> axis;
  axis[0] = Q[1][2]-Q[2][1];
  axis[1] = Q[2][0]-Q[0][2];
  axis[2] = Q[0][1]-Q[1][0];

  return axis;
}

double rbbNAnoder::RotationAngle(vector<VectorGeneric<3>> Q){

  double cost=(Q[0][0]+Q[1][1]+Q[2][2]-1)/2;
  if(-1>cost){cost = -0.99999999999999;}
  if(1<cost){cost = 0.99999999999999;}
  double theta=acos(cost);

  return theta;
}

VectorGeneric<3> rbbNAnoder::RotationVector(vector<VectorGeneric<3>> Q){
  VectorGeneric<3> Axis, axis;
  double theta;

  Axis = RotationAxis(Q);
  axis = Normalize(Axis);
  theta = RotationAngle(Q);

  return theta*axis;
}

vector<VectorGeneric<3>> rbbNAnoder::midFrame(vector<VectorGeneric<3>> frame1, VectorGeneric<3> RotationVector){
  vector<VectorGeneric<3>> mid(3, Vector(0,0,0));
  double theta = RotationVector.modulo();
  VectorGeneric<3> axis = RotationVector/theta;

  for(int i=0; i<3; i++){
    mid[i]=AxisRotation(frame1[i], axis, theta/2); //altered last 13.15 replaced -theta by theta
    mid[i]=mid[i]/mid[i].modulo();
  }

  return mid;
}

vector<VectorGeneric<3>> rbbNAnoder::BaseFrame(VectorGeneric<3> glycosidic, VectorGeneric<3> internal){
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

VectorGeneric<3> rbbNAnoder::RefPoint(VectorGeneric<3> glycosidic, VectorGeneric<3> normal, VectorGeneric<3> nitrogen){
  VectorGeneric<3> displacement, position;

  displacement = AxisRotation( glycosidic/glycosidic.modulo(), normal, 2.4691172928);
  position = nitrogen + 0.4702*(displacement/displacement.modulo());
  return position;

}


//actual calculation
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void rbbNAnoder::calculate(){


 // pbc handling
  if(pbc){makeWhole();}

 // 

 // initialize basic counter and numbers
  int i, j, l, k, n, N, M, x;
  int npairs, nbases, nsteps, natoms;  // global numbers
  int purineAtoms, pyrimidineAtoms; // number of atoms in base dependeing on type
  int base;
 // initialize dummys to set size of other varibles
  VectorGeneric<3> Zero3 = Vector(0, 0, 0);
  vector<double> Zero10(10, 0.0);
  vector<VectorGeneric<3>> Zero2x3(2, Vector(0, 0, 0)), Zero3x3(3, Vector(0, 0, 0)), Zero10x3(10, Vector(0, 0, 0));
  vector<vector<VectorGeneric<3>>> Zero2x2x3(2, Zero2x3), Zero3x3x3(3, Zero3x3), Zero2x3x3(2, Zero3x3);
  vector<vector<vector<VectorGeneric<3>>>> Zero3x3x3x3(3, Zero3x3x3), Zero2x3x3x3(2, Zero3x3x3), Zero2x2x3x3(2, Zero2x3x3);
  vector<vector<vector<vector<VectorGeneric<3>>>>> Zero2x2x3x3x3(2, Zero2x3x3x3);

 //

 // get size of chain for initialization
  nbases = BASES.size();
  npairs = nbases/2;
  nsteps = npairs-1;
 //


 // check if Base completion is used or input is done giving the three relevent atoms alone
  bool BASECOMPLETION = (FITTED or (rbbNAnoder::INPUTMODE==sequence));
  if(BASECOMPLETION){
    purineAtoms = 10;
    pyrimidineAtoms = 7;
  }else{
    purineAtoms = 3;
    pyrimidineAtoms = 3;
  }
 //

 // initialize containers for information regarding the bases

  VectorGeneric<3> glycosidic, internal;  // important intermediate variables to get to the base frames

  vector<int> AtomsInBase(nbases, 0);
  vector<VectorGeneric<3>> BaseRefs(nbases, Vector(0, 0, 0)), BaseCenter(nbases, Vector(0, 0, 0));
  vector<vector<VectorGeneric<3>>> AtomsBases(nbases, Zero10x3), refAtoms(nbases, Zero3x3), Baseframes(nbases, Zero3x3);

 //

 // collect the atoms and all needed info on them
  N=0;
  for(base=0;base<(nbases-1);base=base+2){ // for all pairs in the DNA chain check the first base, second base should be complementary
    switch(BASES[base]){  //statements to load the required atoms depinding on the sequence provided

      case('A'): {
        AtomsInBase[base]=purineAtoms;
        AtomsBases[base].resize(purineAtoms);

        for (n = 0; n < purineAtoms; n++) {
          AtomsBases[base][n] = getPosition(N);
          N++;
        }

        AtomsInBase[base+1]=pyrimidineAtoms;
        AtomsBases[base+1].resize(pyrimidineAtoms);

        for (n = 0; n < pyrimidineAtoms; n++) {
          AtomsBases[base + 1][n] = getPosition(N);
          N++;
        }

        break;
      }

      case('C'): {
        AtomsInBase[base]=pyrimidineAtoms;
        AtomsBases[base].resize(pyrimidineAtoms);

        for (n = 0; n < pyrimidineAtoms; n++) {
          AtomsBases[base][n] = getPosition(N);
          N++;
        }

        AtomsInBase[base+1]=purineAtoms;
        AtomsBases[base+1].resize(purineAtoms);

        for (n = 0; n < purineAtoms; n++) {
          AtomsBases[base + 1][n] = getPosition(N);
          N++;
        }

        break;
      }

      case('G'): {
        AtomsInBase[base]=purineAtoms;
        AtomsBases[base].resize(purineAtoms);

        for (n = 0; n < purineAtoms; n++) {
          AtomsBases[base][n] = getPosition(N);
          N++;
        }

        AtomsInBase[base+1]=pyrimidineAtoms;
        AtomsBases[base+1].resize(pyrimidineAtoms);

        for (n = 0; n < pyrimidineAtoms; n++) {
          AtomsBases[base + 1][n] = getPosition(N);
          N++;
        }

        break;
      }

      case('T'): {
        AtomsInBase[base]=pyrimidineAtoms;
        AtomsBases[base].resize(pyrimidineAtoms);

        for (n = 0; n < pyrimidineAtoms; n++) {
          AtomsBases[base][n] = getPosition(N);
          N++;
        }

        AtomsInBase[base+1]=purineAtoms;
        AtomsBases[base+1].resize(purineAtoms);

        for (n = 0; n < purineAtoms; n++) {
          AtomsBases[base + 1][n] = getPosition(N);
          N++;
        }

        break;
      }

      case('U'): {
        AtomsInBase[base]=pyrimidineAtoms;
        AtomsBases[base].resize(pyrimidineAtoms);

        for (n = 0; n < pyrimidineAtoms; n++) {
          AtomsBases[base][n] = getPosition(N);
          N++;
        }

        AtomsInBase[base+1]=purineAtoms;
        AtomsBases[base+1].resize(purineAtoms);

        for (n = 0; n < purineAtoms; n++) {
          AtomsBases[base + 1][n] = getPosition(N);
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


    if((base%2) == 1){  // deal with the 3'->5', 5'->3' anti-symmetry making the tangents of opposing bases point in the same direction()
      Baseframes[base][2] = -1*Baseframes[base][2]; //making the tangents of opposing bases point in the same direction
      Baseframes[base][1] = -1*Baseframes[base][1]; //also making the e2 vectors (connecting backbones) not point against but "with" each other
    }

  }
 //    
      

  switch(rbbNAnoder::INPUTMODE){

 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 ///////////////////////////////////////////Pair//////////////////////////////////////////////////////////////////////////////////////// 
    case(pair):{
     //initialize

      // transformtion
      vector<VectorGeneric<3>> Pairframe = Zero3x3, TransformationMatrixPair = Zero3x3;
      VectorGeneric<3> RotationvectorPair = Vector(0, 0, 0), pairTranslation = Vector(0, 0, 0);

      
      // output vars
      double buckle, propeller, opening, shear, stretch, stagger;
      //output pointers
      Value* valuebuckle; Value* valuepropeller; Value* valueopening;
      Value* valueshear; Value* valuestretch; Value* valuestagger;

      //calculate

      TransformationMatrixPair = RotationMatrix(Baseframes[0], Baseframes[1]);  // the rotation matrix that transforms the left base to the right (given in the lab frame)
      RotationvectorPair = RotationVector(TransformationMatrixPair); // convert the rotation to axis angle representation
      Pairframe = midFrame(Baseframes[0], RotationvectorPair); // find the frame halfway between the two base frames
      pairTranslation = BaseRefs[1]-BaseRefs[0];  // Calculate the translation vector between the bases


      // final calculation

      buckle = dotProduct(RotationvectorPair,Pairframe[0]);
      propeller = dotProduct(RotationvectorPair,Pairframe[1]);
      opening = dotProduct(RotationvectorPair,Pairframe[2]);

      shear = dotProduct(pairTranslation,Pairframe[0]);
      stretch = dotProduct(pairTranslation,Pairframe[1]);
      stagger = dotProduct(pairTranslation,Pairframe[2]);


      // output handling
      valuebuckle=getPntrToComponent("buckle");
      valuepropeller=getPntrToComponent("propeller");
      valueopening=getPntrToComponent("opening");

      valueshear=getPntrToComponent("shear");
      valuestretch=getPntrToComponent("stretch");
      valuestagger=getPntrToComponent("stagger");


      valuebuckle->set(buckle); 
      valuepropeller->set(propeller); 
      valueopening->set(opening);
 
      valueshear->set(shear); 
      valuestretch->set(stretch); 
      valuestagger->set(stagger);

      break;
    }
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 ///////////////////////////////////////////Step////////////////////////////////////////////////////////////////////////////////////////      
    case(step):{  // to be redone (simplified) from Sequence mode
      //initialize
      int pair;  // counters

      double tilt, roll, twist, shift, slide, rise;  // output Steps
      // information related to bases
      vector<int> AtomsInPair(2,0);

      // information on Pairs and Steps
      vector<vector<VectorGeneric<3>>> Pairframes(2, Zero3x3), TransformationMatrixPair(2, Zero3x3);
      vector<VectorGeneric<3>> Stepframe(3, Vector(0, 0, 0)), TransformationMatrixStep(3, Vector(0, 0, 0)), PairCenter(2, Vector(0, 0, 0)), RotationvectorPair(2, Vector(0, 0, 0));
      VectorGeneric<3> RotationvectorStep = Vector(0, 0, 0) ; VectorGeneric<3> stepTranslations = Vector(0, 0, 0) ;

      
     // initialize pointers for output
      Value* valuetilt;Value* valueroll;Value* valuetwist;
      Value* valueshift;Value* valueslide;Value* valuerise;

     //calculate
      for(pair=0;pair<2;pair++){ //input of loop is two base frames, output are the rotation and translation vector between them and their mid frame, appropriate derivatives are also set
        TransformationMatrixPair[pair] = RotationMatrix(Baseframes[2*pair], Baseframes[2*pair+1]);  // the rotation matrix that transforms the left base to the right (given in the lab frame)
        RotationvectorPair[pair] = RotationVector(TransformationMatrixPair[pair]); // convert the rotation to axis angle representation
        Pairframes[pair] = midFrame(Baseframes[2*pair], RotationvectorPair[pair]); // find the frame halfway between the two base frames
      }
   
      TransformationMatrixStep = RotationMatrix(Pairframes[0], Pairframes[1]);
      RotationvectorStep = RotationVector(TransformationMatrixStep);
      Stepframe = midFrame(Pairframes[0], RotationvectorStep);
      stepTranslations = 0.5*(BaseRefs[2]+BaseRefs[3] - BaseRefs[0]-BaseRefs[1]);

     //final calculation

      tilt = dotProduct(RotationvectorStep,Stepframe[0]);
      roll = dotProduct(RotationvectorStep,Stepframe[1]);
      twist = dotProduct(RotationvectorStep,Stepframe[2]);

      shift = dotProduct(stepTranslations,Stepframe[0]);
      slide = dotProduct(stepTranslations,Stepframe[1]);
      rise = dotProduct(stepTranslations,Stepframe[2]);
    
     //output
      valuetilt=getPntrToComponent("tilt");
      valueroll=getPntrToComponent("roll");
      valuetwist=getPntrToComponent("twist");
        
      valueshift=getPntrToComponent("shift");
      valueslide=getPntrToComponent("slide");
      valuerise=getPntrToComponent("rise");

      valuetilt->set(tilt); 
      valueroll->set(roll); 
      valuetwist->set(twist);
 
      valueshift->set(shift);         
      valueslide->set(slide); 
      valuerise->set(rise);
      	
      break;
    }
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 /////////////////////////////////////////Sequence//////////////////////////////////////////////////////////////////////////////////////

    case(sequence):{
     //initialize
      //define needed variables

      int pair, base, step;  // counters

     //

     //set size of vector variables and initialize correctly

      // dummys to set size of other varibles

      vector<double> tilt(nsteps, 0), roll(nsteps, 0), twist(nsteps, 0), shift(nsteps, 0), slide(nsteps, 0), rise(nsteps, 0);  // output Steps
      vector<double> buckle(npairs, 0), propeller(npairs, 0), opening(npairs, 0), shear(npairs, 0), stretch(npairs, 0), stagger(npairs, 0); // output Pairs

      // information related to bases

      // information on Pairs and Steps
      vector<vector<VectorGeneric<3>>> Pairframes(npairs, Zero3x3), Stepframes(nsteps, Zero3x3), TransformationMatrixPair(npairs, Zero3x3), TransformationMatrixStep(nsteps, Zero3x3);
      vector<VectorGeneric<3>> RotationvectorPair(npairs, Vector(0, 0, 0)), pairTranslations(npairs, Vector(0, 0, 0)), RotationvectorStep(nsteps, Vector(0, 0, 0)), stepTranslations(nsteps, Vector(0, 0, 0));

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
      }
   
  
      for(step=0;step<nsteps;step++){ //input of loop is two pair frames, output is the rotation and translation vector between them and their mid frame, appropriate derivatives are also set
        TransformationMatrixStep[step] = RotationMatrix(Pairframes[step], Pairframes[step+1]);
        RotationvectorStep[step] = RotationVector(TransformationMatrixStep[step]);
        Stepframes[step] = midFrame(Pairframes[step], RotationvectorStep[step]);
        stepTranslations[step] = 0.5*(BaseRefs[2*step+2]+BaseRefs[2*step+3] - BaseRefs[2*step]-BaseRefs[2*step+1]);
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
     //   output handling

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
      }

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
       
      }

     //
     // end
      break;
    }
  }

}



} // colvar
} // PLMD



