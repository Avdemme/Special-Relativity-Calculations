/*This code looks at physical properties of objects under special relativistic 
 * scenarios*/


#include "stdio.h"
#include <iostream>
#include "math.h"
#include <cmath>

using namespace std;

int main ( int argc, char** argv ) {

  double velocity, time, length, mass,SpeedOfLight,gamma,beta;
  double DilatedTime,ContractedLength,AdjustedMass,Momentum,KineticEnergy,TotalEnergy;
  double LowKineticEnergy,LowTotalEnergy;
  double imaginarygamma, imaginaryDilatedTime, imaginaryContractedLength, imaginaryAdjustedMass;
  double imaginaryKE, imaginaryMomentum,imaginaryTotalEnergySquared,imaginaryTotalEnergySquaredMagnitude;
  double imaginaryTotalEnergy;
  
  SpeedOfLight=299792458.;
  
  cout << " Enter a Velocity (meters/second) (assume a positive velocity as direction does not matter), " <<endl;
  cout << " a period of time (seconds)," << endl;
  cout << " a length (meters),"<< endl;
  cout << " and a mass (kilograms): " << endl;

  cin >> velocity >> time >> length >> mass;
  
  beta=(velocity/SpeedOfLight);
  gamma=1/sqrt(1-pow(beta,2));
 
  DilatedTime=time*gamma;
  ContractedLength=length/gamma;
  AdjustedMass=mass*gamma;
  Momentum=gamma*mass*velocity;
  KineticEnergy=mass*pow(SpeedOfLight,2)*gamma-mass*pow(SpeedOfLight,2);
  TotalEnergy=sqrt(pow(KineticEnergy,2)+pow(mass,2)*pow(SpeedOfLight,4));
 
  LowKineticEnergy=mass*pow(velocity,2)/2;
  LowTotalEnergy=LowKineticEnergy+mass*pow(SpeedOfLight,2);

  imaginarygamma= -1/sqrt((-1)*(1-pow(beta,2)));
  imaginaryDilatedTime=time*imaginarygamma;
  imaginaryContractedLength=-length/imaginarygamma;
  imaginaryAdjustedMass=mass*imaginarygamma;
  imaginaryMomentum=imaginarygamma*mass*velocity;
  imaginaryKE=imaginaryMomentum*SpeedOfLight;
  imaginaryTotalEnergySquared=-pow(imaginaryKE,2)+pow(mass,2)*pow(SpeedOfLight,4);
  imaginaryTotalEnergySquaredMagnitude=abs(imaginaryTotalEnergySquared);
  imaginaryTotalEnergy=sqrt(imaginaryTotalEnergySquaredMagnitude);

  if ((beta>=0.4) && (beta<1.)) {

    cout << "beta=" << "\t" << beta << endl;
    cout << "gamma=" << "\t" << gamma << endl;
    cout << "Time Dilation Amount=" << "\t"  << DilatedTime << "seconds/seconds" << endl;
    cout << "Contracted Length=" << "\t" << ContractedLength << "meters"<< endl;
    cout << "Adjusted Mass=" << "\t" <<  AdjustedMass << "kilograms" << endl;
    cout << "Kinetic energy=" << "\t" << KineticEnergy << "joules" << endl;
    cout << "Total energy=" << "\t" << TotalEnergy << "joules" << endl;

  }

  else if ((beta>=0.) && (beta<0.4)) {
    cout << "beta=" << "\t" << beta << endl;
    cout << "gamma=" << "\t" << gamma << endl;
    cout << "Time Dilation Amount=" << "\t"  << DilatedTime << "seconds/seconds" << endl;
    cout << "Contracted Length=" << "\t" << ContractedLength << "meters"<< endl;
    cout << "Adjusted Mass=" << "\t" <<  AdjustedMass << "kilograms" << endl;
    cout << "Kinetic energy=" << "\t" << LowKineticEnergy << "joules" << endl;
    cout << "Total energy=" << "\t" << LowTotalEnergy << "joules" << endl;

  }

  else if (beta<0.) {
    cout << "I told you to use a positive velocity value. Please go back and input your velocity as a positive number" << endl;
  }

  else if (beta == 1. ) {
    cout << "beta=" << "\t" << 1 << endl;
    cout << "gamma=" << "\t" << "Infinity"  << endl;
    cout << "Time Dilation Amount=" << "\t"  << "infinity"  << "seconds/seconds" << endl;
    cout << "Contracted Length=" << "\t" << "0" << "meters"<< endl;
    cout << "Adjusted Mass=" << "\t" <<  "infinity" << "kilograms" << endl;
    cout << "Kinetic energy=" << "\t" << "infinity"  << "joules" << endl;
    cout << "Total energy=" << "\t" << "infinity"  << "joules" << endl;
  }

  
  else if ( beta>1. ) {
    cout << "beta=" << "\t" << beta  << endl;
    cout << "gamma=" << "\t" << imaginarygamma << "i"  << endl;
    cout << "Time Dilation Amount=" << "\t"  << imaginaryDilatedTime << "i"  << "\t" << "seconds/seconds" << endl;
    cout << "Contracted Length=" << "\t" << imaginaryContractedLength << "i"<< "\t" << "meters"<< endl;
    cout << "Adjusted Mass=" << "\t" <<  imaginaryAdjustedMass <<"i" << "\t" << "kilograms" << endl;
    cout << "Kinetic energy=" << "\t" << "(" << imaginaryKE<< ")" <<  "i"<<"\t" << "joules" << endl;
    cout << "Total energy=" << "\t" << imaginaryTotalEnergy  << "joules" << endl;
    cout << "These results are non-physical. Tachyons are not real.";
    cout << "\t";
  }


  return 0;
}

 

