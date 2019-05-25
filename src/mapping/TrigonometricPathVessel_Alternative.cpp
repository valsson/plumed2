/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2019 The plumed team
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
#include "TrigonometricPathVessel_Alternative.h"
#include "vesselbase/VesselRegister.h"
#include "tools/PDB.h"

namespace PLMD {
namespace mapping {

PLUMED_REGISTER_VESSEL(TrigonometricPathVessel_Alternative,"GPATH_ALT")

void TrigonometricPathVessel_Alternative::registerKeywords( Keywords& keys ) {
  StoreDataVessel::registerKeywords(keys);
}

void TrigonometricPathVessel_Alternative::reserveKeyword( Keywords& keys ) {
  keys.reserve("vessel","GPATH_ALT","calculate the position on the path using trigonometry");
  keys.add("optional","CALC_THRESHOLD","use a threshold on the calcuations of the distances");
  keys.addOutputComponent("gaspath","GPATH_ALT","the position on the path calculated using trigonometry");
  keys.addOutputComponent("gazpath","GPATH_ALT","the distance from the path calculated using trigonometry");
}

TrigonometricPathVessel_Alternative::TrigonometricPathVessel_Alternative( const vesselbase::VesselOptions& da ):
  StoreDataVessel(da),
  projdir(ReferenceConfigurationOptions("DIRECTION")),
  mydpack1( 1, getAction()->getNumberOfDerivatives() ),
  mydpack2( 1, getAction()->getNumberOfDerivatives() ),
  mydpack3( 1, getAction()->getNumberOfDerivatives() ),
  mypack1( 0, 0, mydpack1 ),
  mypack2( 0, 0, mydpack2 ),
  mypack3( 0, 0, mydpack3 ),
  first_time(true),
  task_prepare_thold(0)
{
  mymap=dynamic_cast<Mapping*>( getAction() );
  plumed_massert( mymap, "Trigonometric path vessel can only be used with mappings");
  // Retrieve the index of the property in the underlying mapping
  if( mymap->getNumberOfProperties()!=1 ) error("cannot use trigonometric paths when there are multiple properties");

  for(unsigned i=0; i<mymap->getFullNumberOfTasks(); ++i) {
    if( mymap->getTaskCode(i)!=mymap->getPositionInFullTaskList(i) ) error("mismatched tasks and codes");
  }
  mymap->addComponentWithDerivatives("gspath"); mymap->componentIsNotPeriodic("gspath");
  sp=mymap->copyOutput( mymap->getNumberOfComponents()-1 ); sp->resizeDerivatives( mymap->getNumberOfDerivatives() );
  mymap->addComponentWithDerivatives("gzpath"); mymap->componentIsNotPeriodic("gzpath");
  zp=mymap->copyOutput( mymap->getNumberOfComponents()-1 ); zp->resizeDerivatives( mymap->getNumberOfDerivatives() );

  mymap->parse("CALC_THRESHOLD",task_prepare_thold);
  if(task_prepare_thold>0) {
    mymap->log.printf("  using a threshold on the calcuations of the distances, only distances within %d of the previous node are calculated\n",task_prepare_thold);
  }

  // Check we have PCA
  ReferenceConfiguration* ref0=mymap->getReferenceConfiguration(0);
  for(unsigned i=0; i<mymap->getFullNumberOfTasks(); ++i) {
    if( !(mymap->getReferenceConfiguration(i))->pcaIsEnabledForThisReference() ) error("pca must be implemented in order to use trigometric path");
    if( ref0->getName()!=(mymap->getReferenceConfiguration(i))->getName() ) error("cannot use mixed metrics");
    if( mymap->getNumberOfAtoms()!=(mymap->getReferenceConfiguration(i))->getNumberOfReferencePositions() ) error("all frames must use the same set of atoms");
    if( mymap->getNumberOfArguments()!=(mymap->getReferenceConfiguration(i))->getNumberOfReferenceArguments() ) error("all frames must use the same set of arguments");
  }

  cargs.resize( mymap->getNumberOfArguments() ); std::vector<std::string> argument_names( mymap->getNumberOfArguments() );
  for(unsigned i=0; i<mymap->getNumberOfArguments(); ++i) argument_names[i] = (mymap->getPntrToArgument(i))->getName();
  PDB mypdb; mypdb.setAtomNumbers( mymap->getAbsoluteIndexes() ); mypdb.addBlockEnd( mymap->getAbsoluteIndexes().size() );
  if( argument_names.size()>0 ) mypdb.setArgumentNames( argument_names );
  projdir.read( mypdb );
  mypack1.resize( mymap->getNumberOfArguments(), mymap->getNumberOfAtoms() ); ref0->setupPCAStorage( mypack1 );
  mypack2.resize( mymap->getNumberOfArguments(), mymap->getNumberOfAtoms() ); ref0->setupPCAStorage( mypack2 );
  mypack3.resize( mymap->getNumberOfArguments(), mymap->getNumberOfAtoms() );
  for(unsigned i=0; i<mymap->getNumberOfAtoms(); ++i) { mypack1.setAtomIndex(i,i); mypack2.setAtomIndex(i,i); mypack3.setAtomIndex(i,i); }
  mypack1_stashd_atoms.resize( mymap->getNumberOfAtoms() ); mypack1_stashd_args.resize( mymap->getNumberOfArguments() );
}

std::string TrigonometricPathVessel_Alternative::description() {
  return "values gspath and gzpath contain the position on and distance from the path calculated using trigonometry";
}

void TrigonometricPathVessel_Alternative::resize() {
  StoreDataVessel::resize();
  if( getAction()->derivativesAreRequired() ) {
    unsigned nderivatives=getAction()->getNumberOfDerivatives();
    sp->resizeDerivatives( nderivatives ); zp->resizeDerivatives( nderivatives );
  }
}

void TrigonometricPathVessel_Alternative::finish( const std::vector<double>& buffer ) {
  double mindist1, mindist2, mindist3;
  // Store the data calculated during mpi loop
  StoreDataVessel::finish( buffer );
  // Get current value of all arguments
  for(unsigned i=0; i<cargs.size(); ++i) cargs[i]=mymap->getArgument(i);
  // Determine closest and second closest point to current position
  double lambda=mymap->getLambda();
  std::vector<double> dist1( getNumberOfComponents() ), dist2( getNumberOfComponents() ), dist3( getNumberOfComponents() );

  if( first_time ) { // check all nodes, find closest
    retrieveValueWithIndex( 0, false, dist1 ); mindist1 = dist1[0];
    if( lambda>0.0 ) mindist1=-std::log( mindist1 )  / lambda;
    iclose1=getStoreIndex(0);
    for(unsigned i=1; i<mymap->getFullNumberOfTasks(); ++i) {
      retrieveValueWithIndex( i, false, dist2 ); mindist2 = dist2[0];
      if( lambda>0.0 ) { mindist2 = -std::log( mindist2 ) / lambda; }
      if( mindist2 < mindist1 ) {
        mindist1 = mindist2;
        iclose1  = getStoreIndex(i);
      }
    }
    first_time = false;
  }
  else { // remember closest node from previous step
    iclose1 = iclose1_prev;
  }

  if( iclose1 == 0 || iclose1 == (mymap->getFullNumberOfTasks() -1) ) { // geometric path formula breaks down at end points - exception thrown
    plumed_error()<<"simulation too close to end of path, geometric formula breaks down - use tighter upper/lower walls on gspath\n";
  }

  for( ;; ) { // iterate until closest node is in the middle of three consecuitve nodes
    iclose2 = iclose1 -1;
    iclose3 = iclose1 +1;
    retrieveValueWithIndex( iclose1, false, dist1 ); mindist1 = dist1[0];
    retrieveValueWithIndex( iclose2, false, dist2 ); mindist2 = dist2[0];
    retrieveValueWithIndex( iclose3, false, dist3 ); mindist3 = dist3[0];
    if( lambda>0.0 ) {
      mindist1 = -std::log( mindist1 ) / lambda;
      mindist2 = -std::log( mindist2 ) / lambda;
      mindist3 = -std::log( mindist3 ) / lambda;
    }
    //fprintf( stderr, "***1***  iclose1 = %3d  d1 = %15.10g  iclose2 = %3d  d2 = %15.10g  iclose3 = %3d  d3 = %15.10g\n", iclose1, mindist1, iclose2, mindist2, iclose3, mindist3 );
    if( mindist1 < mindist2 ) {
      if( mindist1 < mindist3 ) { // got the correct triplet of nodes to compute gspath/gzpath
        if( mindist2 > mindist3 ) {
          iclose2 = iclose1 +1;
          iclose3 = iclose1 -1;
        }
        //fprintf( stderr, "***2***  iclose1 = %3d  iclose2 = %3d  iclose3 = %3d\n", iclose1, iclose2, iclose3 );
        break;
      } else { // new iclose1 is the lesser of iclose2 and iclose3
        if ( mindist2 < mindist3 ) { iclose1 = iclose2; }
        else                       { iclose1 = iclose3; }
        //fprintf( stderr, "***3***  iclose1 = %3d  iclose2 = %3d  iclose3 = %3d\n", iclose1, iclose2, iclose3 );
      }
    } else { // new iclose1 is the lesser of iclose2 and iclose3
      if ( mindist2 < mindist3 ) { iclose1 = iclose2; }
      else                       { iclose1 = iclose3; }
      //fprintf( stderr, "***4***  iclose1 = %3d  iclose2 = %3d  iclose3 = %3d\n", iclose1, iclose2, iclose3 );
    }
    if( iclose1 == 0 || iclose1 == (mymap->getFullNumberOfTasks()-1) ) { // geometric path formula breaks down at end points - exception thrown
      plumed_error()<<"simulation too close to end of path, geometric formula breaks down - use tighter upper/lower walls on gspath\n";
    }
    //fprintf( stderr, "***5***  iclose1 = %3d  iclose2 = %3d  iclose3 = %3d\n", iclose1, iclose2, iclose3 );
  }

  iclose1_prev = iclose1;

  // We now have to compute vectors connecting the three closest points to the
  // new point
  double v1v1 = (mymap->getReferenceConfiguration( iclose1 ))->calculate( mymap->getPositions(), mymap->getPbc(), mymap->getArguments(), mypack1, true );
  double v3v3 = (mymap->getReferenceConfiguration( iclose2 ))->calculate( mymap->getPositions(), mymap->getPbc(), mymap->getArguments(), mypack3, true );
  ReferenceConfiguration* conf2=mymap->getReferenceConfiguration( iclose3 );
  double v2v2 = (mymap->getReferenceConfiguration( iclose1 ))->calc( conf2->getReferencePositions(), mymap->getPbc(), mymap->getArguments(), conf2->getReferenceArguments(), mypack2, true );
  (mymap->getReferenceConfiguration( iclose1 ))->extractDisplacementVector( conf2->getReferencePositions(), mymap->getArguments(), conf2->getReferenceArguments(), false, projdir );

  // Stash derivatives of v1v1
  for(unsigned i=0; i<mymap->getNumberOfArguments(); ++i) mypack1_stashd_args[i]=mypack1.getArgumentDerivative(i);
  if( mymap->getNumberOfAtoms()>0 ) {
    ReferenceAtoms* at = dynamic_cast<ReferenceAtoms*>( mymap->getReferenceConfiguration( iclose1 ) );
    const std::vector<double> & displace( at->getDisplace() );
    for(unsigned i=0; i<mymap->getNumberOfAtoms(); ++i) {
      mypack1_stashd_atoms[i]=mypack1.getAtomDerivative(i); mypack1.getAtomsDisplacementVector()[i] /= displace[i];
    }
  }
  // Calculate the dot product of v1 with v2
  double v1v2 = (mymap->getReferenceConfiguration(iclose1))->projectDisplacementOnVector( projdir, mymap->getArguments(), cargs, mypack1 );

  // This computes s value
  double spacing = mymap->getPropertyValue( iclose1, (mymap->property.begin())->first ) - mymap->getPropertyValue( iclose2, (mymap->property.begin())->first );
  double root = sqrt( v1v2*v1v2 - v2v2 * ( v1v1 - v3v3) );
  dx = 0.5 * ( (root + v1v2) / v2v2 - 1.);
  double path_s = mymap->getPropertyValue(iclose1, (mymap->property.begin())->first ) + spacing * dx; sp->set( path_s );
  double fact = 0.25*spacing / v2v2;
  // Derivative of s wrt arguments
  for(unsigned i=0; i<mymap->getNumberOfArguments(); ++i) {
    sp->setDerivative( i, fact*( mypack2.getArgumentDerivative(i) + (v2v2 * (-mypack1_stashd_args[i] + mypack3.getArgumentDerivative(i))
                                 + v1v2*mypack2.getArgumentDerivative(i) )/root ) );
  }
  // Derivative of s wrt atoms
  unsigned narg=mymap->getNumberOfArguments(); Tensor vir; vir.zero(); fact = 0.5*spacing / v2v2;
  if( mymap->getNumberOfAtoms()>0 ) {
    for(unsigned i=0; i<mymap->getNumberOfAtoms(); ++i) {
      Vector ader = fact*(( v1v2*mypack1.getAtomDerivative(i) + 0.5*v2v2*(-mypack1_stashd_atoms[i] + mypack3.getAtomDerivative(i) ) )/root + mypack1.getAtomDerivative(i) );
      for(unsigned k=0; k<3; ++k) sp->setDerivative( narg+3*i+k, ader[k] );
      vir-=Tensor( mymap->getPosition(i), ader );
    }
    // Set the virial
    unsigned nbase=narg+3*mymap->getNumberOfAtoms();
    for(unsigned i=0; i<3; ++i) for(unsigned j=0; j<3; ++j) sp->setDerivative( nbase+3*i+j, vir(i,j) );
  }
  // Now compute z value
  //ReferenceConfiguration* conf2=mymap->getReferenceConfiguration( iclose1 );
  conf2=mymap->getReferenceConfiguration( iclose1 );
  double v4v4=(mymap->getReferenceConfiguration( iclose2 ))->calc( conf2->getReferencePositions(), mymap->getPbc(), mymap->getArguments(),
              conf2->getReferenceArguments(), mypack2, true );
  // Extract vector connecting frames
  (mymap->getReferenceConfiguration( iclose2 ))->extractDisplacementVector( conf2->getReferencePositions(), mymap->getArguments(),
      conf2->getReferenceArguments(), false, projdir );
  // Calculate projection of vector on line connnecting frames
  double proj = (mymap->getReferenceConfiguration(iclose1))->projectDisplacementOnVector( projdir, mymap->getArguments(), cargs, mypack1 );
  double path_z = v1v1 + dx*dx*v4v4 - 2*dx*proj;
  // Derivatives for z path
  path_z = sqrt(path_z); zp->set( path_z ); vir.zero();
  for(unsigned i=0; i<mymap->getNumberOfArguments(); ++i) zp->setDerivative( i, (mypack1_stashd_args[i] - 2*dx*mypack1.getArgumentDerivative(i))/(2.0*path_z) );
  // Derivative wrt atoms
  if( mymap->getNumberOfAtoms()>0 ) {
    for(unsigned i=0; i<mymap->getNumberOfAtoms(); ++i) {
      Vector dxder; for(unsigned k=0; k<3; ++k) dxder[k] = ( 2*v4v4*dx - 2*proj )*spacing*sp->getDerivative( narg + 3*i+k );
      Vector ader = ( mypack1_stashd_atoms[i] - 2.*dx*mypack1.getAtomDerivative(i) + dxder )/ (2.0*path_z);
      for(unsigned k=0; k<3; ++k) zp->setDerivative( narg+3*i+k, ader[k] );
      vir-=Tensor( mymap->getPosition(i), ader );
    }
    // Set the virial
    unsigned nbase=narg+3*mymap->getNumberOfAtoms();
    for(unsigned i=0; i<3; ++i) for(unsigned j=0; j<3; ++j) zp->setDerivative( nbase+3*i+j, vir(i,j) );
  }
  //fprintf( stderr, "path_s = %15.10g  path_z = %15.10g\n", path_s, path_z );
}

bool TrigonometricPathVessel_Alternative::applyForce( std::vector<double>& forces ) {
  std::vector<double> tmpforce( forces.size(), 0.0 );
  forces.assign(forces.size(),0.0); bool wasforced=false;
  if( sp->applyForce( tmpforce ) ) {
    wasforced=true;
    for(unsigned j=0; j<forces.size(); ++j) forces[j]+=tmpforce[j];
  }
  tmpforce.assign(forces.size(),0.0);
  if( zp->applyForce( tmpforce ) ) {
    wasforced=true;
    for(unsigned j=0; j<forces.size(); ++j) forces[j]+=tmpforce[j];
  }
  return wasforced;
}

void TrigonometricPathVessel_Alternative::prepare() {
  if(task_prepare_thold>0 && !first_time) {
    // as this will be run before ::finish above
    mymap->deactivateAllTasks();
    mymap->taskFlags[iclose1_prev] = 1;
    // fprintf( stderr, "iclose1_prev = %3d\n",iclose1_prev);
    for(unsigned int i=0; i<task_prepare_thold; i++) {
      if(iclose1_prev >= (i+1)) {
        mymap->taskFlags[iclose1_prev-(i+1)] = 1;
      }
      if(iclose1_prev+(i+1) <= (mymap->getFullNumberOfTasks()-1)) {
        mymap->taskFlags[iclose1_prev+(i+1)] = 1;
      }
    }
    mymap->lockContributors();
  }
}


}
}
