#include "mass_spring.hpp"

int main()
{
  MassSpringSystem<2> mss;
  mss.setGravity( {0,-9.81} );
  auto fA = mss.addFix( { { 0.0, 0.0 } } );
  auto mA = mss.addMass( { 1, { 1.0, 0.0 } } );
  mss.addSpring ( { 1, 10, { fA, mA } }  );

  auto mB = mss.addMass( { 1, { 2.0, 0.0 } } );
  mss.addSpring ( { 1, 20, { mA, mB } } );
  
  cout << "mss: " << endl << mss << endl;


  double tend = 10;
  double steps = 1000;
  
  Vector<> x(2*mss.masses().size());
  Vector<> dx(2*mss.masses().size());  
  Vector<> ddx(2*mss.masses().size());  

  auto mss_func = make_shared<MSS_Function<2>> (mss);
  auto mass = make_shared<IdentityFunction> (x.size());      

  mss.setState (x, dx, ddx);
  
  SolveODE_Newmark(tend, steps, x, dx,  mss_func, mass,
                   [](double t, VectorView<double> x) { cout << "t = " << t
                                                             << ", x = " << Vec<4>(x) << endl; });
}
