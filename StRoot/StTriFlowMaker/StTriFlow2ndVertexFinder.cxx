#include "StTriFlow2ndVertexFinder.h"
#include "StTriFlowConstants.h"

ClassImp(StTriFlow2ndVertexFinder)

//------------------------------------------------------------------------------------------------------------------
StTriFlow2ndVertexFinder::StTriFlow2ndVertexFinder()
{
}

/*
StTriFlow2ndVertexFinder::~StTriFlow2ndVertexFinder()
{
}
*/
//------------------------------------------------------------------------------------------------------------------
// Alex's Secondary Vertex Finder
Int_t StTriFlow2ndVertexFinder::fDCA_Helix_Estimate(StPhysicalHelixD helixA, StPhysicalHelixD helixB, Float_t &pathA, Float_t &pathB, Float_t &dcaAB)
{

  // Calculates the 2D crossing point, calculates the corresponding 3D point and returns pathA and pathB

  Double_t x1 = helixA.xcenter();
  Double_t y1 = helixA.ycenter();
  Double_t x2 = helixB.xcenter();
  Double_t y2 = helixB.ycenter();
  Double_t c1 = helixA.curvature();
  Double_t c2 = helixB.curvature();
  Double_t r1 = 0.0;
  Double_t r2 = 0.0;
  if(c1 != 0 && c2 != 0)
  {
    r1 = 1.0/c1;
    r2 = 1.0/c2;
  }

  Double_t x1_c = 0.0;
  Double_t y1_c = 0.0;
  Double_t x2_c = 0.0;
  Double_t y2_c = 0.0;

  Int_t bCross_points = fCross_points_Circles(x1,y1,r1,x2,y2,r2,x1_c,y1_c,x2_c,y2_c);

  //cout << "bCross_points = " << bCross_points << ", xyr(1) = {" << x1 << ", " << y1 << ", " << r1
  //    << "}, xyr(2) = {"  << x2 << ", " << y2 << ", " << r2 << "}, p1 = {" << x1_c << ", " << y1_c << "}, p2 = {" << x2_c << ", " << y2_c << "}" << endl;

  if(bCross_points == 0) return 0;

  StThreeVectorF pointA,pointB,pointA1,pointB1,pointA2,pointB2;

  Double_t path_lengthA_c1,path_lengthA_c2,path_lengthB_c1,path_lengthB_c2;

  // first crossing point for helix A
  pair< double, double > path_lengthA = helixA.pathLength(sqrt(x1_c*x1_c+y1_c*y1_c));
  Double_t path_lengthA1 = path_lengthA.first;
  Double_t path_lengthA2 = path_lengthA.second;
  pointA1 = helixA.at(path_lengthA1);
  pointA2 = helixA.at(path_lengthA2);
  if( ((x1_c-pointA1.x())*(x1_c-pointA1.x()) + (y1_c-pointA1.y())*(y1_c-pointA1.y())) <
      ((x1_c-pointA2.x())*(x1_c-pointA2.x()) + (y1_c-pointA2.y())*(y1_c-pointA2.y())))
  {
    path_lengthA_c1 = path_lengthA1;
  }
  else
  {
    path_lengthA_c1 = path_lengthA2;
  }

  // second crossing point for helix A
  path_lengthA = helixA.pathLength(sqrt(x2_c*x2_c+y2_c*y2_c));
  path_lengthA1 = path_lengthA.first;
  path_lengthA2 = path_lengthA.second;
  pointA1 = helixA.at(path_lengthA1);
  pointA2 = helixA.at(path_lengthA2);
  if( ((x2_c-pointA1.x())*(x2_c-pointA1.x()) + (y2_c-pointA1.y())*(y2_c-pointA1.y())) <
      ((x2_c-pointA2.x())*(x2_c-pointA2.x()) + (y2_c-pointA2.y())*(y2_c-pointA2.y())))
  {
    path_lengthA_c2 = path_lengthA1;
  }
  else
  {
    path_lengthA_c2 = path_lengthA2;
  }

  // first crossing point for helix B
  pair< double, double > path_lengthB = helixB.pathLength(sqrt(x1_c*x1_c+y1_c*y1_c));
  Double_t path_lengthB1 = path_lengthB.first;
  Double_t path_lengthB2 = path_lengthB.second;
  pointB1 = helixB.at(path_lengthB1);
  pointB2 = helixB.at(path_lengthB2);
  if( ((x1_c-pointB1.x())*(x1_c-pointB1.x()) + (y1_c-pointB1.y())*(y1_c-pointB1.y())) <
      ((x1_c-pointB2.x())*(x1_c-pointB2.x()) + (y1_c-pointB2.y())*(y1_c-pointB2.y())))
  {
    path_lengthB_c1 = path_lengthB1;
  }
  else
  {
    path_lengthB_c1 = path_lengthB2;
  }

  // second crossing point for helix B
  path_lengthB = helixB.pathLength(sqrt(x2_c*x2_c+y2_c*y2_c));
  path_lengthB1 = path_lengthB.first;
  path_lengthB2 = path_lengthB.second;
  pointB1 = helixB.at(path_lengthB1);
  pointB2 = helixB.at(path_lengthB2);
  if( ((x2_c-pointB1.x())*(x2_c-pointB1.x()) + (y2_c-pointB1.y())*(y2_c-pointB1.y())) <
      ((x2_c-pointB2.x())*(x2_c-pointB2.x()) + (y2_c-pointB2.y())*(y2_c-pointB2.y())))
  {
    path_lengthB_c2 = path_lengthB1;
  }
  else
  {
    path_lengthB_c2 = path_lengthB2;
  }

  pointA1 = helixA.at(path_lengthA_c1);
  pointA2 = helixA.at(path_lengthA_c2);

  pointB1 = helixB.at(path_lengthB_c1);
  pointB2 = helixB.at(path_lengthB_c2);

  if((pointA1-pointB1).mag() < (pointA2-pointB2).mag())
  {
    pathA = path_lengthA_c1;
    pathB = path_lengthB_c1;
    dcaAB = (pointA1-pointB1).mag();
  }
  else
  {
    pathA = path_lengthA_c2;
    pathB = path_lengthB_c2;
    dcaAB = (pointA2-pointB2).mag();
  }

  //cout << "3D-point (A1) = (" << pointA1.x() << ", " << pointA1.y() << ", " << pointA1.z() << "), calc = (" << x1_c << ", " << y1_c << ")" << endl;
  //cout << "3D-point (A2) = (" << pointA2.x() << ", " << pointA2.y() << ", " << pointA2.z() << "), calc = (" << x2_c << ", " << y2_c << ")" << endl;
  //cout << "3D-point (B1) = (" << pointB1.x() << ", " << pointB1.y() << ", " << pointB1.z() << "), calc = (" << x1_c << ", " << y1_c << ")" << endl;
  //cout << "3D-point (B2) = (" << pointB2.x() << ", " << pointB2.y() << ", " << pointB2.z() << "), calc = (" << x2_c << ", " << y2_c << ")" << endl;

  /*
     StThreeVectorF pointA,pointB,pointA1,pointB1,pointA2,pointB2;
     pointA = helixA.at(0); // 3D-vector of helixA at path 0
     pointB = helixB.at(0); // 3D-vector of helixA at path 0
     Double_t Delta_phiA = TMath::ATan2(pointA.y()-y1,pointA.x()-x1); // polar angle for path length 0 -> starting point for helix A
     Double_t Delta_phiB = TMath::ATan2(pointB.y()-y2,pointB.x()-x2); // polar angle for path length 0 -> starting point for helix B

  //Double_t phi_valA = 0.0;
  //Double_t x_new = x1 + r1*TMath::Cos(Delta_phi+phi_valA);
  //Double_t y_new = y1 + r1*TMath::Sin(Delta_phi+phi_valA);
  //cout << "3D-point = (" << pointA.x() << ", " << pointA.y() << "), calc = (" << x_new << ", " << y_new << ")" << endl;

  Double_t path_lengthA1 = (TMath::ATan2((y1_c-y1),(x1_c-x1)) - Delta_phiA)*r1;  // crossing point 1 for helix A
  Double_t path_lengthA2 = (TMath::ATan2((y2_c-y1),(x2_c-x1)) - Delta_phiA)*r1;  // crossing point 2 for helix A

  Double_t path_lengthB1 = (TMath::ATan2((y1_c-y2),(x1_c-x2)) - Delta_phiB)*r2;  // crossing point 1 for helix B
  Double_t path_lengthB2 = (TMath::ATan2((y2_c-y2),(x2_c-x2)) - Delta_phiB)*r2;  // crossing point 2 for helix B

  pointA1 = helixA.at(path_lengthA1);
  pointA2 = helixA.at(path_lengthA2);

  pointB1 = helixB.at(path_lengthB1);
  pointB2 = helixB.at(path_lengthB2);

  Double_t path_lengthTest = (TMath::ATan2((pointA.y()-y1),(pointA.x()-x1)) - Delta_phiA)*r1;  // crossing point 1 for helix A
  cout << "path_lengthTest = " << path_lengthTest << endl;


  cout << "3D-point (A1) = (" << pointA1.x() << ", " << pointA1.y() << ", " << pointA1.z() << "), calc = (" << x1_c << ", " << y1_c << ")" << endl;
  cout << "3D-point (A2) = (" << pointA2.x() << ", " << pointA2.y() << ", " << pointA2.z() << "), calc = (" << x2_c << ", " << y2_c << ")" << endl;
  cout << "3D-point (B1) = (" << pointB1.x() << ", " << pointB1.y() << ", " << pointB1.z() << "), calc = (" << x1_c << ", " << y1_c << ")" << endl;
  cout << "3D-point (B2) = (" << pointB2.x() << ", " << pointB2.y() << ", " << pointB2.z() << "), calc = (" << x2_c << ", " << y2_c << ")" << endl;
  */

  return 1;

}

Int_t StTriFlow2ndVertexFinder::fCross_points_Circles(Double_t x1, Double_t y1, Double_t r1, Double_t x2, Double_t y2, Double_t r2, Double_t &x1_c, Double_t &y1_c, Double_t &x2_c, Double_t &y2_c)
{
    // (x1,y1) -> center of circle 1, r1 -> radius of circle 1
    // (x2,y2) -> center of circle 2, r2 -> radius of circle 2
    // (x1_c,y1_c) -> crossing point 1
    // (x2_c,y2_c) -> crossing point 2
    // Solution see Wikipedia

    Double_t s = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));  // distance between the center of the circles

    x1_c = 0;
    y1_c = 0;
    x2_c = 0;
    y2_c = 0;

    if(x1 != x2 && y1 != y2 && s < (r1+r2))
    {
        Double_t m  = (x1-x2)/(y2-y1);
        Double_t n  = (-r2*r2 + r1*r1 + y2*y2 - y1*y1 + x2*x2 - x1*x1)/(2.0*(y2-y1));
        Double_t p  = (2.0*(-x1 + m*(n-y1)))/(1.0 + m*m);
        Double_t q  = (x1*x1 + (n-y1)*(n-y1) -r1*r1)/(1.0 + m*m);
        Double_t p2 = (p/2.0)*(p/2.0);

        if(p2 >= q)
        {
            x1_c = (-p/2.0) + sqrt(p2 - q);
            x2_c = (-p/2.0) - sqrt(p2 - q);
            y1_c = m*x1_c + n;
            y2_c = m*x2_c + n;
            return 1;
        }
        else return 0;
    }
    else return 0;

}

StThreeVectorF StTriFlow2ndVertexFinder::calcVertexAnalytical(StThreeVectorF &base1, StThreeVectorF &dir1, StThreeVectorF &base2, StThreeVectorF &dir2)
{
  // Calculates the Vertex of two straight lines define by the vectors base and dir
  //
  //      g: x1 = base1 + l * dir1 
  //      h: x2 = base2 + m * dir2 , where l,m are real numbers 
  //
  // 1. are g and h
  //       parallel / identical, i.e. are dir1 and dir2 linear dependent?
  //       
  //                                        /-                               
  //                                        |
  //                                        |   = 0    linear dependent, no unique solution, returning dummy  
  //      => cross product : dir1 x dir2 = -|  
  //                                        |  != 0    linear independent
  //                                        |
  //                                        \\-         
  //
  // 2. are g and h 
  //       skew or do they have a crossing point, i.e are dir1, dir2 and (base1 - base2) linear dependent ?
  //
  //                                                    /-                               
  //                                                    |
  //                                                    |   = 0    linear dependent
  //                                                    |          g and h are intersecting
  //                                                    |          calculating vertex as point of intersection
  //                                                    |
  //    => determinant: det[ dir1, dir2, base1-base2]= -|
  //                                                    |  != 0    linear independent
  //                                                    |          g and h are skew
  //                                                    |          calulating vertex as point of closest approach
  //                                                    |
  //                                                    \\-         
  //  
  // 3.
  //    (a) calculating intersection point
  //    (b) calculating point of closest approach



  // 1. exists a unique solution ?

  if ((dir1.cross(dir2)).mag()> 0.) // dir1 and dir2 linear independent
    {
      // straight lines are either skew or have a cross point

      StThreeVectorF diff = base1;
      diff-=base2; // Difference of two base vectors base1 - base2
      
      // 2. skew or intersecting ?
	
      if (fabs(calcDeterminant(dir2, dir1 ,diff))>0.) 
	{
	  // 3. (b) skew 
	  return StThreeVectorF(calculatePointOfClosestApproach(base1, dir1, base2, dir2));
	}
      else
	{
	  // 3. (a) intersection 
	  return StThreeVectorF(calculateCrossPoint(base1 ,dir1, base2 ,dir2));
	}
    }
  else
    {
      // dir1 and dir2 linear dependent -> g1 and g2 identical or parallel
      return StThreeVectorF(-10000000.,-10000000.,-10000000.);
    }
  return StThreeVectorF(-10000000.,-10000000.,-10000000.);
}

StThreeVectorF StTriFlow2ndVertexFinder::calculatePointOfClosestApproach(StThreeVectorF &base1, StThreeVectorF &dir1, StThreeVectorF &base2, StThreeVectorF &dir2)
{
  //  calculating point of closest approach
  //        
  //        from the equations of the straight lines of g and h 
  //        g: x1 = base1 + l * dir1 
  //        h: x2 = base2 + m * dir2 
  //        
  //        you can construct the following planes:
  //        
  //        E1: e1 = base1  +  a * dir1  +  b * (dir1 x dir2)
  //        E2: e2 = base2  +  s * dir2  +  t * (dir1 x dir2)
  //        
  //        now the intersection point of E1 with g2 = {P1} 
  //        and the intersection point of E2 with g1 = {P2}
  //        
  //        form the base points of the perpendicular to both straight lines.
  //        
  //        The point of closest approach is the middle point between P1 and P2: 
  //        
  //        vertex = (p2 - p1)/2
  // 
  //        E1 ^ g2:
  //
  //           e1 = x2
  //    -->    base1  +  a * dir1  +  b * (dir1 x dir2) = base2 + m * dir2 
  //    -->    base1 - base2 = m * dir2  -  a * dir1  -  b * (dir1 x dir2)       
  //                                          (m)
  //    -->    [ dir2, -dir1, -(dir1 x dir2)] (a) = base1 - base2        
  //                                          (b)
  //           
  //           using CRAMER's RULE you can find the solution for m (a,b, not used)
  //           
  //           using the rules for converting determinants:
  //           
  //           D12 = det [dir2, -dir1, -(dir1 x dir2)] 
  //               = det [dir2,  dir1,  (dir1 x dir2)]
  //           
  //           Dm  = det [base1 - base2, -dir1, -(dir1 x dir2)]
  //               = det [base1 - base2,  dir1,  (dir1 x dir2)]
  //  
  //            m  = Dm/D12
  //           
  //           P1: p1 = x2(m)
  //                  = base2 + Dm/D12 * dir2
  //
  //        E2 ^ g1:
  //
  //           e2 = x1
  //    -->    base2  +  s * dir2  +  t * (dir1 x dir2) = base1 + l * dir1 
  //    -->    base2 - base1 = l * dir1  -  s * dir2  -  t * (dir1 x dir2)       
  //                                          (l)
  //    -->    [ dir1, -dir2, -(dir1 x dir2)] (s) = base2 - base1        
  //                                          (t)
  //           
  //           again using CRAMER's RULE you can find the solution for m (a,b, not used)
  //           
  //           using the rules for converting determinants:
  //           
  //           D21 =  det [dir1, -dir2, -(dir1 x dir2)] 
  //               =  det [dir1,  dir2,  (dir1 x dir2)]
  //               = -det [dir2,  dir1,  (dir1 x dir2)]
  //               = -D12
  //           
  //           Dl  =  det [base2 - base1, -dir2, -(dir1 x dir2)]
  //               =  det [base2 - base1,  dir1,  (dir1 x dir2)]
  //               = -det [base1 - base2,  dir1,  (dir1 x dir2)]
  //
  //            l  =   Dl/D21
  //               = - Dl/D12
  //           
  //           P2: p2 = x1(m)
  //                  = base1 - Dl/D12 * dir1
  //           
  //           
  //           vertex = p1 + 1/2 * (p2 - p1)
  //                  = 1/2 * (p2 + p1)
  //                  = 1/2 *( (base1 + base2) +  1/D12 * ( Dm * dir2 - Dl * dir1) )
  //                      

  StThreeVectorF cross = dir1.cross(dir2); // cross product: dir1 x dir2

  // straight lines are either skew or have a cross point
	      
  StThreeVectorF diff = base1;
  diff-=base2; // Difference of two base vectors base1 - base2
		
  Double_t D;
  D =  calcDeterminant(dir2, dir1 ,cross);

  if (!(fabs(D) > 0.))
    {
      ::Warning(":calculatePointOfClosestApproach","Dirs and cross-product are lin. dependent: returning default Vertex (-20000,-20000,-20000)");
      return StThreeVectorF(-20000.,-20000.,-20000.);
    }

  Double_t Dm =  calcDeterminant(diff , dir1, cross);
  Double_t Dl = -calcDeterminant(diff , dir2, cross);

  StThreeVectorF vertex;
  StThreeVectorF dm;
  StThreeVectorF dl;

  dm = dir2;
  dm *= Dm;

  dl = dir1;
  dl *= Dl;

  vertex = dm - dl;

  vertex *= ((1.)/D);

  vertex+=base1;
  vertex+=base2;
  vertex*=0.5;

  return StThreeVectorF(vertex);
}

StThreeVectorF StTriFlow2ndVertexFinder::calculateCrossPoint(StThreeVectorF &base1, StThreeVectorF &dir1, StThreeVectorF &base2, StThreeVectorF &dir2)
{ 
  // calculating cross point 
  // taking all three equations into account solving the overdetermined set of lin. equations
  // of 
  // base1 + l * dir2 =  base1 + m * dir2 
  //
  // set of lin. equations:
  //  
  //   base1(0) + l * dir1(0) = base2(0) + m * dir2(0) 
  //   base1(1) + l * dir1(1) = base2(1) + m * dir2(1)
  //   base1(2) + l * dir1(2) = base2(2) + m * dir2(2) this line is ignored
  //
  //   written in matrix form
  //
  //        l
  //   M * |   | = base2 - base1
  //       \\ m /
  //
  //   M is a 3x2 matrix
  //     
  // to solve multiply the equation by the transposed Matrix of M from the left: M 
  //     
  //  T      /  l \\                                                               .
  // M * M * |    | = M  * (base2 - base1)
  //         \\ -m /
  // MIND THE '-' of m
  //
  //     / dir1(0) dir2(0) \\                                                      .
  //     |                 |    T   / dir1(0) dir1(1) dir1(2) \\                   .
  // M = | dir1(1) dir2(1) |,  M  = |                         |
  //     |                 |        \\ dir2(0) dir2(1) dir2(2) /                   .
  //     \\ dir1(2) dir2(2) /                                    
  //
  //  T      / (dir1(0)*dir1(0) + dir1(1)*dir1(1) + dir1(2)*dir1(2))   (dir1(0)*dir2(0) + dir1(1)*dir2(1) + dir1(2)*dir2(2))  \\ .
  // M * M = |                                                                                                                |
  //         \\ (dir1(0)*dir2(0) + dir1(1)*dir2(1) + dir1(2)*dir2(2))   (dir2(0)*dir2(0) + dir2(1)*dir2(1) + dir2(2)*dir2(2))  /                        
  //
  //  T       / d1d1 d1d2 \\                           .
  // M  * M = |           |
  //          \\ d1d2 d2d2 /
  //
  // diff = base2 - base1
  //
  //  T           /  (dir1(0)*diff(0) + dir1(1)*diff(1) + dir1(2)*diff(2)) \\         .
  // M  * diff =  |                                                        |
  //              \\  (dir2(0)*diff(0) + dir2(1)*diff(1) + dir2(2)*diff(2)) /
  //
  //  T           /  d1diff  \\                                          .
  // M  * diff =  |          |
  //              \\  d2diff  /
  // 
  // now the new Matrix set is to be solved by CRAMER'S Rule:
  // 
  // / d1d1 d1d2 \\   /  l \\   /  d1diff \\                   .
  // |           | * |    | = |          |
  // \\ d1d2 d2d2 /   \\ -m /   \\  d2diff /
  //
  //     | d1d1 d1d2 |
  // D = |           | = d1d1*d2d2 - d1d2*d1d2;
  //     | d1d2 d2d2 |
  // 
  //     | d1diff d1d2 |
  // Dl= |              | = d1diff*d2d2 - d1d2*d2diff;
  //     | d2diff d2d2 |              
  //
  // l = Dl/D = l_cross
  // 
  // vertex = base1 + l_cross * dir1
  //

  Double_t d1d1 = dir1(0)*dir1(0) + dir1(1)*dir1(1) + dir1(2)*dir1(2);
  Double_t d2d2 = dir2(0)*dir2(0) + dir2(1)*dir2(1) + dir2(2)*dir2(2);
  Double_t d1d2 = dir1(0)*dir2(0) + dir1(1)*dir2(1) + dir1(2)*dir2(2);
  
  Double_t D = d1d1*d2d2 - (d1d2*d1d2);
  
  if (!(fabs(D) > 0.))
    {
      ::Warning("calculateCrossPoint","Error while calculating cross point ... eqns are lin. dependent:returning default Vertex (-20000,-20000,-20000)");
      return StThreeVectorF(-20000.,-20000.,-20000.);
    }

  Double_t d1diff = dir1(0)*(base2(0)-base1(0))+dir1(1)*(base2(1)-base1(1))+dir1(2)*(base2(2)-base1(2));
  Double_t d2diff = dir2(0)*(base2(0)-base1(0))+dir2(1)*(base2(1)-base1(1))+dir2(2)*(base2(2)-base1(2));

  Double_t Dlambda = d1diff*d2d2-d1d2*d2diff;
  
  Double_t lambda = Dlambda/D;
  
  StThreeVectorF vertex;
  vertex += dir1;
  vertex *= lambda;
  vertex += base1;

  //cout << "Cross point calculated" << endl;
  return StThreeVectorF(vertex);

 // return StThreeVectorF(-20000.,-20000.,-20000.);
}

Double_t StTriFlow2ndVertexFinder::calcDeterminant(StThreeVectorF& v1,StThreeVectorF& v2,StThreeVectorF& v3)
{
  // calculating the Determinant of a 3 x 3 Matrix 
  // with the column vectors [v1, v2, v3]
  // using the RULE of SARRUS
  //
  // | v1(0)   v2(0)   v3(0) |      | v1(0) v2(0) v3(0)|v1(0) v2(0)  .
  // |                       |      |  \\     \\     X   |  /     /    . 
  // |                       |      |   \\     \\   / \\  | /     /     . 
  // |                       |      |    \\     \\ /   \\ |/     /      . 
  // |                       |      |     \\     X     \\/     /       . 
  // |                       |      |      \\   / \\    /\\    /        .  
  // |                       |      |       \\ /   \\  / |\\  /         . 
  // | v1(1)   v2(1)   v3(1) |   =  | v1(1) v2(1) v3(1)|v1(1) v2(1)  .
  // |                       |      |       / \\    /\\  | /\\          . 
  // |                       |      |      /   \\  /  \\ |/  \\         . 
  // |                       |      |     /     \\/    \\/    \\        . 
  // |                       |      |    /      /\\    /\\     \\       . 
  // |                       |      |   /      /  \\  / |\\     \\      .  
  // |                       |      |  /      /    \\/  | \\     \\     . 
  // | v1(2)   v2(2)   v3(2) |      | v1(2) v2(2) v3(2)| v1(2) v2(2) .  
  //                                 /      /     /  \\     \\     \\   .
  //                                                                
  //                                -      -     -    +     +     +  .

  return ( v1(0) * v2(1) * v3(2) 
	   + v2(0) * v3(1) * v1(2) 
	   + v3(0) * v1(1) * v2(2) 
	   - v3(0) * v2(1) * v1(2) 
	   - v1(0) * v3(1) * v2(2) 
	   - v2(0) * v1(1) * v3(2)); 
}

Double_t StTriFlow2ndVertexFinder::calculateMinimumDistance(StThreeVectorF &base1, StThreeVectorF &dir1, StThreeVectorF &base2, StThreeVectorF &dir2)
{
  // calculates the minimum distance of two tracks given as parametric straights x = base + n * dir

  StThreeVectorF cross = dir1.cross(dir2);

  StThreeVectorF ab = base1 - base2;

  if ( !( fabs(cross.mag())>0.)) // dir1 || dir2
    {
      return calculateMinimumDistanceStraightToPoint(base1, dir1, base2);
    }
 
  return fabs(ab.dot(cross)/cross.mag());
}

Double_t StTriFlow2ndVertexFinder::calculateMinimumDistanceStraightToPoint(StThreeVectorF &base, StThreeVectorF &dir, StThreeVectorF &point)
{
  // calculates the minimum distance of a point to a straight given as parametric straight x = base + n * dir

  if (!(dir.mag()>0))
    {
      return -1000000.;
    }
  
  StThreeVectorF diff = base-point;

  StThreeVectorF cross = dir.cross(diff);
  
  return cross.mag()/dir.mag();
}

void StTriFlow2ndVertexFinder::fHelixABdca_start_params(StPhysicalHelixD helixA, StPhysicalHelixD helixB, Float_t &pathA, Float_t &pathB, Float_t &dcaAB, Float_t path_in_A, Float_t path_in_B)
{
    Float_t pA[2] = {path_in_A+2,path_in_A-2};
    Float_t pB[2] = {path_in_B+2,path_in_B-2}; // the two start values for pathB, 0.0 is the origin of the helix at the first measured point
    Float_t distarray[2];
    StThreeVectorF testA, testB;
    for(Int_t r = 0; r < 2; r++)
    {
        testB     = helixB.at(pB[r]);  // 3D-vector of helixB point at path pB[r]

        Float_t pathA_dca = -999.0;
        Float_t dcaAB_dca = -999.0;
        fHelixAtoPointdca_start_params(testB,helixA,pathA_dca,dcaAB_dca,path_in_A); // new helix to point dca calculation
        testA = helixA.at(pathA_dca);
        //testA     = helixA.at(helixA.pathLength(testB)); // 3D-vector of helixA point at dca to testB

        distarray[r] = (testA-testB).mag(); // dca between helixA and helixB
    }
    Int_t loopcounter = 0;
    Float_t scale = 1.0;
    Float_t flip  = 1.0; // checks if the minimization direction changed
    Float_t scale_length = 100.0;
    while(fabs(scale_length) > 0.05 && loopcounter < 100) // stops when the length is too small
    {
        //cout << "n = " << loopcounter << ", pB[0] = " << pB[0]
        //    << ", pB[1] = " << pB[1] << ", d[0] = " << distarray[0]
        //    << ", d[1] = " << distarray[1] << ", flip = " << flip
        //    << ", scale_length = " << scale_length << endl;
        if(distarray[0] > distarray[1])
        {
            if(loopcounter != 0)
            {
                if(flip == 1.0) scale = 0.4; // if minimization direction changes -> go back, but only the way * 0.4
                else scale = 0.7; // go on in this direction but only by the way * 0.7
            }
            scale_length = (pB[1]-pB[0])*scale; // the next length interval
            pB[0]     = pB[1] + scale_length; // the new path
            testB     = helixB.at(pB[0]); // new vector testB

            Float_t pathA_dca = -999.0;
            Float_t dcaAB_dca = -999.0;
            fHelixAtoPointdca_start_params(testB,helixA,pathA_dca,dcaAB_dca,path_in_A); // new helix to point dca calculation
            pA[0] = pathA_dca;
            //pA[0]     = helixA.pathLength(testB); // pathA at dca to helixB

            testA     = helixA.at(pA[0]); // new vector testA
            distarray[0] = (testA-testB).mag(); // new dca
            flip = 1.0;
        }
        else
        {
            if(loopcounter != 0)
            {
                if(flip == -1.0) scale = 0.4;
                else scale = 0.7;
            }
            scale_length = (pB[0]-pB[1])*scale;
            pB[1]     = pB[0] + scale_length;
            testB     = helixB.at(pB[1]);

            Float_t pathA_dca = -999.0;
            Float_t dcaAB_dca = -999.0;
            fHelixAtoPointdca_start_params(testB,helixA,pathA_dca,dcaAB_dca,path_in_A); // new helix to point dca calculation
            pA[1] = pathA_dca;
            //pA[1]     = helixA.pathLength(testB); // pathA at dca to helixB

            testA     = helixA.at(pA[1]); // pathA at dca to helixB
            distarray[1] = (testA-testB).mag();
            flip = -1.0;
        }
        loopcounter++;
    }
    if(distarray[0] < distarray[1])
    {
        pathB = pB[0];
        pathA = pA[0];
        dcaAB = distarray[0];
    }
    else
    {
        pathB = pB[1];
        pathA = pA[1];
        dcaAB = distarray[1];
    }
}

void StTriFlow2ndVertexFinder::fHelixABdca(StPhysicalHelixD helixA, StPhysicalHelixD helixB, Float_t &pathA, Float_t &pathB, Float_t &dcaAB)
{
    //cout << "Standard fHelixABdca called..." << endl;
    Float_t pA[2] = {0.0,0.0};
    Float_t pB[2] = {0.0,-70.0}; // the two start values for pathB, 0.0 is the origin of the helix at the first measured point
    Float_t distarray[2];
    StThreeVectorF testA, testB;
    for(Int_t r = 0; r < 2; r++)
    {
        testB     = helixB.at(pB[r]);  // 3D-vector of helixB point at path pB[r]

        Float_t pathA_dca = -999.0;
        Float_t dcaAB_dca = -999.0;
        fHelixAtoPointdca(testB,helixA,pathA_dca,dcaAB_dca); // new helix to point dca calculation
        testA = helixA.at(pathA_dca);
        //testA     = helixA.at(helixA.pathLength(testB)); // 3D-vector of helixA point at dca to testB

        distarray[r] = (testA-testB).mag(); // dca between helixA and helixB
    }
    Int_t loopcounter = 0;
    Float_t scale = 1.0;
    Float_t flip  = 1.0; // checks if the minimization direction changed
    Float_t scale_length = 100.0;
    while(fabs(scale_length) > 0.05 && loopcounter < 100) // stops when the length is too small
    {
        //cout << "n = " << loopcounter << ", pB[0] = " << pB[0]
        //    << ", pB[1] = " << pB[1] << ", d[0] = " << distarray[0]
        //    << ", d[1] = " << distarray[1] << ", flip = " << flip
        //    << ", scale_length = " << scale_length << endl;
        if(distarray[0] > distarray[1])
        {
            if(loopcounter != 0)
            {
                if(flip == 1.0) scale = 0.4; // if minimization direction changes -> go back, but only the way * 0.4
                else scale = 0.7; // go on in this direction but only by the way * 0.7
            }
            scale_length = (pB[1]-pB[0])*scale; // the next length interval
            pB[0]     = pB[1] + scale_length; // the new path
            testB     = helixB.at(pB[0]); // new vector testB

            Float_t pathA_dca = -999.0;
            Float_t dcaAB_dca = -999.0;
            fHelixAtoPointdca(testB,helixA,pathA_dca,dcaAB_dca); // new helix to point dca calculation
            pA[0] = pathA_dca;
            //pA[0]     = helixA.pathLength(testB); // pathA at dca to helixB

            testA     = helixA.at(pA[0]); // new vector testA
            distarray[0] = (testA-testB).mag(); // new dca
            flip = 1.0;
        }
        else
        {
            if(loopcounter != 0)
            {
                if(flip == -1.0) scale = 0.4;
                else scale = 0.7;
            }
            scale_length = (pB[0]-pB[1])*scale;
            pB[1]     = pB[0] + scale_length;
            testB     = helixB.at(pB[1]);

            Float_t pathA_dca = -999.0;
            Float_t dcaAB_dca = -999.0;
            fHelixAtoPointdca(testB,helixA,pathA_dca,dcaAB_dca); // new helix to point dca calculation
            pA[1] = pathA_dca;
            //pA[1]     = helixA.pathLength(testB); // pathA at dca to helixB

            testA     = helixA.at(pA[1]); // pathA at dca to helixB
            distarray[1] = (testA-testB).mag();
            flip = -1.0;
        }
        loopcounter++;
    }
    if(distarray[0] < distarray[1])
    {
        pathB = pB[0];
        pathA = pA[0];
        dcaAB = distarray[0];
    }
    else
    {
        pathB = pB[1];
        pathA = pA[1];
        dcaAB = distarray[1];
    }
}

void StTriFlow2ndVertexFinder::fHelixAtoPointdca_start_params(StThreeVectorF space_vec, StPhysicalHelixD helixA, Float_t &pathA, Float_t &dcaAB, Float_t path_in_A)
{
    Float_t pA[2] = {path_in_A+2,path_in_A-2}; // the two start values for pathB, 0.0 is the origin of the helix at the first measured point
    Float_t distarray[2];
    StThreeVectorF testA;
    for(Int_t r = 0; r < 2; r++)
    {
        testA     = helixA.at(pA[r]); // 3D-vector of helixA at path pA[r]
        distarray[r] = (testA-space_vec).mag(); // dca between helixA and helixB
    }
    Int_t loopcounter = 0;
    Float_t scale = 1.0;
    Float_t flip  = 1.0; // checks if the minimization direction changed
    Float_t scale_length = 100.0;
    while(fabs(scale_length) > 0.05 && loopcounter < 100) // stops when the length is too small
    {
        //cout << "n = " << loopcounter << ", pA[0] = " << pA[0]
        //    << ", pA[1] = " << pA[1] << ", d[0] = " << distarray[0]
        //    << ", d[1] = " << distarray[1] << ", flip = " << flip
        //    << ", scale_length = " << scale_length << endl;
        if(distarray[0] > distarray[1])
        {
            if(loopcounter != 0)
            {
                if(flip == 1.0) scale = 0.4; // if minimization direction changes -> go back, but only the way * 0.4
                else scale = 0.7; // go on in this direction but only by the way * 0.7
            }
            scale_length = (pA[1]-pA[0])*scale; // the next length interval
            pA[0]     = pA[1] + scale_length; // the new path
            testA     = helixA.at(pA[0]); // 3D-vector of helixA at path pA[0]
            distarray[0] = (testA-space_vec).mag(); // new dca
            flip = 1.0;
        }
        else
        {
            if(loopcounter != 0)
            {
                if(flip == -1.0) scale = 0.4;
                else scale = 0.7;
            }
            scale_length = (pA[0]-pA[1])*scale;
            pA[1]     = pA[0] + scale_length;
            testA     =  helixA.at(pA[1]); // 3D-vector of helixA at path pA[0]
            distarray[1] = (testA-space_vec).mag();
            flip = -1.0;
        }
        loopcounter++;
    }
    if(distarray[0] < distarray[1])
    {
        pathA = pA[0];
        dcaAB = distarray[0];
    }
    else
    {
        pathA = pA[1];
        dcaAB = distarray[1];
    }
}

void StTriFlow2ndVertexFinder::fHelixAtoPointdca(StThreeVectorF space_vec, StPhysicalHelixD helixA, Float_t &pathA, Float_t &dcaAB)
{
    // V1.1
    Float_t pA[2] = {0.0,-100.0}; // the two start values for pathB, 0.0 is the origin of the helix at the first measured point
    Float_t distarray[2];
    StThreeVectorF testA;
    for(Int_t r = 0; r < 2; r++)
    {
        testA     = helixA.at(pA[r]); // 3D-vector of helixA at path pA[r]
        distarray[r] = (testA-space_vec).mag(); // dca between helixA and helixB
    }
    Int_t loopcounter = 0;
    Float_t scale = 1.0;
    Float_t flip  = 1.0; // checks if the minimization direction changed
    Float_t scale_length = 100.0;
    while(fabs(scale_length) > 0.01 && loopcounter < 100) // stops when the length is too small
    {
        //cout << "n = " << loopcounter << ", pA[0] = " << pA[0]
        //    << ", pA[1] = " << pA[1] << ", d[0] = " << distarray[0]
        //    << ", d[1] = " << distarray[1] << ", flip = " << flip
        //    << ", scale_length = " << scale_length << endl;
        if(distarray[0] > distarray[1])
        {
            if(loopcounter != 0)
            {
                if(flip == 1.0) scale = 0.4; // if minimization direction changes -> go back, but only the way * 0.4
                else scale = 0.7; // go on in this direction but only by the way * 0.7
            }
            scale_length = (pA[1]-pA[0])*scale; // the next length interval
            pA[0]     = pA[1] + scale_length; // the new path
            testA     = helixA.at(pA[0]); // 3D-vector of helixA at path pA[0]
            distarray[0] = (testA-space_vec).mag(); // new dca
            flip = 1.0;
        }
        else
        {
            if(loopcounter != 0)
            {
                if(flip == -1.0) scale = 0.4;
                else scale = 0.7;
            }
            scale_length = (pA[0]-pA[1])*scale;
            pA[1]     = pA[0] + scale_length;
            testA     =  helixA.at(pA[1]); // 3D-vector of helixA at path pA[0]
            distarray[1] = (testA-space_vec).mag();
            flip = -1.0;
        }
        loopcounter++;
    }
    if(distarray[0] < distarray[1])
    {
        pathA = pA[0];
        dcaAB = distarray[0];
    }
    else
    {
        pathA = pA[1];
        dcaAB = distarray[1];
    }
    //cout << "pathA = " << pathA << ", dcaAB = " << dcaAB << endl;
}
//------------------------------------------------------------------------------------------------------------------

// for Lambda A is proton, B is pi_minus
// for antiLambda A is proton, B is pi_minus
void StTriFlow2ndVertexFinder::Find2ndVertex(StPhysicalHelixD helixA, StPhysicalHelixD helixB, StThreeVectorF vectorprim, Float_t MomentumA, Float_t MomentumB, TLorentzVector &lTrackA, TLorentzVector &lTrackB, Float_t &VerdistX, Float_t &VerdistY, StThreeVectorF &VecterAB, Float_t &DcaAB, Int_t mode)
{
  lTrackA.SetXYZM(0.0,0.0,0.0,-999.9);
  lTrackB.SetXYZM(0.0,0.0,0.0,-999.9);
  VerdistX = -999.9;
  VerdistY = -999.9;
  VecterAB.set(-999.9,-999.9,-999.9);
  DcaAB = -999.9;

  Float_t pathA_f, pathB_f, dcaAB_f;
  Float_t pathA_test = 0.0;
  Float_t pathB_test = 0.0;
  Int_t fDCA_Helix_out = fDCA_Helix_Estimate(helixA,helixB,pathA_test,pathB_test,dcaAB_f);
  
  StThreeVectorF vectoratsA = helixA.at(pathA_test);  // space vector of helixA at dca to helixB
  StThreeVectorF vectoratsB = helixB.at(pathB_test);  // space vector of helixB at dca to helixA
  StThreeVectorF vectorAB   = vectoratsA+vectoratsB;
  vectorAB   = vectorAB/2.0; // decay vertex


  StThreeVectorF baseA,dirA,baseB,dirB;
  baseA = helixA.at(pathA_test);
  baseB = helixB.at(pathB_test);
  dirA  = helixA.at(pathA_test-2.0) - helixA.at(pathA_test+2.0);
  dirB  = helixB.at(pathB_test-2.0) - helixB.at(pathB_test+2.0);

  StThreeVectorF vectorAB_lin  = calcVertexAnalytical(baseA,dirA,baseB,dirB); // vertex of the two tracks

  Double_t dcaAB_lin = calculateMinimumDistance(baseA,dirA,baseB,dirB);       // minimum distance between the two tracks
  StThreeVectorF vectorABtoPrim_lin = vectorAB_lin - vectorprim; // vector primary vertex to decay vertex
  Float_t VerdistX_lin = vectorABtoPrim_lin.mag(); // distance between primary vertex and decay vertex

  // calculate the scalar product with the approximated secondary vertex position
  StThreeVectorF vectornewA_lin     = helixA.cat(pathA_test); // direction vector at dca for helixA
  StThreeVectorF vectornewB_lin     = helixB.cat(pathB_test); // direction vector at dca for helixB
  vectornewA_lin = MomentumA*vectornewA_lin/vectornewA_lin.mag(); // new momentum vector at decay vertex
  vectornewB_lin = MomentumB*vectornewB_lin/vectornewB_lin.mag(); // new momentum vector at decay vertex

  TLorentzVector ltrackA_lin, ltrackB_lin;
  if(mode == 0 || mode == 1) // Lambda
  {
    ltrackA_lin.SetXYZM(vectornewA_lin.x(),vectornewA_lin.y(),vectornewA_lin.z(),TriFlow::mMassProton); // set Lorentz vector for proton
    ltrackB_lin.SetXYZM(vectornewB_lin.x(),vectornewB_lin.y(),vectornewB_lin.z(),TriFlow::mMassPion); // set Lorentz vector for pion
  }
  if(mode == 2) // K0S
  {
    ltrackA_lin.SetXYZM(vectornewA_lin.x(),vectornewA_lin.y(),vectornewA_lin.z(),TriFlow::mMassPion); // set Lorentz vector for pi+
    ltrackB_lin.SetXYZM(vectornewB_lin.x(),vectornewB_lin.y(),vectornewB_lin.z(),TriFlow::mMassPion); // set Lorentz vector for pi-
  }
  TLorentzVector trackAB_lin      = ltrackA_lin+ltrackB_lin; // mother particle

  StThreeVectorF dirY_lin(trackAB_lin.Px(),trackAB_lin.Py(),trackAB_lin.Pz());
  dirY_lin = dirY_lin/dirY_lin.mag();
  Double_t scalarProduct_lin = dirY_lin.dot(vectorABtoPrim_lin/vectorABtoPrim_lin.mag());

  if( VerdistX_lin > 1.8 && dcaAB_lin < 2.0 && scalarProduct_lin > 0.0 )
  {
    if(fDCA_Helix_out == 1)
    {
      fHelixABdca_start_params(helixA,helixB,pathA_f,pathB_f,dcaAB_f,pathA_test,pathB_test); // calculate dca between two helices
    }
    else
    {
      fHelixABdca(helixA,helixB,pathA_f,pathB_f,dcaAB_f); // calculate dca between two helices
    }

    vectoratsA     = helixA.at(pathA_f);  // space vector of helixA at dca to helixB
    vectoratsB     = helixB.at(pathB_f);  // space vector of helixB at dca to helixA
    vectorAB       = vectoratsA+vectoratsB;
    vectorAB       = vectorAB/2.0; // decay vertex

    StThreeVectorF vectorABtoPrim = vectorAB - vectorprim; // vector primary vertex to decay vertex
    Float_t verdistX = vectorABtoPrim.mag(); // distance between primary vertex and decay vertex

    StThreeVectorF vectornewA     = helixA.cat(pathA_f); // direction vector at dca for helixA
    StThreeVectorF vectornewB     = helixB.cat(pathB_f); // direction vector at dca for helixB

    vectornewA = MomentumA*vectornewA/vectornewA.mag(); // new momentum vector at decay vertex
    vectornewB = MomentumB*vectornewB/vectornewB.mag(); // new momentum vector at decay vertex

    TLorentzVector ltrackA, ltrackB;
    if(mode == 0 || mode == 1) // Lambda || anti-Lambda
    {
      ltrackA.SetXYZM(vectornewA.x(),vectornewA.y(),vectornewA.z(),TriFlow::mMassProton);
      ltrackB.SetXYZM(vectornewB.x(),vectornewB.y(),vectornewB.z(),TriFlow::mMassPion);
    }
    if(mode == 2) // K0S
    {
      ltrackA.SetXYZM(vectornewA.x(),vectornewA.y(),vectornewA.z(),TriFlow::mMassPion);
      ltrackB.SetXYZM(vectornewB.x(),vectornewB.y(),vectornewB.z(),TriFlow::mMassPion);
    }

    // Invariant mass calculations
    TLorentzVector trackAB = ltrackA+ltrackB; // mother particle

    StThreeVectorF dirY(trackAB.Px(),trackAB.Py(),trackAB.Pz());
    dirY = dirY/dirY.mag();
    Double_t scalarProduct = dirY.dot(vectorABtoPrim/vectorABtoPrim.mag());

    StThreeVectorF baseY = vectorAB;
    Double_t  verdistY  = calculateMinimumDistanceStraightToPoint(vectorAB,dirY,vectorprim);

    // Set TLorentzVector and Distence between primary to secondary vertex and dca from primary vertex to reconstract mother particle
    if(scalarProduct > 0)
    {
      lTrackA  = ltrackA;
      lTrackB  = ltrackB;
      VerdistX = verdistX;
      VerdistY = verdistY;
      VecterAB = vectorAB;
      DcaAB    = dcaAB_f;
    }
  }
}
