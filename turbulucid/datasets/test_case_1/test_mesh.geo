// Gmsh project created on Thu Dec 22 14:16:48 2016
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {1, 1, 0, 1.0};
//+
Point(4) = {0, 1, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Point(5) = {1, -1, 0, 1.0};
//+
Point(6) = {2, -1, 0, 1.0};
//+
Circle(5) = {2, 5, 6};
//+
Point(7) = {2, 1, 0, 1.0};
//+
Line(6) = {6, 7};
//+
Line(7) = {7, 3};
//+
Point(8) = {2.5, -1, 0, 1.0};
//+
Point(9) = {4.5, -1, 0, 1.0};
//+
Point(10) = {4.5, 1, 0, 1.0};
//+
Point(11) = {6.5, 1, 0, 1.0};
//+
Line(8) = {7, 10};
//+
Line(9) = {10, 11};
//+
Line(10) = {11, 9};
//+
Line(11) = {9, 8};
//+
Line(12) = {8, 6};
//+
Line(13) = {8, 10};
//+
Line Loop(14) = {3, 4, 1, 2};
//+
Plane Surface(15) = {14};
//+
Line Loop(16) = {6, 7, -2, 5};
//+
Plane Surface(17) = {16};
//+
Line Loop(18) = {6, 8, -13, 12};
//+
Plane Surface(19) = {18};
//+
Line Loop(20) = {9, 10, 11, 13};
//+
Plane Surface(21) = {20};
//+
Transfinite Line {2, 4, 6, 13, 10} = 10 Using Progression 1;
//+
Transfinite Line {3, 1, 7, 5, 8, 12, 11, 9} = 10 Using Progression 1;
//+
Transfinite Surface {15};
//+
Transfinite Surface {17};
//+
Transfinite Surface {21};

Recombine Surface {15, 17, 21};

Extrude {0, 0, 1} {Surface{15, 17, 19, 21}; Layers{10}; Recombine; }
//+
Delete {
  Point{5};
}
//+
Physical Surface("right") = {43, 65, 87, 109};
//+
Physical Surface("left") = {15, 17, 19, 21};
//+
Physical Surface("inlet") = {34};
//+
Physical Surface("outlet") = {100};

Physical Volume("internal") = {1, 2, 3, 4};
//+
Physical Surface("top") = {30, 56, 78, 96};
//+
Physical Surface("botOrthoHex") = {38};
//+
Physical Surface("botCurved") = {64};
//+
Physical Surface("botPrism") = {86};
//+
Physical Surface("botSkewedHex") = {104};
