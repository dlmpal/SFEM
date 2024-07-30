//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {-1, 0, 0, 1.0};
//+
Point(3) = {0, 1, 0, 1.0};
//+
Circle(1) = {2, 1, 3};
//+
Point(4) = {-2, 0, 0, 1.0};
//+
Point(5) = {-2, 2, 0, 1.0};
//+
Point(6) = {0, 2, 0, 1.0};
//+
Line(2) = {4, 2};
//+
Line(3) = {5, 4};
//+
Line(4) = {6, 5};
//+
Line(5) = {3, 6};
//+
Curve Loop(1) = {3, 4, 5, 1, 2};
//+
Plane Surface(1) = {1};
//+
Physical Curve("Left", 6) = {3};
//+
Physical Curve("Top", 7) = {4};
//+
Physical Curve("Right", 8) = {5};
//+
Physical Curve("Cylinder", 9) = {1};
//+
Physical Curve("Bottom", 10) = {2};
//+
Physical Surface("Fluid", 11) = {1};
