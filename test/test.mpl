read "helpers/doublePlusOneLift.mpl";

A := Matrix([[3, 0], [0, 9]]);
X := 1000;
#A := Matrix([[47, 31], [29, 74]]);
#X := 2;
A0, R, M := DoublePlusOneLift(A,X, 2, 3);
i := Matrix([[1,0],[0,1]]);
print(1/X * (i - A . modp(A0, X)));

R0 := convert(R[0], Matrix);
M0 := convert(M[0], Matrix);
R1 := convert(R[1], Matrix);
M1 := convert(M[1], Matrix);
R2 := convert(R[2], Matrix);
M2 := convert(M[2], Matrix);
X0 := doublePlusOne(X, 0);
X1 := doublePlusOne(X, 1);
X2 := doublePlusOne(X, 2);
X3 := doublePlusOne(X, 3);


B1 := modp(A0 . (i + R0 *X0) + M0 * X0^2, X1);
print(B1);
print(modp(A^(-1), X1));

B2 := modp(B1 . (i + R1 * X1) + M1 * X1^2, X2);
print(B2);
print(modp(A^(-1), X2));

#print(modp(B1 . (i + R1 * X1) , X1^2));
#print(modp(A^(-1), X1^2));

#B3 := modp(B2 . (i + R2 * X2) + M2 * X2^2, X3);
#print(B3);
#print(modp(A^(-1), X2));
