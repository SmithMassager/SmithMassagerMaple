read "SmithMassager/SmithMassager.mpl";
read "helpers/diagHelpers.mpl";

# TrivialLowerHermiteForm(B, n)
#
# Input:
#   n, positive integer.
#   B, n x n nonsingular integer matrix.
#
# Output: vector(h, f) given the hermite form is trivial else return false.
#   Lower triangular hermite form of B with shape:
#   [h        ]
#   [f I_(n-1)]
#
TrivialLowerHermiteForm := proc(B, n)
  local ret, M, S, _, s, h;

  ret := SmithMassager(B, n);

  if (ret = false) then return false; end if;

  _, M, _, S := ret;
  s := S[n];
  if (S[n-1] <> 1 or igcd(s,M[1,n]) <> 1) then return false end if;

  # Calculate Rem(-M[2..n,n]/M[1,n], s).
  h := Vector(n);
  igcdex(M[1,n], s, 'inv');
  ModMultiply(s, modp(-inv, s), M[1..n, n], h);
  h[1] := s;

  return h;
end proc;

# SmithFormMultipliers(A, n)
#
# Input:
#   A, n x n nonsingular integer matrix.
#   n, positive integer.
#
# Output: Return (S, U, V) If the multipliers are found else return false.
#   S, n x n integer matrix in smith form.
#   U, V, n x n integer unimodular matrix.
#   AV = US
#
SmithFormMultipliers := proc(A, n)
  local ret, M, TwoS, S, detTwoS, lambda, R, B, H, V, hInv, _, U;

  ret := SmithMassager(2 * A, n);

  if (ret = false) then return false; end if;

  _, M, _, TwoS := ret;
  S := TwoS / 2;

  detTwoS := DiagDet(TwoS, n);

  lambda := 105 * max(n, ceil(root(detTwoS, n)));
  R := RandomMatrix(n, n, generator = 0..lambda-1);
  # Calculate B := M + R TwoS
  R := DiagMult(R, TwoS, n, n);
  B := Add(M, R);
  
  H := TrivialLowerHermiteForm(B, n);
  if (H = false) then return false; end if;

  # Calculate V := BH^-1 by exploiting the structure of H.
  V := Matrix(n, n);
  hInv := Vector(n, (i) -> if i = 1 then 1/H[1] else -H[i]/H[1] end if);
  V[1..n, 1] := Multiply(B, hInv);
  V[1..n, 2..n] := B[1..n, 2..n];

  U := DiagMult(Multiply(A, V), DiagInv(S, n), n, n);

  return S, U, V;
end proc;
