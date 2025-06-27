macro(Dimension=LinearAlgebra:-Dimension);
macro(Copy=LinearAlgebra:-Copy);
macro(HermiteForm=LinearAlgebra:-HermiteForm);
macro(RandomVector=LinearAlgebra:-RandomVector);
macro(RandomMatrix=LinearAlgebra:-RandomMatrix);
macro(MatrixInverse=LinearAlgebra:-MatrixInverse);

# Datastrcture for M for nontrivial S[i,i]

# myGcdex(a, b)
#
# Input: integers a, b
#
# Output: integers s,t,u,v,g such that
#
#     [ s  t ] [ a  0 ]   [ g  t ]
#     [      ] [      ] = [      ]  with g = igcd(a,b), v > 0 and 0 <= t < v
#     [ u  v ] [ b  1 ]   [ 0  v ]
#
myGcdex := proc(a,b)
  local g, u, v, s, t;

   if a=0 and b=0 then return 1,0,0,1,0 fi;
   if a=0 and b<>0 then return 0,sign(b),1,0,abs(b) fi;
   g := igcdex(b,a,'t');
   v := abs(iquo(a,g));
   t := modp(t,v);
   return iquo(g-t*b,a),t,-iquo(v*b,a),v,g

end:

# compactHcol(w, d)
#
# Input: w, an integer column vector of dimension n
#        d, a positive integer modulus 
#
# Output: n x i matrix H', i vector nonTrivialIdx, i
#   H' represent the matrices H below but only with nontrivial columns.
#   nonTrivialIdx represents the indices that is nontrivial.
#
#        [ d |         ]     [ g |         ]
#        [---+---------]     [---+---------]
#        [   |         ] --> [   |         ]
#        [ w |   I_n   ]     [   |    H    ]
#        [   |         ]     [   |         ]
# Where g := gcd(w, d)
#
compactHcol := proc(w,dd)
   local n,H,i,g,s,t,k,x,v,d, nonTrivialIdx, numNonTrivialIdx, tmp;

   n := Dimension(w);
   H := Matrix(n,n);
   hDiag := Vector(n);
   x := Copy(w);

   d := dd;
   g := d;
   numNonTrivialIdx := 0;
   nonTrivialIdx := [];
   for i from n by -1 to 1 do
      s[i],t[i],v[i],hDiag[i],g := myGcdex(g,x[i]); 
      if hDiag[i] <> 1 then
        ++numNonTrivialIdx;
        nonTrivialIdx := [op(nonTrivialIdx), i]
      end if;
   od;

  # Resort the nonTrivialIdx and discard the trivial indices.
  nonTrivialIdx := Vector(numNonTrivialIdx, (i) -> nonTrivialIdx[numNonTrivialIdx - i + 1]);
  H := Matrix(n, numNonTrivialIdx, (i, j) -> if (i = nonTrivialIdx[j]) then hDiag[nonTrivialIdx[j]] else 0 end if);

   # Construct T[1,1], T[1,2],T[2,2], T[1,3],T[2,3],T[3,3], ... and apply
   for k to numNonTrivialIdx do
    d := iquo(d, hDiag[k]);
    idx := nonTrivialIdx[k];
    for i to idx - 1 do
      H[i,k] := modp(-t[idx]*modp(x[i],H[idx,k]),H[idx,k]); 
      x[i] := modp(iquo(x[i]+H[i,k]*x[idx],H[idx,k]),d);
    end do;
    for i from idx+1 to n do x[i] := iquo(x[i], H[idx,k]) end do;
  end do;
   
  return H, nonTrivialIdx, numNonTrivialIdx;
end:


# hcol(w, d)
#
#  Input: w, an integer column vector of dimension n
#         d, a positive integer modulus 
#
# Output: n x n matrix H such that
#
#        [ d |         ]     [ g |         ]
#        [---+---------]     [---+---------]
#        [   |         ] --> [   |         ],
#        [ w |   I_n   ]     [   |    H    ]
#        [   |         ]     [   |         ]
#         
#       with H in (row) Hermite form, and the matrix on
#       the right of "-->" being left equivalent over Z
#       to the matrix on the left.
#
hcol := proc(w,dd)
   local n,H,i,g,s,t,k,x,v,d;

   n := Dimension(w);
   H := Matrix(n,n);
   x := Copy(w);

   d := dd;
   g := d;

   for i from n by -1 to 1 do
      s[i],t[i],v[i],H[i,i],g := myGcdex(g,x[i]); 
   od;

   # ARNE: Correction made Sat Nov 30, 2024.
   #       - code is optimized to use the assumption gcd(d,w)=1
   #       - to correct, remove any content from [d,w]
   if g>1 then 
       x := map(iquo,x,g);   
       d := iquo(d,g);
   fi;

   # Construct T[1,1], T[1,2],T[2,2], T[1,3],T[2,3],T[3,3], ... and apply
   for k to n do
      if H[k,k]=1 then next fi; # col(T,k) is identity vector
      d := iquo(d,H[k,k]);
      for i to k-1 do
         H[i,k] := modp(-t[k]*modp(x[i],H[k,k]),H[k,k]); 
         x[i] := modp(iquo(x[i]+H[i,k]*x[k],H[k,k]),d);
      od;
      for i from k+1 to n do x[i] := iquo(x[i],H[k,k]) od;
   od;

   return H;
end:
