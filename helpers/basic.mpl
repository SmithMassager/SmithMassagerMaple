# --> myGcdex(a,b)
#
# Input: integers a, b
#
# Output: integers s,t,u,v,g such that
#
#     [ s  t ] [ a    ]   [ g  t ]
#     [      ] [      ] = [      ]  with g = igcd(a,b), v > 0 and 0 <= t < v
#     [ u  v ] [ b  1 ]   [    v ]
#
myGcdex := proc(a,b)
  local g, u, v, s, t;

   if a=0 then 
      if b=0 then 
         return 1,0,0,1,0
      else 
         return 0,sign(b),1,0,abs(b)
      fi;
   fi;
   g := igcdex(b,a,'t');
   v := abs(iquo(a,g));
   t := modp(t,v);
   return iquo(g-t*b,a),t,-iquo(v*b,a),v,g

end:

# --> Ann(a,N)
#
# Input: integers a and N, N > 0
#
# Output: a principal generator of the ideal { b | b a = 0 } in Z/(N)
#
Ann := proc(a,N)
    local g,s,t,u,v;

    s,t,u,v,g := myGcdex(a,N);
    if s=0 then
        return 1
    else
        return modp(igcd(u,N),N)
    fi;

end:

# NOTE: This alg is wrong, the fix is presented on section 4 of 2005 'The modular extended gcd ...'
#
# Conditioner(a, b, N, d, q)
#
# Input:
#   a, b, N, q, positive integers.
#   d, list of factorizations of N of size q, not necessarily unique.
#
# Pre:
#   gcd(a, b) = 1
#
# Output:
#   Either (1) finds a minimal c s.t. (a + cb, N) = (a, b, N).
#	   (2) finds a larger factorization of N compared to d.
#
Conditioner := proc(a, b, N, d, q)
  local K, i, g, T, A, B, t, j;
  K := floor(2 * log2(N)^(3/2));
  T := Array(0..K, 1);
  A := Vector(q);
  B := Vector(q);

  for i to q do
    A[i] := modp(a, d[i]);
    B[i] := modp(b, d[i]);
    g := gcd(B[i], d[i]);
    if (2 <= g and g < d[i]) then 
      ret := [op(d), g];
      ret[i] := iquo(d[i], g);
      return ret;
    end if;
  end do;

  for i to q do
    if (B[i] <> 0) then
      s := modp(-A[i]/B[i], d[i]);
      if (s > K) then next end if;
      for j from 0 to iquo((K-s), d[i]) do
	T[s + j * d[i]] := 0;
      end do;
    end if;
  end do;

  for t from 0 to K do
    if (T[t] = 1) then break end if;
  end do;

  for i to q do
    g := gcd(A[i] + t*B[i], d[i]);
    if (2 <= g and g < d[i]) then 
      ret := [op(d), g];
      ret[i] := iquo(d[i], g);
      return ret;
    end if;
  end do;

  return t;
end proc;

# --> Stab(a, b, N)
#
# Input: integers a, b, N, N > 0
#
# Output: nonegative integer c such that
#   gcd(a, b, N) = gcd(a + c*b, N)
#   c is minimal, i.e. c <= ceil(2 * log(N2)^3/2)
#
Stab := proc(a, b, N)
 local g, g2, N2, a2, c, e:

 # Recall: (a/(g*g2) + cb/(g*g2), N/g) = (1) iff (a + cb, N) = (g).

 # Take the gcd of all the arguments
 g := igcd(a, b, N):
 if g = 0 then return 0 fi:

 # Take the gcd of a, b
 g2 := igcd(iquo(a,g), iquo(b,g)):
 if g2 = 0 then return 0 fi:

 # Divide all with g and update a, so that gcd(a, b) = 1
 N2 := iquo(N, g);
 a2 := iquo(a, g * g2);
 b2 := iquo(b, g * g2);

 s := a2;
 if (igcd(s, N2) = 1) then return 0 end if;
 for c from 1 to N2-1 do
  s += b2;
  if (igcd(s, N2) = 1) then return c end if;
 end do;

 #if (c > 2 * ceil(log2(N2)^(3/2))) then print("c in stab is too big"); ASSERT(false); end if;
 #return c;

 #d := [N2];
 #q := 1;
 #while(true) do
 # d := Conditioner(a2, b2, N2, d, q);
 # if (type(d, 'integer')) then return d end if;
 # q := nops(d);
 #end do;

 ## Compute c, by taking the largest part of N2 which has no common factors with a2
 #e := 2;
 #while (e < N) do
 # e := e * e;
 # a2 := modp(a2 * a2, N2);
 #od;

 #c := modp(iquo(N2, igcd(a2, N2)), N2);

 #return c:
end:

# --> Modinv(a, N)
#
# Input: integers a, N, N > 0
#
# Output: integer e such that modp(e * a, N) = gcd(m, N)
#
Modinv := proc(a, N)
 local g, u, v:

 g := igcdex(a, N, 'u');
 v := iquo(N,g);
 return u + Stab(u, v, N) * v;

end:

# --> ExtendedStab(n, a, m, j)
#
# Input:
#   a, n vector mod m.
#   m, positive integer modulus.
#
# Output: c
#  c, n vector mod m s.t. (c1a1 + ... + cnan, m) = (a1, ..., an, m).
#     c_j = 1
#
ExtendedStab := proc(n, a, m, j := 1)
  local c, prev, i;
  c := Vector(n);
  c[j] := 1;
  prev := a[j];

  for i from 1 to n do
    if (i = j) then
      next;
    end if;
    c[i] := Stab(prev, a[i], m);
    prev := modp(prev + c[i] * a[i], m);
    if (prev = 1) then break end if;
  end do;

  return c;
end proc;

#
# MinLargestPower(a, X)
#
# Input:
#   a, non-negative integer.
#   X, positive integer.
#
# Output: q
#   q is the smallest integer s.t. X^(q) > a 
#
MinLargestPower := proc(a, X)
  return nops(convert(a, base, X));
end proc;

#
# getDT(a)
#
# Input:
#   a, non-negative integer.
#
# Output: dt
#  dt := float[8] if a < 2^25-1
#        anything  else
#
getDT := proc(a);
  if (a < 33554432) then
    return float[8];
  else
    return anything;
  end if;
end proc;
