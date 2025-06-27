#  Input: n - a positive integer >= 10
#
# Output: an n x n integer matrix with a row Hermite form that hopefully
#         has many nontrivial diagonal entries
#
randHerm := proc(n)
     local A,f,g,h,i,r1,r2,s;

     A := Matrix(n,n):
     f := rand(1..floor(n/10)):
     g := rand(1..n):
     h := rand(0..1):
     for i to n do A[i,i] := f() od:
     A[1..floor(n/2),1..floor(n/2)] := LinearAlgebra[RandomMatrix](floor(n/2),floor(n/2)):
     for i to 10*n do
          r1 := g(); r2 := g(); while r2 = r1 do r2 := g() od; s := 1-2*h();
          A[r1,1..-1] := A[r1,1..-1]+s*A[r2,1..-1];
          r1 := g(); r2 := g(); while r2 = r1 do r2 := g() od; s := 1-2*h();
          A[1..-1,r1] := A[1..-1,r1]+s*A[1..-1,r2];
     od:

    return A;

end:

#  Input: n - a prime
#
# Output: A - an n x n nonisngular integer matrix with interesting Smith form
#
interestHerm := proc(n)
    local A,i,j;

    A := Matrix(n,n):
    for i to n do for j to n do A[i,j] := modp((i-1)^(j-1),n) od od:

    return A;

end:
