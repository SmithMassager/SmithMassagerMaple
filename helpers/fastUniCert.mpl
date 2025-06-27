read "macros.mpl";

# UniCert(A, n)
#
# Input:
#   A, n x n integer matrix in Z, nonsingular.
#
# Output:
#   True if |det(A)| = 1 else false
#
fastUniCert := proc(A, n)
  local Abar, _;
  if (not type(A, 'Matrix'(integer))) then
    return false;
  end if;

  #Abar, _, _ := RowPL(A, n, n);

  #return uniCertHelper(Abar);
  return uniCertHelper(A);
end proc; 
