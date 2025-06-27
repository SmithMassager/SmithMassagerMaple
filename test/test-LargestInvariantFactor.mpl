read "helpers/largestInvariantFactor.mpl";


#------------------------ Begin LargestInvariantFactor ------------------------
checkLargestInvariantFactor := proc()
  local A, s;
  randomize(1);
  A := Matrix([[1,0], [0,3]]);
  s := LargestInvariantFactor(A, 2);
  ASSERT(s = 3);
end proc;
#-------------------------- End LargestInvariantFactor ------------------------

checkLargestInvariantFactor();
