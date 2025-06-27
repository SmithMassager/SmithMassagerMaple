#interface(quiet=true);
read("SmithMassager/SmithMassager.mpl");
read("HermiteForm/genherm.mpl");
kernelopts(numcpus=60);

randomize(3);

n := nextprime(100);
resetprofile();
resetGlobalTime();
profile(SmithMassager, IntCertificate, UniCert, IndexMassager, RowPL, SmithFormModS, getbAdic);
A := interestHerm(n);
U, M, T, S := SmithMassager(A, n);
showprofile();
printGlobalTime();
#lprint(A);
#lprint(U);
#lprint(M);
#lprint(T);
#lprint(S);
