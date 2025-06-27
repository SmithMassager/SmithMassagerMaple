interface(quiet=true);
read("SmithMassager/SmithMassager.mpl");
read("HermiteForm/genherm.mpl");
kernelopts(numcpus=60);

randomize(3);

#n := nextprime(500);
#A := interestHerm(n);
#
#profile(SmithMassager, IntCertificate, UniCert, IndexMassager);
#U, M, T, S := SmithMassager(A, n);
#showprofile();
#writeto("smithMassager1-new.out");
#lprint(A);
#lprint(U);
#lprint(M);
#lprint(T);
#lprint(S);
#writeto(terminal):

#A := 2 * A;
#
#resetprofile();
#profile(SmithMassager, IntCertificate, UniCert, IndexMassager);
#U, M, T, S := SmithMassager(A, n);
#showprofile();
#writeto("smithMassager2-new.out");
#lprint(A);
#lprint(U);
#lprint(M);
#lprint(T);
#lprint(S);
#writeto(terminal):
#
#n := nextprime(250);
#resetprofile();
#profile(SmithMassager, IntCertificate, UniCert, IndexMassager);
#A := interestHerm(n);
#U, M, T, S := SmithMassager(A, n);
#showprofile();
#writeto("smithMassager3-new.out");
#lprint(A);
#lprint(U);
#lprint(M);
#lprint(T);
#lprint(S);
#writeto(terminal):
#
#A := 2 * A;
#
#resetprofile();
#profile(SmithMassager, IntCertificate, UniCert, IndexMassager);
#U, M, T, S := SmithMassager(A, n);
#showprofile();
#writeto("smithMassager4-new.out");
#lprint(A);
#lprint(U);
#lprint(M);
#lprint(T);
#lprint(S);
#writeto(terminal):


#resetprofile();
#profile(SmithMassager, IntCertificate, UniCert, IndexMassager);
#n := nextprime(750);
#A := interestHerm(n);
#U, M, T, S := SmithMassager(A, n);
#showprofile();
#writeto("smithMassager5-new.out");
#lprint(A);
#lprint(U);
#lprint(M);
#lprint(T);
#lprint(S);
#writeto(terminal):

#A := interestHerm(n);
#A := 2 * A;
#
#resetprofile();
#profile(SmithMassager, IntCertificate, UniCert, IndexMassager);
#U, M, T, S := SmithMassager(A, n);
#showprofile();
#writeto("smithMassager6-new.out");
#lprint(A);
#lprint(U);
#lprint(M);
#lprint(T);
#lprint(S);
#writeto(terminal):


#n := nextprime(1200);
#resetprofile();
#profile(SmithMassager, IntCertificate, UniCert, IndexMassager, RowPL, SmithFormModS, getbAdic);
#A := interestHerm(n);
#U, M, T, S := SmithMassager(A, n);
#showprofile();
#writeto("smithMassager11-new.out");
#lprint(A);
#lprint(U);
#lprint(M);
#lprint(T);
#lprint(S);
#writeto(terminal):
#
#A := 2 * A;
#
#resetprofile();
#profile(SmithMassager, IntCertificate, UniCert, IndexMassager, RowPL, SmithFormModS, getbAdic);
#U, M, T, S := SmithMassager(A, n);
#showprofile();
#writeto("smithMassager12-new.out");
#lprint(A);
#lprint(U);
#lprint(M);
#lprint(T);
#lprint(S);
#writeto(terminal):
#
#n := nextprime(1500);
#resetprofile();
#profile(SmithMassager, IntCertificate, UniCert, IndexMassager, RowPL, SmithFormModS, getbAdic);
#A := interestHerm(n);
#U, M, T, S := SmithMassager(A, n);
#showprofile();
#writeto("smithMassager13-new.out");
#lprint(A);
#lprint(U);
#lprint(M);
#lprint(T);
#lprint(S);
#writeto(terminal):
#
#A := 2 * A;
#
#resetprofile();
#profile(SmithMassager, IntCertificate, UniCert, IndexMassager, RowPL, SmithFormModS, getbAdic);
#U, M, T, S := SmithMassager(A, n);
#showprofile();
#writeto("smithMassager14-new.out");
#lprint(A);
#lprint(U);
#lprint(M);
#lprint(T);
#lprint(S);
#writeto(terminal):

#n := nextprime(2200);
#resetprofile();
#profile(SmithMassager, IntCertificate, UniCert, IndexMassager, RowPL, SmithFormModS, getbAdic);
#A := interestHerm(n);
#U, M, T, S := SmithMassager(A, n);
#showprofile();
#writeto("smithMassager15-new.out");
#lprint(A);
#lprint(U);
#lprint(M);
#lprint(T);
#lprint(S);
#writeto(terminal):
#
#A := 2 * A;
#
#resetprofile();
#profile(SmithMassager, IntCertificate, UniCert, IndexMassager, RowPL, SmithFormModS, getbAdic);
#U, M, T, S := SmithMassager(A, n);
#showprofile();
#writeto("smithMassager16-new.out");
#lprint(A);
#lprint(U);
#lprint(M);
#lprint(T);
#lprint(S);
#writeto(terminal):


#n := nextprime(4700);
n := nextprime(1700);
#n := nextprime(400);
resetprofile();
resetGlobalTime();
profile(SmithMassager, IntCertificate, UniCert, IndexMassager, RowPL, SmithFormModS, getbAdic, imlSolveHelper, highOrderResidue, extractMatrixGCD, extractRowGcd, reScaleEntry, eliminateCol, computeProjBasis, LinearAlgebra:-Modular:-Inverse);
A := interestHerm(n);
A := A;
zz = time[real]();
U, M, T, S := SmithMassager(A, n);
print(max(abs(U)));
zzz = time[real]();
a = zzz-zz;
print("SmithMassager time:", a);
showprofile();
printGlobalTime();
writeto("SmithMassager/smithMassager18-new.out");
lprint(A);
lprint(U);
lprint(M);
lprint(T);
lprint(S);
writeto(terminal):

#n := nextprime(2500);
#resetprofile();
#resetGlobalTime();
#profile(SmithMassager, IntCertificate, UniCert, IndexMassager, RowPL, SmithFormModS, getbAdic, imlSolveHelper, highOrderResidue, extractMatrixGCD, extractRowGcd, reScaleEntry, eliminateCol);
#A := interestHerm(n);
#U, M, T, S := SmithMassager(A, n);
#showprofile();
#printGlobalTime();
#writeto("SmithMassager/smithMassager18-new.out");
#lprint(A);
#lprint(U);
#lprint(M);
#lprint(T);
#lprint(S);
#writeto(terminal):
