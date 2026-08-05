BeginPackage["Quantum`"]

RandomGaussianComplex::usage = "RandomGaussianComplex[] returns a complex number \
with independent standard-normal real and imaginary parts. \
RandomGaussianComplex[mu, sigma] shifts and scales by mu and sigma."

RandomGaussian::usage = "RandomGaussian[] returns a standard-normal real random number."

RandomState::usage = "RandomState[dim] returns a Haar-random pure state vector of \
dimension dim (complex vector drawn from the unit sphere)."

TwoRandomOrhtogonalStates::usage = "TwoRandomOrhtogonalStates[dim] returns two \
Haar-random orthonormal pure states via Gram-Schmidt."

Proyector::usage = "Proyector[psi] returns the projector |psi><psi|. \
Proyector[psi, phi] returns the operator |psi><phi|."

GUEMember::usage = "GUEMember[n] returns an n x n GUE random Hermitian matrix. \
Option Normalization -> \"Default\" (factor 1), \"SizeIndependent\" (0.5), \
or \"Kohler\" ((n/2)/Pi^2)."

GOEMember::usage = "GOEMember[n] returns an n x n GOE random real symmetric matrix. \
Option Normalization -> \"Default\" (factor 1), \"f90\" (1/Sqrt[2]), \
\"SizeIndependent\" (0.5), or \"Kohler\" ((n/2)/Pi^2)."

CUEMember::usage = "CUEMember[n] returns a Haar-random n x n unitary matrix (CUE)."

PUEMember::usage = "PUEMember[n] returns an n x n random Hermitian matrix with \
eigenvalues uniform on [-1/2, 1/2] and CUE-random eigenvectors."

RandomHermitianMatrix::usage = "RandomHermitianMatrix[qubits] returns a GUE random \
Hermitian matrix of dimension 2^qubits. Option Normalization -> \"Default\", \
\"SizeIndependent\", or \"Kohler\"."

RandomDensityMatrix::usage = "RandomDensityMatrix[n] returns a random n x n density \
matrix drawn from the Hilbert-Schmidt measure."

RandomMixedState::usage = "RandomMixedState[n, k] returns a random mixed state of \
dimension n as a convex combination of k Haar-random pure states."

BasisState::usage = "BasisState[k, dim] returns the k-th computational basis vector of \
dimension dim (0-indexed: k runs from 0 to dim-1)."

Bell::usage = "Bell[n] returns the maximally entangled state (1/Sqrt[n]) Sum_k |k,k> \
as a vector of length n^2. Bell[] is Bell[2] = (|00> + |11>)/Sqrt[2]."

GHZ::usage = "GHZ[n] returns the n-qubit GHZ state (|00...0> + |11...1>)/Sqrt[2]."

Wstate::usage = "Wstate[n] returns the n-qubit W state \
(|10...0> + |01...0> + ... + |00...1>)/Sqrt[n]."

CoherentState::usage = "CoherentState[theta, phi, n] returns the n-qubit spin coherent \
state, i.e. the n-fold tensor product of the single-qubit state \
{Cos[theta/2], Exp[I phi] Sin[theta/2]}."

BlochSphereToState::usage = "BlochSphereToState[{x, y, z}] converts a Bloch sphere point \
in Cartesian coordinates to a qubit pure state. BlochSphereToState[{theta, phi}] uses \
spherical angles directly."

OrthonormalBasisContaningVector::usage = "OrthonormalBasisContaningVector[psi] returns an \
orthonormal basis for the Hilbert space such that the first vector is psi (up to a global \
phase chosen to make the first component real and positive)."

Begin["`Private`"]

(* {{{ Random complex and real scalars *)

(* {{{ *) RandomGaussianComplex[] :=
  RandomReal[NormalDistribution[]] + I RandomReal[NormalDistribution[]]

RandomGaussianComplex[mu_, sigma_] :=
  mu + sigma RandomGaussianComplex[]
(* }}} *)

(* {{{ *) RandomGaussian[] := RandomReal[NormalDistribution[]]
(* }}} *)

(* }}} *)

(* {{{ Random states and matrices *)

(* {{{ *) RandomState[dim_Integer] :=
  Normalize[Table[RandomGaussianComplex[], {dim}]]
(* }}} *)

(* {{{ *) TwoRandomOrhtogonalStates[dim_Integer] :=
  Module[{psi1, psi2},
    psi1 = RandomState[dim];
    psi2 = RandomState[dim];
    psi2 = psi2 - (Conjugate[psi1] . psi2) psi1;
    {psi1, Normalize[psi2]}]
(* }}} *)

(* {{{ *) Proyector[psi_?VectorQ, phi_?VectorQ] := Outer[Times, psi, Conjugate[phi]]
Proyector[psi_?VectorQ] := Outer[Times, psi, Conjugate[psi]]
(* }}} *)

(* {{{ *) Options[GUEMember] = {Normalization -> "Default"};
GUEMember[n_Integer, opts : OptionsPattern[]] :=
  Module[{m, scale},
    m     = Table[RandomGaussianComplex[], {n}, {n}];
    scale = Switch[OptionValue[Normalization],
      "Default",         1,
      "SizeIndependent", 1/2,
      "Kohler",          (n/2) / Pi^2];
    scale (m + ConjugateTranspose[m])]
(* }}} *)

(* {{{ *) Options[GOEMember] = {Normalization -> "Default"};
GOEMember[n_Integer, opts : OptionsPattern[]] :=
  Module[{m, scale},
    m     = Table[RandomGaussian[], {n}, {n}];
    scale = Switch[OptionValue[Normalization],
      "Default",         1,
      "f90",             1/Sqrt[2],
      "SizeIndependent", 1/2,
      "Kohler",          (n/2) / Pi^2];
    scale (m + Transpose[m])]
(* }}} *)

(* {{{ *) CUEMember[n_Integer] := RandomVariate[CircularUnitaryMatrixDistribution[n]]
(* }}} *)

(* {{{ *) PUEMember[n_Integer] :=
  Module[{u},
    u = CUEMember[n];
    Chop[u . DiagonalMatrix[RandomReal[{-1/2, 1/2}, n]] . ConjugateTranspose[u]]]
(* }}} *)

(* {{{ *) Options[RandomHermitianMatrix] = {Normalization -> "Default"};
RandomHermitianMatrix[qubits_Integer, opts : OptionsPattern[]] :=
  Module[{n, m, scale},
    n     = 2^qubits;
    m     = Table[RandomGaussianComplex[], {n}, {n}];
    scale = Switch[OptionValue[Normalization],
      "Default",         1,
      "SizeIndependent", 1/2,
      "Kohler",          2^(qubits - 1) / Pi^2];
    scale (m + ConjugateTranspose[m])]
(* }}} *)

(* {{{ *) RandomDensityMatrix[n_Integer] :=
  Module[{m},
    m = GUEMember[n];
    # / Tr[#] & [m . ConjugateTranspose[m]]]
(* }}} *)

(* {{{ *) RandomMixedState[n_Integer, combinations_Integer] :=
  Module[{states, weights},
    states  = Table[Proyector[RandomState[n]], {combinations}];
    weights = # / Total[#] & [RandomReal[{0, 1}, combinations]];
    Chop[weights . states]]
(* }}} *)

(* }}} *)

(* {{{ Named states *)

(* {{{ *) BasisState[k_Integer, dim_Integer] := UnitVector[dim, k + 1]
(* }}} *)

(* {{{ *) Bell[n_Integer] := Flatten[IdentityMatrix[n]] / Sqrt[n]
Bell[] := Bell[2]
(* }}} *)

(* {{{ *) GHZ[n_Integer] :=
  Module[{dim = 2^n},
    (UnitVector[dim, 1] + UnitVector[dim, dim]) / Sqrt[2]]
(* }}} *)

(* {{{ *) Wstate[n_Integer] := Sum[BasisState[2^m, 2^n], {m, 0, n-1}] / Sqrt[n]
(* }}} *)

(* {{{ *) CoherentState[theta_, phi_, n_Integer] :=
  Flatten[tensorPower[N[{Cos[theta/2], Exp[I phi] Sin[theta/2]}], n]]
(* }}} *)

(* }}} *)

(* {{{ Bloch sphere *)

(* {{{ *) BlochSphereToState[{x_, y_, z_}] :=
  Module[{sph = Quiet[CoordinateTransform["Cartesian" -> "Spherical", {x, y, z}], ArcTan::indet]},
    BlochSphereToState[{sph[[2]], sph[[3]] /. Indeterminate -> 0}]]

BlochSphereToState[{theta_, phi_}] := {Cos[theta/2], Exp[I phi] Sin[theta/2]}
(* }}} *)

(* }}} *)

(* {{{ Basis construction *)

(* {{{ *) OrthonormalBasisContaningVector[psi_?VectorQ] :=
  Module[{n = Length[psi], ri, Ri, F, eks},
    F = Sum[
      ri = Normalize[Table[RandomGaussianComplex[], {n}]];
      Ri = ri - (Conjugate[psi] . ri) psi;
      Proyector[Ri], {n - 1}] + Proyector[psi];
    eks = Transpose[Sort[Transpose[Eigensystem[F]],
      Abs[#1[[1]] - 1] < Abs[#2[[1]] - 1] &]][[2]];
    Exp[I (Arg[psi[[1]]] - Arg[eks[[1, 1]]])] Normalize[#] & /@ eks]
(* }}} *)

(* }}} *)

End[]

EndPackage[]
