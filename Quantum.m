(* ::Package:: *)

(* {{{ *) BeginPackage["Quantum`"]
(* {{{ TODO
Crear un swap matrix en el mismo espiritu que la CNOT. 
}}} *)
(* {{{ Bitwise manipulation and basic math*)
ExtractDigits::usage = "Extract the digits of a number in its binary form. Same in spirity as exdig in math.f90
   with numin is NumberIn and nwhich is Location digits,
   This routine takes numin and puts two numbers
   n1out and n2out which result from the digits
   of numin that are marked with the 1 bits
   of the number nwhich
   nwhich=   0 1 0 1 0 0 1 = 42
   numin=    0 1 1 0 1 1 1 = 55
   n1out=      1   0     1 = 5
   n2out=    0   1   1 1   = 7, {n1out,n2out}"
BitsOn::usare = "The number of 1s in a binary representation, say 
   n = 5 = 1 0 1
   BitsOn[n] = 2"
MergeTwoIntegers::usage = "MergeTwoIntegers[na_, nb_, ndigits_]
  Merge two intergers based on the positions given by the 1s of the 
  third parameter. Functions similar to the one written in fortran, merge_two_integers.
  In this routine we merge two numbers. It is useful for doing the tensor
  product. 'a' and 'b' are the input number wheares as usual ndigits
  indicates the position. Example
  ndigits            = 0 0 1 0 1 0 0 1 = 41
  a                  =     1   0     1 = 5
  b                  = 1 0   0   1 0   = 18
  merge_two_integers = 1 0 1 0 0 1 0 1 = 165"
xlog2x::usage = "Calculates x Log[2, x], but if x is 0, takes the correct limit"
(* }}} *)
(* {{{ Random States and matrices*)
TwoRandomOrhtogonalStates::usage = "TwoRandomOrhtogonalStates[dim_] creates two random states orthogonal to each other using gram schmidt process. It is useful for creating states that provide a maximally mixed state in a qubit"
RandomState::usage = "A random state RandomState[dim_Integer]"
RandomGaussianComplex::usage = "RandomGaussianComplex[] gives a complex random number with Gaussian distribution
   centered at 0 and with width 1"
RandomGaussian::usage = "RandomGaussian[] gives a random number with Gaussian distribution
   centered at 0 and with width 1"
RandomHermitianMatrix::usage = "RandomHermitianMatrix[ ] To generate a GUE Random Hermitian Matrix, various normalizatios are possible"
PUEMember::usage = "PUEMember[n_] To generate a PUE Random Hermitian Matrix, with spectral span from -1/2 to 1/2. n is the dimension"
GUEMember::usage = "GUEMember[n_] To generate a GUE Random Hermitian Matrix, various normalizatios are possible. n is the dimension"
GOEMember::usage = "GOEMember[n_] To generate a GOE Random Hermitian Matrix, various normalizatios are possible. n is the dimension"
CUEMember::usage = "To generate a CUE Random Unitary Matrix: CUEMember[n_]. n is the dimension"
RandomDensityMatrix::usage = "RandomDensityMatrix[n_] to generate a random density matrix of length n"
Normalization::usage = "default sucks, SizeIndependent yields  Hij Hkl = delta il delta jk, kohler gives
                          = N/Pi^2 delta il delta jk" 
(* }}} *)
(* {{{ Matrix manipulation and advanced linear algebra*)
BaseForHermitianMatrices::usage = "Base element for Hermitian matrices"
Commutator::usage = "Makes the commutator between A and B: Commutator[A,B]=A.B-B.A"
PartialTrace::usage = "Yields the partial trace of a given state. The second argument indicates the digits to leave. PartialTrace[Rho_?MatrixQ, DigitsToLeave_] or PartialTrace[Psi_?VectorQ, LocationDigits_]"
PartialTranspose::usage = "Takes the partial transpososition with  respect
  to the indices specified, PartialTranspose[A_, TransposeIndices_]"
PartialTransposeFirst::usage = "Transpose the first part of a bipartite system with equal dimensions"
PartialTransposeSecond::usage = "Transpose the Second part of a bipartite system with equal dimensions"
DirectSum::usage = "DirectSum[A1, A2] A1 \[CirclePlus] A2 or DirectSum[{A1, A2, ..., An}] = A1 \[CirclePlus] A2\[CirclePlus] ... \[CirclePlus] An"
tensorProduct::usage = "Creates the tensor product for qubits, like TensorProduct_real_routine in linear.f90
      I think that the inidices in the number refer to where to place the second matrix."
tensorPower::usage = "Creates the tensor product of n matrices"
OrthonormalBasisContaningVector::usage=" OrthonormalBasisContaningVector[psi_?VectorQ] will create an orthonorlam basis that contains the given vector"
GetMatrixForm[Gate_, Qubits_] := Gate /@ IdentityMatrix[Power[2, Qubits]]
(* }}} *)
(* {{{ Quantum information routines *)
ApplyControlGate::usage = "ApplyControlGate[U_, State_, TargetQubit_, ControlQubit_]"
ApplyGate::usage = "ApplyGate[U_, State_, TargetQubit_]"
ApplyChannel::usage = "Apply Quantum Channel to an operator"
ControlNotMatrix::usage = "Get the control not matrix ControlNotMatrix[QubitControl_Integer, QubitTarget_Integer, NumberOfQubits_Integer]"
QuantumDotProduct::usage = "Dot Product in which the first thing is congutate. Yields the true Psi1 dot Psi2" 
ToSuperOperatorSpaceComputationalBasis::usage = "ToSuperOperatorSpaceComputationalBasis[rho_] converts a density matrix rho to superoperator space, where it is a vector, and the basis in the superoperator space is the computational basis, which is the same (or at least similar) to the Choi basis"
FromSuperOperatorSpaceComputationalBasis::usage = "FromSuperOperatorSpaceComputationalBasis[rho_?VectorQ] converts a density matrix rho from superoperator space, where it is a vector, and the basis in the superoperator space is the computational basis, which is the same (or at least similar) to the Choi basis to the normal space, where it is a matrix"
FromSuperOperatorSpacePauliBasis::usage = "FromSuperOperatorSpacePauliBasis[r_?VectorQ] maps a vetor to a density, as the inverse of ToSuperOperatorSpacePauliBasis"
ToSuperOperatorSpacePauliBasis::usage = "ToSuperOperatorSpacePauliBasis[r_?MatrixQ] maps a density
matrix to superoperator space using the Pauli matrices and tensor products"
BlochNormalFormVectors::usage = "Vectors to calculate normal form, see arXiv:1103.3189v1"
BlochNormalFormFromVectors::usage = "Density matrix in normal form, see arXiv:1103.3189v1"
StateToBlochBasisRepresentation::usage = "The so called Bloch matrix, see arXiv:1103.3189v1"
CoherentState::usage = "CoherentState[\[Theta]_, \[Phi]_, n_] Makes a Coherent Spin state where the inputs are the Bloch angles and the number of qubits"
Pauli::usage = "Pauli[0-3] gives Pauli Matrices according to wikipedia, and Pauli[{i1,i2,...,in}] gives Pauli[i1] \[CircleTimes]Pauli[i2] \[CircleTimes] ... \[CircleTimes] Pauli[in]"
Paulib::usage = "Pauli[b_] gives b.Sigma, if b is a 3 entry vector"
SumSigmaX::usage = "Matrix correspoding to sum_j Pauli[1]_j"
SumSigmaY::usage = "Matrix correspoding to sum_j Pauli[2]_j"
SumSigmaZ::usage = "Matrix correspoding to sum_j Pauli[3]_j"
ValidDensityMatrix::usage = "Test whether a given density matrix is valid"
Dagger::usage = "Hermitian Conjugate"
Concurrence::usage = "Concurrence of a 2 qubit density matrix density matrix"
MultipartiteConcurrence::usage = "Concurrence of multipartite entanglement for pure states"
Purity::usage = "Purity of a 2 qubit density matrix density matrix"
Proyector::usage = "Gives the density matrix rho=|Psi><Psi| corresponging to |Psi>, or rho=|Phi><Psi| se se le dan dos argumentos: Proyector[Phi, Psi]"
KrausOperatorsForSuperUnitaryEvolution::usage = "Gives the Kraus Operators for a given unitary and a given state of the environement"
ApplyKrausOperators::usage = "Apply a set of Kraus Operors, (for example fro the output of KrausOperatorsForSuperUnitaryEvolution to a state. Synthaxis is ApplyKrausOperators[Operators_, rho_]. "
vonNeumannEntropy::usage = "von Neumann entropy of a density Matrix, or a list of vectors. Usage vonNeumannEntropy[r_]"
Bell::usage = "Bell[n] state gives a maximally entangled state of the type 
sum_i^n |ii\>. Bell[ ] gives |00> + |11>" 
StateToBlochSphere::usage="Get the Bloch sphere point of a mixed state"
BlochSphereToState::usage="BlochSphereToState[CartesianCoordinatesOfPoint_List] transforms the points in the Bloch Sphere to a mixed state. Alternatively, one can give the angles, so BlochSphereToState[{\[Theta]_, \[Phi]_}]"
QuantumMutualInformationReducedEntropies::usage="QuantumMutualInformationReducedEntropies[r_?MatrixQ] calcula la informacion mutia con entropias"
QuantumMutualInformationMeasurements::usage="QuantumMutualInformationMeasurements[r_?MatrixQ] calcula la informacion mutua cuando 
permitimos unla entropia condicional de una matriz de dos qubits cuando hacemos una medicion sobre un qubit. "
Discord::usage="Discord[r_?MatrixQ] calcula el quantum discord de una matriz de dos qubits."
BasisState::usage="BasisState[BasisNumber_,Dimension_] gives a state of the computational basis. It is 
numbered from 0 to Dimension-1"
Hadamard::usage="Hadamard[] gate, Hadamard[QubitToApply, Total"
(* }}} *)
(* {{{ Quantum channels and basis transformations *)
JamiolkowskiStateToOperatorChoi::usage = "JamiolkowskiStateToOperatorChoi[Rho] applies the Jamiolkovski isomorphism as is understood in Geometry of Quantum States pgs. 241, 243, 266"
JamiolkowskiOperatorChoiToState::usage = "JamiolkowskiOperatorChoiToState[O] applies the inverse Jamiolkovski isomorphism as is understood in Geometry of Quantum States"
TransformationMatrixPauliBasisToComputationalBasis::usage = "The matrix that allows to tranform, in superoperator space, from the Pauli basis (the GellMann basis for dimension 2 modulo order) to the computational basis, aka the Choi basis"
Reshuffle::usage = "Apply the reshufle operation as undestood in Geometry of Quantum States, pags 260-262, and 264"
RandomTracePreservingMapChoiBasis::usage = "Creates a singlequbit random trace preserving map"
AveragePurityChannelPauliBasis::usage = "Calculates the average final purity given that the initial states are pure, and chosen with the Haar measure. "
BlochEllipsoid::usage = "BlochEllipsoid[Cha_] Show the deformation of BlochSphere for a qubit channel in the pauli basis"
EvolvGate::usage="EvolvGate[Gate_, steps_, env_, state_]... Evoluciona cualquier estado aplicando un numero steps_ de veces la compuerta Gate de  la forma Gate[#, otherparameters_] donde debe ponerse # en el lugar donde Gate toma el estado"
MakeQuantumChannel::usage="MakeQuantumChannel[Gate_, steps_, env_] Donde Gate va de la forma Gate[#, otherparameters_]  donde debe ponerse # en el lugar donde Gate toma el estado"
SumPositiveDerivatives::usage="SumPositiveDerivatives[list_] Suma todas las contribuciones list(max)-list(min) sucesivos cuando la derivada es positiva"
GHZ::usage="GHZ[qu_] Creates a GHZ state, |000000\[RightAngleBracket]+|111111\[RightAngleBracket]"
(* }}} *)
(* }}} *)
Begin["Private`"] 
(* {{{ Bitwise manipulation and basic math*)
(* {{{ *) ExtractDigits[NumberIn_, LocationDigits_] := Module[{AuxArray, NumberOfDigits}, 
	NumberOfDigits = IntegerLength[NumberIn, 2]; 
	AuxArray = Transpose[{IntegerDigits[LocationDigits, 2, NumberOfDigits], 
	IntegerDigits[NumberIn, 2, NumberOfDigits]}];
	{FromDigits[Select[AuxArray, #[[1]] != 1 &][[All, 2]], 2], 
	FromDigits[Select[AuxArray, #[[1]] == 1 &][[All, 2]], 2]}] 

(* }}} *)
(* {{{ *) BitsOn[n_] := Count[IntegerDigits[n, 2], 1]
(* }}} *)
(* {{{ *) MergeTwoIntegers[na_, nb_, ndigits_] := 
 Module[{LongitudTotal, Digits01s, Result}, 
  LongitudTotal = 
   Max[Count[IntegerDigits[ndigits, 2], 0], BitLength[nb]] + 
    BitsOn[ndigits];
  Digits01s = 
   PadRight[Reverse[IntegerDigits[ndigits, 2]], LongitudTotal];
  Result = PadRight[{}, LongitudTotal];
  Result[[Flatten[Position[Digits01s, 1]]]] = 
   Reverse[IntegerDigits[na, 2, Count[Digits01s, 1]]];
  Result[[Flatten[Position[Digits01s, 0]]]] = 
   Reverse[IntegerDigits[nb, 2, Count[Digits01s, 0]]];
  FromDigits[Reverse[Result], 2]]


(* }}} *)
(* {{{ *) xlog2x[x_] := If[x == 0, 0, x Log[2, x]]

(* }}} *)
(* }}} *)
(* {{{ Random States and matrices*)
(* {{{ *) TwoRandomOrhtogonalStates[dim_] := Module[{psi1, psi2, prepsi2},
  psi1 = RandomState[dim];
  prepsi2 = RandomState[dim];
  prepsi2 = prepsi2 - Dot[Conjugate[psi1], prepsi2] psi1;
  psi2 = #/Norm[#] &[prepsi2];
  {psi1, psi2}]
(* }}} *)
(* {{{ *) RandomState[dim_Integer] := #/(Sqrt[Conjugate[#].#])&[Table[ RandomGaussianComplex[], {i, dim}]]; 
(* }}} *)
(* {{{ *) RandomGaussianComplex[] := #[[1]] + #[[2]] I &[RandomReal[NormalDistribution[], {2}]];
(* }}} *)
(* {{{ *) RandomGaussianComplex[A_,sigma_] := sigma(#[[1]] + #[[2]] I &[RandomReal[NormalDistribution[], {2}]])+A;
(* }}} *)
(* {{{ *) RandomGaussian[] :=RandomReal[NormalDistribution[0, 1]];
(* }}} *)
(* {{{ *) RandomHermitianMatrix[qubits_Integer,   OptionsPattern[]] :=
	Switch[OptionValue[Normalization],"Default",1,"SizeIndependent",0.5, "Kohler", 
	Power[2,qubits-1]/Power[Pi,2]] ((Transpose[Conjugate[#]] + #)&) @@ 
	{Table[ RandomGaussianComplex[], {i, Power[2, qubits]}, {j, Power[2, qubits]}]};
Options[RandomHermitianMatrix] = {Normalization -> "Default"};

(* }}} *)
(* {{{ *) PUEMember[n_Integer] := Module[{U},
	U = CUEMember[n];
	Chop[U.DiagonalMatrix[Table[Random[], {n}]-0.5].Dagger[U]]]
(* }}} *)
(* {{{ *) GUEMember[n_Integer,   OptionsPattern[]] :=
	Times[Switch[OptionValue[Normalization],"Default",1,"SizeIndependent",0.5, "Kohler", (n/2)/Power[Pi,2]], 
	((Transpose[Conjugate[#]] + #)&) @@ {Table[ RandomGaussianComplex[], {i, n}, {j, n}]}];
Options[GUEMember] = {Normalization -> "Default"};

(* }}} *)
(* {{{ *) GOEMember[n_Integer,   OptionsPattern[]] :=
(* With the normalization "f90" we get consistent values as with the module *)
(* In particular, we obtain that <W_{ij} W_{kl}>=\delta_{il}\delta_{jk}+\delta{ik}\delta_{jl} *)
	Times[Switch[OptionValue[Normalization],"Default",1,"f90",1/Sqrt[2],"SizeIndependent",0.5, "Kohler", (n/2)/Power[Pi,2]], 
	((Transpose[#] + #)&) @@ {Table[ RandomGaussian[], {i, n}, {j, n}]}];
Options[GOEMember] = {Normalization -> "Default"};

(* }}} *)
(* {{{ *) CUEMember[n_] := Transpose[ Inner[Times, Table[Exp[I Random[Real, 2 \[Pi]]], {n}], 
	Eigenvectors[GUEMember[n]], List]]


(* }}} *)
(* {{{ *) RandomDensityMatrix[n_] := #/Tr[#] &[(#.Transpose[Conjugate[#]]) &[GUEMember[n]]]

(* }}} *)
(* }}} *)
(* {{{ Matrix manipulation and advanced linear algebra*)
(* {{{ *) BaseForHermitianMatrices[j_Integer, Ncen_Integer] := 
 If[j <= Ncen, SparseArray[{j, j} -> 1, {Ncen, Ncen}], 
  SparseArray[{{Ncen - Ceiling[TriangularRoot[j - Ncen]], 
      1 + Ncen - (j - Ncen - Triangular[Ceiling[TriangularRoot[j - Ncen]] - 1])}, 
      {1 + Ncen - (j - Ncen - Triangular[Ceiling[TriangularRoot[j - Ncen]] - 1]), 
      Ncen - Ceiling[TriangularRoot[j - Ncen]]}} -> 1/Sqrt[2], {Ncen, Ncen}]]
(* {{{ TriangularRoot and Triangular*) 
TriangularRoot[n_] := (-1 + Sqrt[8 n + 1])/2
Triangular[n_] := n (n + 1)/2
(* }}} *)

(* }}} *)
(* {{{ *) Commutator[A_?MatrixQ, B_?MatrixQ] := A.B - B.A
(* }}} *)
(* {{{ *) PartialTrace[Rho_?MatrixQ, DigitsToLeave_] := Module[{ab1, ab2, na, nb},
	nb = Power[2, BitsOn[DigitsToLeave]];
	na = Length[Rho]/nb;
	Table[
		Sum[
		(*Print[a,b1,n,DigitsToLeave,{MergeTwoIntegers[b1-1,a-1,
		 DigitsToLeave],MergeTwoIntegers[a-1,b1-1,n]},{MergeTwoIntegers[
		 b2-1,a-1,DigitsToLeave],MergeTwoIntegers[a-1,b2-1,n]}];*)

		ab1 = MergeTwoIntegers[b1 - 1, a - 1, DigitsToLeave];
	ab2 = MergeTwoIntegers[b2 - 1, a - 1, DigitsToLeave];
	Rho[[ab1 + 1, ab2 + 1]], {a, na}],
		{b1, nb}, {b2, nb}]]
(* }}} *)
(* {{{ *) PartialTrace[Psi_?VectorQ, LocationDigits_] := Module[{DimHCentral, MatrixState, i}, 
	DimHCentral = Power[2, DigitCount[LocationDigits, 2, 1]];
	MatrixState = SparseArray[Table[1 + ExtractDigits[i - 1, LocationDigits] -> Psi[[i]], {i,Length[Psi]}]];
	Table[QuantumDotProduct[MatrixState[[All, j2]], MatrixState[[All, j1]]], 
		{j1, DimHCentral}, {j2, DimHCentral}]]


(* }}} *)
(* {{{ *) PartialTranspose[A_?MatrixQ, TransposeIndices_Integer] := 
 Module[{i, j, k, l, il, kj},
  Table[
   {i, j} = {BitAnd[ij - 1, BitNot[TransposeIndices]], 
     BitAnd[ij - 1, TransposeIndices]};
   {k, l} = {BitAnd[kl - 1, BitNot[TransposeIndices]], 
     BitAnd[kl - 1, TransposeIndices]};
   {il, kj} = {i + l, k + j};
   A[[il + 1, kj + 1]], {ij, Length[A]}, {kl, Length[A]}]]
(* }}} *)
(* {{{ *) PartialTransposeFirst[A_?MatrixQ] := 
 Module[{n, i, j, k, l, ij, kl, kj, il},
  n = Sqrt[Length[A]];
  Table[
   {i, j} = IntegerDigits[ij, n, 2];
   {k, l} = IntegerDigits[kl, n, 2];
   kj = FromDigits[{k, j}, n];
   il = FromDigits[{i, l}, n];
   A[[kj + 1, il + 1]], {ij, 0, n^2 - 1}, {kl, 0, n^2 - 1}]
  ]
(* }}} *)
(* {{{ *) PartialTransposeSecond[A_?MatrixQ] := 
 Module[{n, i, j, k, l, ij, kl, kj, il},
  n = Sqrt[Length[A]];
  Table[
   {i, j} = IntegerDigits[ij, n, 2];
   {k, l} = IntegerDigits[kl, n, 2];
   kj = FromDigits[{k, j}, n];
   il = FromDigits[{i, l}, n];
   A[[il + 1, kj + 1]], {ij, 0, n^2 - 1}, {kl, 0, n^2 - 1}]
  ]
(* }}} *)
(* {{{ *) DirectSum[MatrixList_List] := Fold[DirectSum, MatrixList[[1]], Drop[MatrixList, 1]]
(* }}} *)
(* {{{ *) DirectSum[A1_, A2_] := Module[{dims},
  dims = Dimensions /@ {A1, A2};
  ArrayFlatten[{{A1, Table[0, {dims[[1, 1]]}, {dims[[2, 2]]}]},
    {Table[0, {dims[[2, 1]]}, {dims[[1, 2]]}], A2}}]]
(* }}} *)
(* {{{ *) tensorProduct[Matrix1_?MatrixQ, Matrix2_?MatrixQ] := KroneckerProduct[Matrix1, Matrix2]
(* }}} *)
(* {{{ *) tensorProduct[LocationDigits_Integer, Matrix1_?MatrixQ, Matrix2_?MatrixQ] := 
	Module[{Indices, iRow, iCol, L1, L2}, 
		{L1,L2}=Length/@{Matrix1,Matrix2};
		Normal[SparseArray[Flatten[Table[
			Indices = {1 + ExtractDigits[iRow - 1, LocationDigits], 
			1 + ExtractDigits[iCol - 1, LocationDigits]};
			{iRow, iCol} -> Part @@ Join[{Matrix1}, Indices[[All, 1]]] Part @@ Join[{Matrix2}, 
			Indices[[All, 2]]], {iRow, L1 L2}, {iCol, L1 L2}]]]]];
(* }}} *)
(* {{{ *) tensorProduct[LocationDigits_, State1_?VectorQ, State2_?VectorQ] := 
	Module[{Index, i, L1, L2}, 
		{L1, L2} = Length /@ {State1, State2};
		Normal[SparseArray[Table[Index = 1 + ExtractDigits[i - 1, LocationDigits]; 
		i -> State1[[Index[[1]]]] State2[[Index[[2]]]], {i, L1 L2}]]]];
(* }}} *)
(* {{{ *) tensorProduct[State1_?VectorQ, State2_?VectorQ] := Flatten[KroneckerProduct[State1, State2]]
(* }}} *)
(* {{{ *) tensorPower[A_, n_] := Nest[KroneckerProduct[A, #] &, A, n - 1]
(* }}} *)
(* {{{ *) OrthonormalBasisContaningVector[psi_?VectorQ] := 
(*El Ri generado es orthogonal a psi y tiene norma menor que uno*)
 Module[{n, ri, i, Ri, F, eks},
  n = Length[psi];
  F = Sum[
     ri = #/(Sqrt[Conjugate[#].#]) &[
       Table[RandomGaussianComplex[], {i, n}]];
     Ri = ri - (Conjugate[psi].ri) psi;
     {Conjugate[Ri].psi, Conjugate[Ri].Ri};
     Proyector[Ri], {n - 1}] + Proyector[psi];
  eks = Transpose[
     Sort[Transpose[Eigensystem[F]], 
      Abs[#1[[1]] - 1] < Abs[#2[[1]] - 1] &]][[2]];
  Exp[I (Arg[psi[[1]]] - Arg[eks[[1, 1]]])] #/Norm[#] & /@
    eks]
(* }}} *)
(* {{{  *) GetMatrixForm[Gate_, Qubits_] := Transpose[Gate /@ IdentityMatrix[Power[2, Qubits]]]
(* }}} *)
(* }}} *)
(* {{{ Quantum information routines *)
(* {{{ *) ApplyChannel[Es_, rho_] := 
 Sum[Es[[k]].rho.Dagger[Es[[k]]], {k, Length[Es]}]
(* }}} *)
(* {{{ *) Hadamard[] := 
 {{1,1},{1,-1}}/Sqrt[2]
(* }}} *)
(* {{{ *) ControlNotMatrix[QubitControl_Integer, QubitTarget_Integer, NumberOfQubits_Integer] :=
 SparseArray[
  Table[{NumeroACambiar, ControlNot[QubitControl, QubitTarget, NumeroACambiar - 1] + 1}, 
  	{NumeroACambiar, Power[2, NumberOfQubits]}] -> 1]
ControlNot[QubitControl_Integer, QubitTarget_Integer, NumeroACambiar_Integer] := 
 Module[{NumeroEnBinario, DigitoControl, DigitoTarget, n},
  NumeroEnBinario = IntegerDigits[NumeroACambiar, 2];
  n = Max[Length[NumeroEnBinario], QubitControl + 1, 
    QubitTarget + 1];
  NumeroEnBinario = IntegerDigits[NumeroACambiar, 2, n];
  DigitoControl = NumeroEnBinario[[-QubitControl - 1]];
  DigitoTarget = NumeroEnBinario[[-QubitTarget - 1]];
  NumeroEnBinario[[-QubitTarget - 1]] = 
   Mod[DigitoControl + DigitoTarget, 2];
  FromDigits[NumeroEnBinario, 2]]
(* }}} *)
(* {{{ *) QuantumDotProduct[Psi1_?VectorQ, Psi2_?VectorQ] :=  Dot[Conjugate[Psi1], Psi2]
(* }}} *)
(* {{{ *) ToSuperOperatorSpaceComputationalBasis[rho_?MatrixQ] := 
            Flatten[rho]
(* }}} *)
(* {{{ *) FromSuperOperatorSpaceComputationalBasis[rho_?VectorQ] := 
 Partition[rho, Sqrt[Length[rho]]]
(* }}} *)
(* {{{ *) ToSuperOperatorSpacePauliBasis[r_?MatrixQ] := Module[{Qubits},
  Qubits = Log[2, Length[r]];
  Chop[Table[ Tr[r.Pauli[IntegerDigits[i, 4, Qubits]]/Power[2.,Qubits/2]], {i, 0, Power[Length[r], 2] - 1}]]]
(* }}} *)
(* {{{ *) FromSuperOperatorSpacePauliBasis[r_?VectorQ] := Module[{Qubits},
  Qubits = Log[2, Length[r]]/2; 
  Sum[r[[i + 1]] Pauli[IntegerDigits[i, 4, Qubits]]/
	  Power[2., Qubits/2], {i, 0, Length[r] - 1}]
   ]
(* }}} *)
(* {{{  *) BlochNormalFormVectors[r_?MatrixQ] := Module[{R, T, oa, cd, ob, A, a, b, c},
  R = Chop[StateToBlochBasisRepresentation[r]];
  T = R[[2 ;;, 2 ;;]];
  {oa, cd, ob} = SingularValueDecomposition[T];
  A = Chop[
    DirectSum[{{1}}, Transpose[oa/Det[oa]]].R.DirectSum[{{1}}, 
      ob/Det[ob]]];
  {a = A[[2 ;;, 1]], b = A[[1, 2 ;;]], c = Diagonal[A][[2 ;;]]}
  ]
(* }}} *)
(* {{{  *) BlochNormalFormFromVectors[{a_, b_, c_}] := Module[{i},
  (IdentityMatrix[4] + Sum[b[[i]] Pauli[{0, i}], {i, 1, 3}] + 
   Sum[a[[i]] Pauli[{i, 0}], {i, 1, 3}] + 
   Sum[c[[i]] Pauli[{i, i}], {i, 1, 3}])/4
  ]


(* }}} *)
(* {{{  *) StateToBlochBasisRepresentation[r_?MatrixQ] := Table[Tr[r.Pauli[{i, j}]], {i, 0, 3}, {j, 0, 3}]
(* }}} *)
(* {{{ *) CoherentState[\[Theta]_, \[Phi]_, n_] := 
 N[Flatten[
   tensorPower[{Cos[\[Theta]/2], Exp[I \[Phi]] Sin[\[Theta]/2]}, n], 
   1]]
(* }}} *)
(* {{{ Pauli matrices*)  
Pauli[0]=Pauli[{0}]={{1,0}, {0,1}}; 
Pauli[1]=Pauli[{1}]={{0,1}, {1,0}}; 
Pauli[2]=Pauli[{2}]={{0,-I},{I,0}}; 
Pauli[3]=Pauli[{3}]={{1,0}, {0,-1}};
Pauli[Indices_List] := KroneckerProduct @@ (Pauli /@ Indices)
(* }}} *)
(* {{{ *) Paulib[{b1_,b2_,b3_}] := b1 Pauli[1]+ b2 Pauli[2]+ b3 Pauli[3]
(* }}} *)
(* {{{ *) SumSigmaX[Qubits_] := SumSigmaX[Qubits] = Table[If[DigitCount[BitXor[i - 1, j - 1], 2, 1] == 1, 1, 0], {i, Power[2, Qubits]}, {j, Power[2, Qubits]}]
(* }}} *)
(* {{{ *) SumSigmaY[Qubits_] := SumSigmaY[Qubits] = Table[ If[DigitCount[BitXor[i - 1, j - 1], 2, 1] == 1, If[i > j, I, -I], 0], {i, Power[2, Qubits]}, {j, Power[2, Qubits]}]
(* }}} *)
(* {{{ *) SumSigmaZ[Qubits_] := SumSigmaZ[Qubits] = DiagonalMatrix[ Table[Qubits - 2 DigitCount[i, 2, 1], {i, 0, Power[2, Qubits] - 1}]]
(* }}} *)
(* {{{ *) ValidDensityMatrix[Rho_?MatrixQ] := (Abs[Tr[Rho] - 1]<Power[10,-13] && 
	Length[Select[# >= 0 & /@ Chop[Eigenvalues[Rho]], ! # &]] == 0)
(* }}} *)
(* {{{ *) Dagger[H_]:=Conjugate[Transpose[H]]

(* }}} *)
(* {{{ *) Concurrence[rho_?MatrixQ]:=Module[{lambda}, 
	lambda=Sqrt[Abs[Eigenvalues[rho.sigmaysigmay.Conjugate[rho].sigmaysigmay]]]; 
	Max[{2*Max[lambda]-Plus@@lambda,0}]]
Concurrence[\[Psi]_?VectorQ]:=Module[{}, 
Abs[QuantumDotProduct[sigmaysigmay.Conjugate[\[Psi]], \[Psi]]]]

(* {{{ *) sigmaysigmay={{0, 0, 0, -1}, {0, 0, 1, 0}, {0, 1, 0, 0}, {-1, 0, 0, 0}}; (* }}} *)

(* {{{ *) MultipartiteConcurrence[\[Psi]_?VectorQ] := Module[{M},M=Length[\[Psi]]-2;
Sqrt[M-Sum[Purity[PartialTrace[\[Psi],i]],{i,1,M}]]*(2/Sqrt[M+2])](* }}} *)
(* }}} *)
(* {{{ *) Purity[rho_]:= Tr[rho.rho]

(* }}} *)
(* {{{ Proyector *)
Proyector[psi_, phi_] := Outer[Times, psi, Conjugate[phi]]
Proyector[psi_] := Proyector[psi, psi]
(* }}} *)
(* {{{ *) ApplyKrausOperators[Operators_, rho_] := 
 Total[#.rho.Dagger[#] & /@ Operators]
(* }}} *)
(* {{{ *) KrausOperatorsForSuperUnitaryEvolution[psienv_, U_] := 
  Module[{Nenv, Ntotal, Nsub},
   (*
   See Chuang and Nielsen p. 360
   So, We have that E(rho)=sum_k Ek rho Ek(dagger);
   One can see that Ek=<ek|U|psienv>, 
   that is a matrix with matrix elements
   Ek_{ij}=<i|Ek|j> = <i ek|U|psienv j>;
   Metiendo una identidad obrenemos que 
   Ek_{ij}=sum_m <ek i|U|m j > <m|psienv>;
   *)
   Nenv = Length[psienv];
   Ntotal = Length[U];
   Nsub = Ntotal/Nenv;
   Table[(*Print["k="<>ToString[k]];*)
    
    Table[(*Print["i="<>ToString[i]<>" j="<>ToString[j]];*)
     Sum[
      (*If[k==2,Print["k="<>ToString[k]<>" i="<>ToString[i]<>" j="<>
      ToString[j]<>" m="<>ToString[m]<>" m="<>ToString[{Nsub k+i+1,
      Nsub m+j+1,m+1}]],];*)
      
      U[[Nsub k + i + 1, Nsub m + j + 1]] psienv[[m + 1]], {m, 0, 
       Nenv - 1}], {i, 0, Nsub - 1}, {j, 0, Nsub - 1}],
    {k, 0, Nenv - 1}]];
(* }}} *)
(* {{{ *) vonNeumannEntropy[r_?MatrixQ] := vonNeumannEntropy[Eigenvalues[r]]
(* }}} *)
(* {{{ *) vonNeumannEntropy[r_?ListQ] := -Total[If[# == 0, 0, # Log[2, #]] & /@ r]
(* }}} *)
(* {{{ *) Bell[n_] := (1/Sqrt[n]) Table[If[Mod[i - 1, n + 1] == 0, 1, 0], {i, n^2}]
(* }}} *)
(* {{{ *) Bell[] = Bell[2]

(* }}} *)
(* {{{  *) StateToBlochSphere[R_?MatrixQ] := Module[{sigma}, sigma={Pauli[1], Pauli[2], Pauli[3]}; Chop[Tr /@ (sigma.R)]]
(* }}} *)
(* {{{  *) BlochSphereToState[CartesianCoordinatesOfPoint_List] :=  
 Module[{r, \[Theta], \[Phi]}, 
  r = CoordinateTransform["Cartesian"->"Spherical", CartesianCoordinatesOfPoint]; 
	\[Theta] = r[[2]]; \[Phi] = r[[3]]; {Cos[\[Theta]/2], Exp[I \[Phi]] Sin[\[Theta]/2]}]
BlochSphereToState[{\[Theta]_, \[Phi]_}] :=  {Cos[\[Theta]/2], Exp[I \[Phi]] Sin[\[Theta]/2]}
(* }}} *)
(* {{{  *) QuantumMutualInformationReducedEntropies[r_?MatrixQ] := Module[{ra, rb},
  ra = PartialTrace[r, 1]; rb = PartialTrace[r, 2];
  vonNeumannEntropy[ra] + vonNeumannEntropy[rb] -  vonNeumannEntropy[r]]
(* }}} *)
(* {{{  *) QuantumMutualInformationMeasurements[r_?MatrixQ] := 
 Module[{th, phi}, Maximize[ QuantumMutualInformationMeasurements[ r, {th, phi}], {th, phi}][[1]]]


(* {{{  *) QuantumMutualInformationMeasurements[ r_?MatrixQ, {\[Theta]_, \[CurlyPhi]_}] :=
(*Defined in eq 11 of Luo*)
 vonNeumannEntropy[PartialTrace[r, 2]] - 
  ConditionalEntropy[r, {\[Theta], \[CurlyPhi]}]
(* }}} *)
(* {{{  *) ConditionalEntropy[\[Rho]_, {\[Theta]_, \[CurlyPhi]_}] := 
 Module[{X, a, b, c, mm, mp, ps, \[Lambda]s},
  X = {2 Cos[\[Theta]] Sin[\[Theta]] Cos[\[CurlyPhi]], 
    2 Cos[\[Theta]] Sin[\[Theta]] Sin[\[CurlyPhi]], 
    2 Cos[\[Theta]]^2 - 1};
  {a, b, c} = BlochNormalFormVectors[\[Rho]];
  {mp = a + c X, mm = a - c X};
  ps = {(1 + b.X)/2, (1 - b.X)/2};
  \[Lambda]s = {{1/2 (1 + Norm[mp]/(1 + b.X)), 
     1/2 (1 - Norm[mp]/(1 + b.X))}, {1/2 (1 + Norm[mm]/(1 - b.X)), 
     1/2 (1 - Norm[mm]/(1 - b.X))}};
  ps[[1]] vonNeumannEntropy[\[Lambda]s[[1]]] + 
   ps[[2]] vonNeumannEntropy[\[Lambda]s[[2]]]
  ]
(* }}} *)
 
(* }}} *)
(* {{{  *) Discord[r_?MatrixQ] := QuantumMutualInformationReducedEntropies[r] - QuantumMutualInformationMeasurements[r]
(* }}} *)
(* {{{  *) ApplyGate[U_, State_, TargetQubits_] := Module[{StateOut, i, ie},
  StateOut = Table[Null, {Length[State]}];
  For[i = 0, i < Length[State]/2, i++,
   ie = MergeTwoIntegers[#, i, 
       Power[2, TargetQubits]] & /@ (Range[Length[U]] - 1);
   StateOut[[ie + 1]] = U.State[[ie + 1]];
   ];
  StateOut]
(* }}} *)
(* {{{  *) ApplyControlGate[U_, State_, TargetQubit_, ControlQubit_] := 
 Module[{StateOut, iwithcontrol, ie, NormalizedControlQubit},
  StateOut = State;
  If[ControlQubit > TargetQubit, 
   NormalizedControlQubit = ControlQubit - 1, 
   NormalizedControlQubit = ControlQubit];
  For[i = 0, i < Length[State]/4, i++,
   iwithcontrol = MergeTwoIntegers[1, i, Power[2, NormalizedControlQubit]];
   (*Print[{i,iwithcontrol}];*)
   
   ie = MergeTwoIntegers[#, iwithcontrol, 
       Power[2, TargetQubit]] & /@ {0, 1};
   StateOut[[ie + 1]] = U.State[[ie + 1]];
   ];
  StateOut]
(* }}} *)
(* {{{  *) BasisState[BasisNumber_,Dimension_] := Table[If[i == BasisNumber, 1, 0], {i, 0, Dimension - 1}]
(* }}} *)
(* }}} *)
(* {{{ Quantum channels and basis transformations *)
JamiolkowskiStateToOperatorChoi[Rho_?MatrixQ] := Sqrt[Length[Rho]] Reshuffle[Rho]
JamiolkowskiOperatorChoiToState[O_?MatrixQ] := Reshuffle[O]/Sqrt[Length[O]]
TransformationMatrixPauliBasisToComputationalBasis[] := {{1/Sqrt[2], 0, 0, 1/Sqrt[2]}, {0, 1/Sqrt[2], (-I)/Sqrt[2], 0}, {0, 1/Sqrt[2], I/Sqrt[2], 0}, 
 {1/Sqrt[2], 0, 0, -(1/Sqrt[2])}};

Reshuffle[Phi_?MatrixQ] := Module[{Dim, mn, MuNu, m, Mu, n, Nu},
   Dim = Sqrt[Length[Phi]];
   Table[ {m, n} = IntegerDigits[mn, Dim, 2];
    {Mu, Nu} = IntegerDigits[MuNu, Dim, 2];
    Phi[[FromDigits[{m, Mu}, Dim] + 1, 
     FromDigits[{n, Nu}, Dim] + 1]], {mn, 0, 
     Dim^2 - 1}, {MuNu , 0, Dim^2 - 1}]];


RandomTracePreservingMapChoiBasis[] := Module[{psi},
  psi = Total[
     MapThread[
      tensorProduct, {TwoRandomOrhtogonalStates[8], 
       TwoRandomOrhtogonalStates[2]}]]/Sqrt[2];
  Reshuffle[2 PartialTrace[Proyector[psi], 3]]
  ]
AveragePurityChannelPauliBasis[Lambda_] := Total[Power[Lambda[[All, 1]], 2]]/2 + Power[Norm[Lambda[[2 ;;, 2 ;;]], "Frobenius"], 2]/6



BlochEllipsoid[Cha_]:=Module[{center,T,coord,vecs,x,y,z,vecs0},
T=Cha[[{2,3,4},{2,3,4}]];
center={Cha[[2,1]],Cha[[3,1]],Cha[[4,1]]};
vecs0=Graphics3D[{{Dashing[0.01],Opacity[0.5],Red,Arrow[{center,1.3*Normalize[{1,0,0}]}],{Dashing[0.01],Opacity[0.5],Blue,Arrow[{center,1.3*Normalize[{0,1,0}]}]},{Dashing[0.01],Opacity[0.5],Green,Arrow[{center,1.3*Normalize[{0,0,1}]}]} } }];
vecs=Graphics3D[{{Red,Arrow[{center,1.3*Normalize[T.{1,0,0}]}],{Blue,Arrow[{center,1.3*Normalize[T.{0,1,0}]}]},{Green,Arrow[{center,1.3*Normalize[T.{0,0,1}]}]} }}];
coord={x,y,z}-center;
coord=Inverse[T].coord;
Show[ContourPlot3D[{coord[[1]]^2+coord[[2]]^2+coord[[3]]^2==1,x^2+y^2+z^2==1},{x,-1,1},{y,-1,1},{z,-1,1},AxesLabel->{"X","Y","Z"},ContourStyle->{Automatic,Opacity[0.3]},Mesh->None],vecs,vecs0,PlotRange->1.6]
]
(* {{{ *) EvolvGate[Gate_, steps_, env_, state_]:=
 Module[{statefinal, list, gate, j},
  statefinal = tensorProduct[env,state];
  gate[statefinal_] := First[Gate & /@ {statefinal}];
  list = Table[statefinal = gate[statefinal];
    PartialTrace[statefinal, 1], {j, steps}
    ];
  list
  ]
(*}}}*)
(* {{{ *) MakeQuantumChannel[Gate_, steps_, env_] := 
 Module[{\[Sigma], channel, Y, X, cero, uno},
  Y = EvolvGate[Gate, steps, env, {-I, 1}];
  cero = EvolvGate[Gate, steps, env, {0, 1}];
  uno = EvolvGate[Gate, steps, env, {1, 0}];
  X = EvolvGate[Gate, steps, env, {1, 1}];
  \[Sigma][2] = Y - cero - uno;
  \[Sigma][1] = X - cero - uno;
  \[Sigma][0] = cero + uno;
  \[Sigma][3] = uno - cero;
  Table[channel[i, j] = 
    1/2 Table[Tr[PauliMatrix[i].\[Sigma][j][[k]]], {k, steps}], {i, 0,
     3}, {j, 0, 3}];
  Chop[Table[
    Table[channel[i, j][[k]], {i, 0, 3}, {j, 0, 3}], {k, steps}]]]
(* }}} *)
(* {{{ *)SumPositiveDerivatives[list_]:=Module[{sum},
sum=0;
Table[
If[list[[i+1]]>list[[i]],
sum=sum+list[[i+1]]-list[[i]];
],
{i,1,Length[list]-1}];
sum]
(* }}} *)
(* {{{ *)
GHZ[qu_]:=1/Sqrt[2]Table[If[i==0||i==2^qu-1,1,0],{i,0,2^qu-1}]
(* }}} *)
(* }}} *)
End[] 
EndPackage[]
(* }}} *)
