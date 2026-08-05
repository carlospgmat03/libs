(* Tests for Utilities.wl *)

(* BitsOn *)
VerificationTest[BitsOn[0], 0, TestID -> "BitsOn-zero"]
VerificationTest[BitsOn[5], 2, TestID -> "BitsOn-5=101-has-two-ones"]
VerificationTest[BitsOn[7], 3, TestID -> "BitsOn-7=111-has-three-ones"]

(* MergeTwoIntegers — example from the original usage comment:
   ndigits = 41 = 0 0 1 0 1 0 0 1
   na = 5        =     1   0     1  (fills the 1-bit positions)
   nb = 18       = 1 0   0   1 0    (fills the 0-bit positions)
   result        = 1 0 1 0 0 1 0 1 = 165 *)
VerificationTest[MergeTwoIntegers[5, 18, 41], 165,
  TestID -> "MergeTwoIntegers-usage-example-na=5-nb=18-ndigits=41"]

(* na occupies only the high bit (ndigits=4=100), nb the low two bits.
   na=1 -> bit 2 set; nb=3=11 -> bits 1,0 set. Result = 111 = 7 *)
VerificationTest[MergeTwoIntegers[1, 3, 4], 7,
  TestID -> "MergeTwoIntegers-na-in-high-bit-only"]

(* na occupies only the low bit (ndigits=1=001), nb the high two bits.
   na=0 -> bit 0 clear; nb=3=11 -> bits 2,1 set. Result = 110 = 6 *)
VerificationTest[MergeTwoIntegers[0, 3, 1], 6,
  TestID -> "MergeTwoIntegers-na-in-low-bit-only-na=0"]

(* na occupies only the low bit; na=1 -> bit 0 set; nb=0 -> bits 2,1 clear. Result = 001 = 1 *)
VerificationTest[MergeTwoIntegers[1, 0, 1], 1,
  TestID -> "MergeTwoIntegers-na-in-low-bit-only-na=1"]

(* ndigits = 7 = 111: all slots go to na, nb is ignored entirely. Result = na = 5 *)
VerificationTest[MergeTwoIntegers[5, 0, 7], 5,
  TestID -> "MergeTwoIntegers-all-bits-to-na"]

(* ndigits = 0: all slots go to nb, na is ignored entirely. Result = nb = 5 *)
VerificationTest[MergeTwoIntegers[0, 5, 0], 5,
  TestID -> "MergeTwoIntegers-all-bits-to-nb"]

(* ExtractDigits — inverse of MergeTwoIntegers:
   MergeTwoIntegers[5, 18, 41] = 165, so ExtractDigits[165, 41] = {18, 5} *)
VerificationTest[ExtractDigits[165, 41], {18, 5},
  TestID -> "ExtractDigits-inverse-of-MergeTwoIntegers"]

VerificationTest[ExtractDigits[0, 41], {0, 0},
  TestID -> "ExtractDigits-zero"]

(* tensorProduct: matrices — Z \[CircleTimes] I *)
VerificationTest[
  tensorProduct[{{1,0},{0,-1}}, {{1,0},{0,1}}],
  {{1,0,0,0},{0,1,0,0},{0,0,-1,0},{0,0,0,-1}},
  TestID -> "tensorProduct-matrices-Z-otimes-I"]

(* tensorProduct: vectors — |0>\[CircleTimes]|1> = |01> *)
VerificationTest[
  tensorProduct[{1,0}, {0,1}],
  {0,1,0,0},
  TestID -> "tensorProduct-vectors-|0>-otimes-|1>"]

(* tensorProduct[loc, A, B]: loc=1 (=01) puts B on the low bit, A on the high bit
   -> standard order A\[CircleTimes]B *)
VerificationTest[
  tensorProduct[1, {{1,0},{0,-1}}, {{0,1},{1,0}}],
  KroneckerProduct[{{1,0},{0,-1}}, {{0,1},{1,0}}],
  TestID -> "tensorProduct-loc=1-gives-standard-order"]

(* tensorProduct[loc, A, B]: loc=2 (=10) puts B on the high bit, A on the low bit
   -> swapped order B\[CircleTimes]A *)
VerificationTest[
  tensorProduct[2, {{1,0},{0,-1}}, {{0,1},{1,0}}],
  KroneckerProduct[{{0,1},{1,0}}, {{1,0},{0,-1}}],
  TestID -> "tensorProduct-loc=2-gives-swapped-order"]

(* tensorPower: identity matrix squares to itself *)
VerificationTest[
  tensorPower[{{1,0},{0,1}}, 2],
  {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}},
  TestID -> "tensorPower-identity-2-fold"]

(* tensorPower: n=1 returns the matrix unchanged *)
VerificationTest[
  tensorPower[{{0,1},{1,0}}, 1],
  {{0,1},{1,0}},
  TestID -> "tensorPower-n=1-is-identity"]

(* QuantumDotProduct *)
VerificationTest[QuantumDotProduct[{1,0},{0,1}], 0,
  TestID -> "QuantumDotProduct-orthogonal-states"]

VerificationTest[QuantumDotProduct[{1,0},{1,0}], 1,
  TestID -> "QuantumDotProduct-norm-squared"]

VerificationTest[QuantumDotProduct[{1,0},{I,0}], I,
  TestID -> "QuantumDotProduct-complex"]

(* PartialTrace: vector input — Bell state traced over qubit 1, keeping qubit 0 *)
VerificationTest[
  PartialTrace[{1,0,0,1}/Sqrt[2], 1],
  {{1/2,0},{0,1/2}},
  TestID -> "PartialTrace-vector-Bell-gives-maximally-mixed"]

(* PartialTrace: vector input — |01> traced over qubit 1, keeping qubit 0 -> |1><1| *)
VerificationTest[
  PartialTrace[{0,1,0,0}, 1],
  {{0,0},{0,1}},
  TestID -> "PartialTrace-vector-|01>-keep-qubit0"]

(* PartialTrace: matrix input — |01><01| traced over qubit 1, keeping qubit 0 -> |1><1| *)
VerificationTest[
  PartialTrace[{{0,0,0,0},{0,1,0,0},{0,0,0,0},{0,0,0,0}}, 1],
  {{0,0},{0,1}},
  TestID -> "PartialTrace-matrix-|01><01|-keep-qubit0"]

(* PartialTrace: matrix input — |01><01| traced over qubit 0, keeping qubit 1 -> |0><0| *)
VerificationTest[
  PartialTrace[{{0,0,0,0},{0,1,0,0},{0,0,0,0},{0,0,0,0}}, 2],
  {{1,0},{0,0}},
  TestID -> "PartialTrace-matrix-|01><01|-keep-qubit1"]

(* PartialTranspose: Bell state has negative PT eigenvalue (PPT criterion) *)
VerificationTest[
  Min[Re[Eigenvalues[PartialTranspose[Proyector[{1,0,0,1}/Sqrt[2]], 1]]]] < 0,
  True,
  TestID -> "PartialTranspose-Bell-fails-PPT"]

(* PartialTranspose: separable state |00><00| has non-negative PT eigenvalues *)
VerificationTest[
  Min[Re[Eigenvalues[PartialTranspose[Proyector[{1,0,0,0}], 1]]]] >= 0,
  True,
  TestID -> "PartialTranspose-separable-passes-PPT"]

(* PartialTransposeFirst/Second agree on Bell state *)
VerificationTest[
  Sort[Re[Eigenvalues[PartialTransposeFirst[Proyector[{1,0,0,1}/Sqrt[2]]]]]],
  {-1/2, 1/2, 1/2, 1/2},
  TestID -> "PartialTransposeFirst-Bell-eigenvalues"]

VerificationTest[
  Sort[Re[Eigenvalues[PartialTransposeSecond[Proyector[{1,0,0,1}/Sqrt[2]]]]]],
  {-1/2, 1/2, 1/2, 1/2},
  TestID -> "PartialTransposeSecond-Bell-eigenvalues"]

(* ApplyMultiQubitGate: vector state *)
(* X on single qubit |0> -> |1> *)
VerificationTest[
  ApplyMultiQubitGate[{1,0}, {{0,1},{1,0}}, 1],
  {0,1},
  TestID -> "ApplyMultiQubitGate-X-on-qubit0-of-|0>"]

(* X on qubit 0 (LSB) of |00> -> |01> *)
VerificationTest[
  ApplyMultiQubitGate[{1,0,0,0}, {{0,1},{1,0}}, 1],
  {0,1,0,0},
  TestID -> "ApplyMultiQubitGate-X-on-qubit0-of-|00>"]

(* X on qubit 1 of |00> -> |10> *)
VerificationTest[
  ApplyMultiQubitGate[{1,0,0,0}, {{0,1},{1,0}}, 2],
  {0,0,1,0},
  TestID -> "ApplyMultiQubitGate-X-on-qubit1-of-|00>"]

(* CNOT (control=qubit1, target=qubit0) on |10> -> |11> *)
VerificationTest[
  ApplyMultiQubitGate[{0,0,1,0}, {{1,0,0,0},{0,1,0,0},{0,0,0,1},{0,0,1,0}}, 3],
  {0,0,0,1},
  TestID -> "ApplyMultiQubitGate-CNOT-on-|10>"]

(* ApplyMultiQubitGate: density matrix state *)
(* X on |0><0| -> |1><1| *)
VerificationTest[
  ApplyMultiQubitGate[{{1,0},{0,0}}, {{0,1},{1,0}}, 1],
  {{0,0},{0,1}},
  TestID -> "ApplyMultiQubitGate-X-on-|0><0|-density-matrix"]

(* bad targets -> $Failed *)
VerificationTest[
  ApplyMultiQubitGate[{1,0,0,0}, {{0,1},{1,0}}, 3],
  $Failed,
  {ApplyMultiQubitGate::dimerr},
  TestID -> "ApplyMultiQubitGate-bad-targets-returns-Failed"]

(* ApplySwap: vector state, 0-indexed from LSB *)
(* swap qubits 0 and 1 of |01> -> |10> *)
VerificationTest[
  ApplySwap[{0,1,0,0}, 0, 1],
  {0,0,1,0},
  TestID -> "ApplySwap-vector-|01>-swap01->|10>"]

(* swap qubits 0 and 2 in 3-qubit system: |001> -> |100> *)
VerificationTest[
  ApplySwap[{0,1,0,0,0,0,0,0}, 0, 2],
  {0,0,0,0,1,0,0,0},
  TestID -> "ApplySwap-vector-3qubit-swap02"]

(* swap qubit with itself is identity *)
VerificationTest[
  ApplySwap[{1,2,3,4}, 0, 0],
  {1,2,3,4},
  TestID -> "ApplySwap-self-swap-is-identity"]

(* list form matches two-arg form *)
VerificationTest[
  ApplySwap[{0,1,0,0}, {0,1}],
  ApplySwap[{0,1,0,0}, 0, 1],
  TestID -> "ApplySwap-list-form-matches-two-arg"]

(* ApplySwap: density matrix |01><01| -> |10><10| after swap qubits 0,1 *)
VerificationTest[
  ApplySwap[{{0,0,0,0},{0,1,0,0},{0,0,0,0},{0,0,0,0}}, 0, 1],
  {{0,0,0,0},{0,0,0,0},{0,0,1,0},{0,0,0,0}},
  TestID -> "ApplySwap-density-matrix-|01><01|->|10><10|"]
