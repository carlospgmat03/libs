(* Tests for new Gates.wl additions *)

(* Hadamard *)
VerificationTest[
  Hadamard[],
  {{1, 1}, {1, -1}} / Sqrt[2],
  TestID -> "Hadamard-definition"]

VerificationTest[
  Chop[Hadamard[] . ConjugateTranspose[Hadamard[]] - IdentityMatrix[2]],
  {{0, 0}, {0, 0}},
  TestID -> "Hadamard-is-unitary"]

VerificationTest[
  Chop[Hadamard[] . Hadamard[] - IdentityMatrix[2]],
  {{0, 0}, {0, 0}},
  TestID -> "Hadamard-is-involution"]

(* LocalOneQubitOperator: 1-indexed from MSB *)
(* k=1 in 2-qubit system places gate on first (high) qubit: a \[CircleTimes] I *)
VerificationTest[
  Normal[LocalOneQubitOperator[2, 1, {{0,1},{1,0}}]],
  Normal[KroneckerProduct[{{0,1},{1,0}}, IdentityMatrix[2]]],
  TestID -> "LocalOneQubitOperator-k=1-is-gate-otimes-I"]

(* k=2 in 2-qubit system places gate on second (low) qubit: I \[CircleTimes] a *)
VerificationTest[
  Normal[LocalOneQubitOperator[2, 2, {{0,1},{1,0}}]],
  Normal[KroneckerProduct[IdentityMatrix[2], {{0,1},{1,0}}]],
  TestID -> "LocalOneQubitOperator-k=2-is-I-otimes-gate"]

VerificationTest[
  Dimensions[LocalOneQubitOperator[3, 2, {{0,1},{1,0}}]],
  {8, 8},
  TestID -> "LocalOneQubitOperator-3qubit-is-8x8"]

(* ControlNotMatrix: 0-indexed from LSB (consistent with ApplyMultiQubitGate) *)
(* control=qubit1(MSB), target=qubit0(LSB) in 2-qubit system *)
VerificationTest[
  Normal[ControlNotMatrix[1, 0, 2]],
  {{1,0,0,0},{0,1,0,0},{0,0,0,1},{0,0,1,0}},
  TestID -> "ControlNotMatrix-standard-CNOT-2qubit"]

(* CNOT on |10> gives |11> *)
VerificationTest[
  Normal[ControlNotMatrix[1, 0, 2]] . {0, 0, 1, 0},
  {0, 0, 0, 1},
  TestID -> "ControlNotMatrix-|10>->|11>"]

(* CNOT is its own inverse *)
VerificationTest[
  Chop[Normal[ControlNotMatrix[1, 0, 2]] . Normal[ControlNotMatrix[1, 0, 2]] - IdentityMatrix[4]],
  Table[0, {4}, {4}],
  TestID -> "ControlNotMatrix-is-involution"]

(* SWAPMatrix: 1-indexed from MSB *)
(* SWAP in 2-qubit: SWAPMatrix[2, {1,2}] *)
VerificationTest[
  Normal[SWAPMatrix[2, {1, 2}]],
  {{1,0,0,0},{0,0,1,0},{0,1,0,0},{0,0,0,1}},
  TestID -> "SWAPMatrix-2qubit"]

VerificationTest[
  SWAPMatrix[2, {1, 1}],
  IdentityMatrix[4, SparseArray],
  TestID -> "SWAPMatrix-same-qubit-is-identity"]

(* SumSigmaX/Y/Z *)
VerificationTest[
  Dimensions[SumSigmaX[3]], {8, 8},
  TestID -> "SumSigmaX-3qubit-is-8x8"]

(* SumSigmaX for 1 qubit should be PauliX *)
VerificationTest[
  SumSigmaX[1],
  {{0, 1}, {1, 0}},
  TestID -> "SumSigmaX-1qubit-is-PauliX"]

VerificationTest[
  SumSigmaZ[1],
  {{1, 0}, {0, -1}},
  TestID -> "SumSigmaZ-1qubit-is-PauliZ"]

(* For 2 qubits, total Sz: eigenvalues should be -2,-1,0,1 -> no, for 2 spin-1/2 with Sz=-1/2,1/2:
   total Sz ranges from -1 to 1. Specifically: {-2, 0, 0, 2} for Sum sigma_z (factor 2x from sigma not S) *)
VerificationTest[
  Sort[Eigenvalues[SumSigmaZ[2]]],
  {-2, 0, 0, 2},
  TestID -> "SumSigmaZ-2qubit-eigenvalues"]

(* ApplyGate: X on qubit 0 of |00> -> |01> *)
VerificationTest[
  ApplyGate[{{0,1},{1,0}}, {1,0,0,0}, 0],
  {0,1,0,0},
  TestID -> "ApplyGate-X-on-qubit0-of-|00>"]

(* ApplyGate: X on qubit 1 of |00> -> |10> *)
VerificationTest[
  ApplyGate[{{0,1},{1,0}}, {1,0,0,0}, 1],
  {0,0,1,0},
  TestID -> "ApplyGate-X-on-qubit1-of-|00>"]

(* ApplyControlGate: U on target=0, control=1: on |10> -> |11> *)
VerificationTest[
  ApplyControlGate[{{0,1},{1,0}}, {0,0,1,0}, 0, 1],
  {0,0,0,1},
  TestID -> "ApplyControlGate-X-controlled-on-qubit1-|10>->|11>"]

(* ApplyControlGate with control bit 0 leaves state unchanged *)
VerificationTest[
  ApplyControlGate[{{0,1},{1,0}}, {0,1,0,0}, 0, 1],
  {0,1,0,0},
  TestID -> "ApplyControlGate-no-action-when-control-is-0"]
