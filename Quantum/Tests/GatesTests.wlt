(* Tests for Gates.wl *)

(* Pauli[i] — single matrices *)
VerificationTest[Normal[Pauli[0]], {{1,0},{0,1}},
  TestID -> "Pauli-0-is-identity"]

VerificationTest[Normal[Pauli[1]], {{0,1},{1,0}},
  TestID -> "Pauli-1-is-X"]

VerificationTest[Normal[Pauli[2]], {{0,-I},{I,0}},
  TestID -> "Pauli-2-is-Y"]

VerificationTest[Normal[Pauli[3]], {{1,0},{0,-1}},
  TestID -> "Pauli-3-is-Z"]

(* Pauli[{i1,i2,...}] — tensor products *)
VerificationTest[
  Normal[Pauli[{1, 2}]],
  {{0,0,0,-I},{0,0,I,0},{0,-I,0,0},{I,0,0,0}},
  TestID -> "Pauli-tensor-X-otimes-Y"]

VerificationTest[
  Normal[Pauli[{3, 3}]],
  {{1,0,0,0},{0,-1,0,0},{0,0,-1,0},{0,0,0,1}},
  TestID -> "Pauli-tensor-Z-otimes-Z"]

(* Pauli matrices are Hermitian: Pauli[i] == ConjugateTranspose[Pauli[i]] *)
VerificationTest[
  Normal[Pauli[2]], Normal[ConjugateTranspose[Pauli[2]]],
  TestID -> "Pauli-2-is-hermitian"]

(* Paulib *)
VerificationTest[Normal[Paulib[{1,0,0}]], Normal[Pauli[1]],
  TestID -> "Paulib-unit-x-gives-PauliX"]

VerificationTest[Normal[Paulib[{0,0,1}]], Normal[Pauli[3]],
  TestID -> "Paulib-unit-z-gives-PauliZ"]

VerificationTest[Normal[Paulib[{0,0,0}]], {{0,0},{0,0}},
  TestID -> "Paulib-zero-vector-gives-zero-matrix"]
