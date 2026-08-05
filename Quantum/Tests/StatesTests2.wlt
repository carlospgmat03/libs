(* Tests for States.wl new additions *)

(* OrthonormalBasisContaningVector *)
VerificationTest[
  Module[{psi = {1, 0, 0, 0}, basis},
    basis = OrthonormalBasisContaningVector[psi];
    Length[basis]],
  4,
  TestID -> "OrthonormalBasisContaningVector-correct-size"]

(* First vector should be the input (up to global phase) *)
VerificationTest[
  Module[{psi = {1, 0, 0, 0}, basis},
    basis = OrthonormalBasisContaningVector[psi];
    Abs[Norm[basis[[1]] - psi]] < 10^-10],
  True,
  TestID -> "OrthonormalBasisContaningVector-first-is-psi"]

(* All basis vectors should be normalized *)
VerificationTest[
  Module[{psi = Normalize[{1, 1, 0, 0}], basis},
    basis = OrthonormalBasisContaningVector[psi];
    AllTrue[Abs[Norm[#] - 1] & /@ basis, # < 10^-10 &]],
  True,
  TestID -> "OrthonormalBasisContaningVector-all-normalized"]

(* Basis should be orthogonal *)
VerificationTest[
  Module[{psi = Normalize[{1, 1, 0, 0}], basis, gram},
    basis = OrthonormalBasisContaningVector[psi];
    gram = Table[Abs[Conjugate[basis[[i]]] . basis[[j]]], {i, 4}, {j, 4}];
    Max[Abs[gram - IdentityMatrix[4]]] < 10^-10],
  True,
  TestID -> "OrthonormalBasisContaningVector-orthonormal"]
