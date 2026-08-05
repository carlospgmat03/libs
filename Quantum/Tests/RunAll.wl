Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "Kernel", "init.m"}]]

report = TestReport[FileNames["*.wlt", DirectoryName[$InputFileName]]]

Do[
  t = report["Results"][[i]];
  status = If[t["Outcome"] === "Success", "PASS", "FAIL"];
  Print[status, "  ", t["TestID"]];
  If[t["Outcome"] =!= "Success",
    Print["      expected: ", t["ExpectedOutput"]];
    Print["      got:      ", t["ActualOutput"]]],
  {i, Length[report["Results"]]}]

succeeded = Length[Select[report["Results"], #["Outcome"] === "Success" &]];
total     = Length[report["Results"]];
Print["\n", succeeded, "/", total, " tests passed"]

If[report["ReportSucceeded"], Quit[0], Quit[1]]
