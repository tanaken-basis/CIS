## CIS (Conditional Independence Statements)
## tanaken (Kentaro TANAKA, 2016.02-)
## Use this program at your own risk.

source("./CIS_sub.R")

n <- 4
known_CIS_abc_vec <- c("<a,b|cd>","<c,d|a>","<c,d|b>","<c,d|>")
unknown_CIS_abc_txt <- c("<c,d|ab>")

### testing ["<a,b|cd>","<c,d|a>","<c,d|b>","<c,d|>"] => "<c,d|ab>" ?
# new method
method <- 1 # positive probability function
test_by_TSTS_restriction_abc(n, known_CIS_abc_vec, unknown_CIS_abc_txt, method = method)

# new method
method <- 0 # not necessarily positive case
test_by_TSTS_restriction_abc(n, known_CIS_abc_vec, unknown_CIS_abc_txt, method = method)

# existing method
test_by_BHLS_abc(n, known_CIS_abc_vec, unknown_CIS_abc_txt)

### closures of ["<a,b|cd>","<c,d|a>","<c,d|b>","<c,d|>"]
CIS_five_types_of_closures(n, known_CIS_abc_vec)
