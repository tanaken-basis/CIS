## CIS (Conditional Independence Statements)
## tanaken (Kentaro TANAKA, 2016.02-)
## Use this program at your own risk.

library(lpSolve)
library(MASS)

####################
#### subfunc
## transform integer into binary vector
i_to_binvec <- function(i, d){ # i >= 0, d:digit
 if(i%%1 != 0){
  return(NULL)
 }
 binvec <- c()
 for(j in (d:1)){
  p <- 2**(j-1)
  q <- i %/% p
  r <- i %% p
  binvec <- c(q, binvec)
  i <- r
 }
 return(binvec)
}
## num --> abc
replace_num_to_abc <- function(txt){
  txt <- gsub("\\(0\\)", "a", txt)
  txt <- gsub("\\(1\\)", "b", txt)
  txt <- gsub("\\(2\\)", "c", txt)
  txt <- gsub("\\(3\\)", "d", txt)
  txt <- gsub("\\(4\\)", "e", txt)
  txt <- gsub("\\(5\\)", "f", txt)
  txt <- gsub("\\(6\\)", "g", txt)
  txt <- gsub("\\(7\\)", "h", txt)
  txt <- gsub("\\(8\\)", "i", txt)
  txt <- gsub("\\(9\\)", "j", txt)
  txt <- gsub("\\(10\\)", "k", txt)
  txt <- gsub("\\(11\\)", "l", txt)
  txt <- gsub("\\(12\\)", "m", txt)
  txt <- gsub("\\(13\\)", "n", txt)
  txt <- gsub("\\(14\\)", "o", txt)
  txt <- gsub("\\(15\\)", "p", txt)
  txt <- gsub("\\(16\\)", "q", txt)
  txt <- gsub("\\(17\\)", "r", txt)
  txt <- gsub("\\(18\\)", "s", txt)
  txt <- gsub("\\(19\\)", "t", txt)
  txt <- gsub("\\(20\\)", "u", txt)
  txt <- gsub("\\(21\\)", "v", txt)
  txt <- gsub("\\(22\\)", "w", txt)
  txt <- gsub("\\(23\\)", "x", txt)
  txt <- gsub("\\(24\\)", "y", txt)
  txt <- gsub("\\(25\\)", "z", txt)
  return(txt)
}
## abc --> num
replace_abc_to_num <- function(txt){
  txt <- gsub("a", "(0)", txt)
  txt <- gsub("b", "(1)", txt)
  txt <- gsub("c", "(2)", txt)
  txt <- gsub("d", "(3)", txt)
  txt <- gsub("e", "(4)", txt)
  txt <- gsub("f", "(5)", txt)
  txt <- gsub("g", "(6)", txt)
  txt <- gsub("h", "(7)", txt)
  txt <- gsub("i", "(8)", txt)
  txt <- gsub("j", "(9)", txt)
  txt <- gsub("k", "(10)", txt)
  txt <- gsub("l", "(11)", txt)
  txt <- gsub("m", "(12)", txt)
  txt <- gsub("n", "(13)", txt)
  txt <- gsub("o", "(14)", txt)
  txt <- gsub("p", "(15)", txt)
  txt <- gsub("q", "(16)", txt)
  txt <- gsub("r", "(17)", txt)
  txt <- gsub("s", "(18)", txt)
  txt <- gsub("t", "(19)", txt)
  txt <- gsub("u", "(20)", txt)
  txt <- gsub("v", "(21)", txt)
  txt <- gsub("w", "(22)", txt)
  txt <- gsub("x", "(23)", txt)
  txt <- gsub("y", "(24)", txt)
  txt <- gsub("z", "(25)", txt)
  return(txt)
}
replace_num_to_abc_matrix_colnum <- function(mat){
  colnames(mat) <- replace_num_to_abc(colnames(mat))
  rownames(mat) <- replace_num_to_abc(rownames(mat))
  return(mat)
}
## generate elmentary imsets
set_elementary_imset <- function(a2, b2, c2, n2){
  elementary_imset <- rep(0, n2+1)
  elementary_imset[c2+1] <- 1
  elementary_imset[a2+c2+1] <- -1
  elementary_imset[b2+c2+1] <- -1
  elementary_imset[a2+b2+c2+1] <- 1
  return(elementary_imset)
}
gen_elementary_imsets <- function(n1, n2){ # all of the elementary_imsets
  n1 <- n-1
  n2 <- (2**n)-1
  elementary_imsets_array <- c()
  elementary_imset_eq_vec <- c()
  elementary_imset_eq_C_vec <- c()
  elementary_imset_eq_C_vec <- c()
  for(c2 in (0:n2)){
    c2_binvec <- i_to_binvec(c2, n)
    elementary_imset_eq_C = ""
    for(c in (0:n1)){
      if(c2_binvec[c+1]==1){
        elementary_imset_eq_C = paste(elementary_imset_eq_C, "(", c, ")", sep="")
        next
      }
    }
    elementary_imset_eq_C_vec <- c(elementary_imset_eq_C_vec, elementary_imset_eq_C)
    for(a in (0:n1)){
      if(c2_binvec[a+1]==1){
        next
      }
      elementary_imset_eq_A = paste("(", a, ")", sep="")
      a2 <- 2**a
      if(a+1 <= n1){
        for(b in ((a+1):n1)){
          if(c2_binvec[b+1]==1){
            next
          }
          elementary_imset_eq_B = paste("(", b, ")", sep="")
          b2 <- 2**b
          elementary_imset <- set_elementary_imset(a2, b2, c2, n2)
          elementary_imsets_array <- cbind(elementary_imsets_array, elementary_imset)
          elementary_imset_eq <- paste("<", elementary_imset_eq_A, ",", elementary_imset_eq_B, "|", elementary_imset_eq_C, ">", sep="")
          elementary_imset_eq_vec <- c(elementary_imset_eq_vec, elementary_imset_eq)
        }
      }
    }
  }
  colnames(elementary_imsets_array) <- elementary_imset_eq_vec
  rownames(elementary_imsets_array) <- elementary_imset_eq_C_vec
  return(elementary_imsets_array)
}
## set operators for CIS
div_set_to_vec <- function(set_txt){
  set_txt <- gsub("\\(", "", set_txt)
  return(unlist(strsplit(set_txt, "\\)")))
}
sort_elements <- function(set_txt){ # (0)(3)(2)(1) -> (0)(1)(2)(3)
  set_vec <- div_set_to_vec(set_txt)
  if(any(duplicated(set_vec))){
    warning("tanaken_warning: duplicated!")
  }
  set_vec <- sort(set_vec)
  set_txt <- paste(set_vec, sep="", collapse=")(")
  if(set_txt!=""){
    set_txt <- paste("(", set_txt, ")", sep="")
  }
  return(set_txt)
}
extract_A_and_B_and_C_from_CIS <- function(CIS_txt){# <A, B| C> -> A_B_C_vec[1]=A, A_B_C_vec[2]=B, A_B_C_vec[3]=C
  CIS_txt <- gsub(">", "", CIS_txt)
  CIS_txt <- gsub("<", "", CIS_txt)
  CIS_txt <- gsub("[[:space:]]", "", CIS_txt)
  A_B_C_vec <- unlist(strsplit(CIS_txt, "\\|"))
  A_B_C_vec <- unlist(strsplit(A_B_C_vec, ","))
  if(length(A_B_C_vec) == 2){
    A_B_C_vec <- c(A_B_C_vec, "")
  }
  return(A_B_C_vec)
}
extract_AB_from_CIS <- function(CIS_txt){# <A, B| C> -> AB
  return(sort_elements(paste((extract_A_and_B_and_C_from_CIS(CIS_txt))[1:2], collapse="")))
}
## CIS_txt -> (AC, BC, C, ABC)
div_CIS_into_four_parts <- function(CIS_txt){ # <A, B| C> -> four_parts_vec[1]=AC, four_parts_vec[2]=BC, four_parts_vec[3]=C, four_parts_vec[4]=ABC
  four_parts_vec <- extract_A_and_B_and_C_from_CIS(CIS_txt)
  four_parts_vec <- c(four_parts_vec, paste(four_parts_vec, sep="", collapse="")) # ABC
  four_parts_vec[1] <- paste(four_parts_vec[c(1,3)], sep="", collapse="") # AC
  four_parts_vec[2] <- paste(four_parts_vec[c(2,3)], sep="", collapse="") # BC
  four_parts_vec <- apply(as.matrix(four_parts_vec), 1, sort_elements)
  return(four_parts_vec)
}
## CIS_txt -> semi_elementary_imsets
convert_CIS_into_semi_elementary_imset <- function(CIS_txt, n2, elementary_imsets_array_rownames){
  four_parts_vec <- div_CIS_into_four_parts(CIS_txt)
  semi_elementary_imset_vec <- rep(0, n2+1)
  semi_elementary_imset_vec[elementary_imsets_array_rownames==four_parts_vec[1]] = -1
  semi_elementary_imset_vec[elementary_imsets_array_rownames==four_parts_vec[2]] = -1
  semi_elementary_imset_vec[elementary_imsets_array_rownames==four_parts_vec[3]] = 1
  semi_elementary_imset_vec[elementary_imsets_array_rownames==four_parts_vec[4]] = 1
  semi_elementary_imset_vec <- as.matrix(semi_elementary_imset_vec)
  rownames(semi_elementary_imset_vec) <- elementary_imsets_array_rownames
  colnames(semi_elementary_imset_vec) <- CIS_txt
  return(semi_elementary_imset_vec)
}
## CIS_vec -> semi_elementary_imsets (array)
convert_CIS_vec_into_semi_elementary_imsets_array <- function(CIS_vec, n2, elementary_imsets_array_rownames){
  CIS_vec_length <- length(CIS_vec)
  semi_elementary_imsets_array <- c()
  for( i in (1:CIS_vec_length) ){
    semi_elementary_imsets_array <- cbind(semi_elementary_imsets_array, convert_CIS_into_semi_elementary_imset(CIS_vec[i], n2, elementary_imsets_array_rownames))
  }
  return(semi_elementary_imsets_array)
}
## CIS_txt -> semi_elementary_cut_imsets
convert_CIS_into_semi_elementary_cut_imset <- function(CIS_txt, n2, elementary_imsets_array_rownames){
  four_parts_vec <- div_CIS_into_four_parts(CIS_txt)
  semi_elementary_imset_vec <- rep(0, n2+1)
  semi_elementary_imset_vec[elementary_imsets_array_rownames==four_parts_vec[1]] = -1
  semi_elementary_imset_vec[elementary_imsets_array_rownames==four_parts_vec[2]] = -1
  semi_elementary_imset_vec[elementary_imsets_array_rownames==four_parts_vec[4]] = 1
  semi_elementary_imset_vec <- as.matrix(semi_elementary_imset_vec)
  rownames(semi_elementary_imset_vec) <- elementary_imsets_array_rownames
  colnames(semi_elementary_imset_vec) <- CIS_txt
  return(semi_elementary_imset_vec)
}
## CIS_vec -> semi_elementary_cut_imsets (array)
convert_CIS_vec_into_semi_elementary_cut_imsets_array <- function(CIS_vec, n2, elementary_imsets_array_rownames){
  CIS_vec_length <- length(CIS_vec)
  semi_elementary_imsets_array <- c()
  for( i in (1:CIS_vec_length) ){
    semi_elementary_imsets_array <- cbind(semi_elementary_imsets_array, convert_CIS_into_semi_elementary_cut_imset(CIS_vec[i], n2, elementary_imsets_array_rownames))
  }
  return(semi_elementary_imsets_array)
}
## set operators
check_intersection <- function(set1_txt, set2_txt){
  set1_vec <- div_set_to_vec(set1_txt)
  set2_vec <- div_set_to_vec(set2_txt)
  if(any(duplicated(c(set1_vec, set2_vec)))){
    return(TRUE)
  } else {
    return(FALSE)
  }
}
check_inclusion_vecs <- function(set1_vec, set2_vec){ # set1 < set2
  set1_vec_length <- length(set1_vec)
  for(i in (1:set1_vec_length)){
    if(!(any(set2_vec == set1_vec[i]))){
      return(FALSE)
    }
  }
  return(TRUE)
}
check_inclusion_txt <- function(set1_txt, set2_txt){ # set1 < set2
  set1_vec <- div_set_to_vec(set1_txt)
  set2_vec <- div_set_to_vec(set2_txt)
  set1_vec_length <- length(set1_vec)
  for(i in (1:set1_vec_length)){
    if(!(any(set2_vec == set1_vec[i]))){
      return(FALSE)
    }
  }
  return(TRUE)
}
## get extra CIS
get_extra_CIS_elementary_imsets_array <- function(unknown_CIS_txt, elementary_imsets_array, method=1){
  # method = 0({EFG} \subset ({AC} \cup {BC})
  # method = 1({EFG} \subset ({ACD} \cup {BCD})
  # method = 2({EF} \subset ({ACD} \cup {BCD})
  elementary_imsets_array_colnames <- colnames(elementary_imsets_array)
  elementary_imsets_array_ncol <- dim(elementary_imsets_array)[2]
  extra_CIS_elementary_imsets_array <- c()
  extra_CIS_elementary_imsets_array_colnames <- c()
  if(method == 0){
    unknown_CIS_four_parts_vec <- div_CIS_into_four_parts(unknown_CIS_txt)
    for(i in (1:elementary_imsets_array_ncol)){
      elementary_imset_four_parts_vec <- div_CIS_into_four_parts(elementary_imsets_array_colnames[i])
      if( (check_inclusion_txt(elementary_imset_four_parts_vec[4], unknown_CIS_four_parts_vec[1])) || (check_inclusion_txt(elementary_imset_four_parts_vec[4], unknown_CIS_four_parts_vec[2])) ){
        extra_CIS_elementary_imsets_array <- cbind(extra_CIS_elementary_imsets_array, elementary_imsets_array[,i])
        extra_CIS_elementary_imsets_array_colnames <- c(extra_CIS_elementary_imsets_array_colnames, elementary_imsets_array_colnames[i])
      }
    }
  } else {
    unknown_CIS_A_B_C_vec <- extract_A_and_B_and_C_from_CIS(unknown_CIS_txt)
    if(method == 2){
      for(i in (1:elementary_imsets_array_ncol)){
        candidate_vec <- extract_AB_from_CIS(elementary_imsets_array_colnames[i]) # {EF}
        if( !(check_intersection(unknown_CIS_A_B_C_vec[1], candidate_vec)) || !(check_intersection(unknown_CIS_A_B_C_vec[2], candidate_vec)) ){ # non-bridging
          extra_CIS_elementary_imsets_array <- cbind(extra_CIS_elementary_imsets_array, elementary_imsets_array[,i])
          extra_CIS_elementary_imsets_array_colnames <- c(extra_CIS_elementary_imsets_array_colnames, elementary_imsets_array_colnames[i])
        }
      }
    } else {
      for(i in (1:elementary_imsets_array_ncol)){
        candidate_vec <- div_CIS_into_four_parts(elementary_imsets_array_colnames[i])[4] # {EFG}
        if( !(check_intersection(unknown_CIS_A_B_C_vec[1], candidate_vec)) || !(check_intersection(unknown_CIS_A_B_C_vec[2], candidate_vec)) ){ # non-bridging
          extra_CIS_elementary_imsets_array <- cbind(extra_CIS_elementary_imsets_array, elementary_imsets_array[,i])
          extra_CIS_elementary_imsets_array_colnames <- c(extra_CIS_elementary_imsets_array_colnames, elementary_imsets_array_colnames[i])
        }
      }
    }
  }
  colnames(extra_CIS_elementary_imsets_array) <- extra_CIS_elementary_imsets_array_colnames
  return(extra_CIS_elementary_imsets_array)
}
## testing true or not
test_by_BHLS <- function(n2, known_CIS_semi_elementary_imsets_array, unknown_CIS_semi_elementary_imset_array, elementary_imsets_array, tol = 1e-8, known_CIS_vec=NULL){
  ## elementary imsets
  elementary_imsets_array_rownames <- rownames(elementary_imsets_array)
  elementary_imsets_array_colnames <- colnames(elementary_imsets_array)
  ## known and unknown
  if(is.null(known_CIS_vec)){
    known_CIS_vec <- colnames(known_CIS_semi_elementary_imsets_array)
  }
  known_CIS_vec <- replace_num_to_abc(known_CIS_vec)
  unknown_CIS_txt <- replace_num_to_abc(colnames(unknown_CIS_semi_elementary_imset_array))
  ## Amat xvec = bvec ?
  Amat <- cbind(known_CIS_semi_elementary_imsets_array, -elementary_imsets_array)
  bvec <- as.matrix(unknown_CIS_semi_elementary_imset_array)
  bvec_length <- dim(bvec)[1]
  for(i in (1:bvec_length)){
    if(bvec[i] < 0){
      bvec[i] <- -bvec[i]
      Amat[i,] <- -Amat[i,]
    }
  }
  ## two-phase simplex method
  ## min (0 | 1)^T (xvec^T | zvec^T)^T
  ## s.t. AAmat (xvec^T | zvec^T)^T = bvec, xvec >= 0, zvec >= 0
  nrow.Amat <- nrow(Amat)
  ncol.Amat <- ncol(Amat)
  cvec <- c( rep(0, ncol.Amat), rep(1, nrow.Amat) )
  Amat <- cbind(Amat, diag(nrow.Amat))
  res <- lp(direction = "min", objective.in=cvec, const.mat=Amat, const.dir=rep("=", nrow.Amat), const.rhs=bvec)
  opt <- res$objval
  sol <- as.matrix(res$solution)
  rownames(sol) <- (c(known_CIS_vec, paste("-", replace_num_to_abc(elementary_imsets_array_colnames), sep=""), paste("dummy", (1:nrow.Amat), sep="")))
  colnames(sol) <- (unknown_CIS_txt)
  solved <- FALSE
  if(abs(opt) < tol){
    solved <- TRUE
  }
  return(list(solved=solved, sol=sol, res=res))
}
test_by_BHLS_abc <- function(n, known_CIS_abc_vec, unknown_CIS_abc_txt, tol = 1e-8){ # ex) known_CIS_abc_vec = c(<a,b|c>,<a,c|>,...)
  n1 <- n-1
  n2 <- (2**n)-1
  elementary_imsets_array <- gen_elementary_imsets(n1, n2)
  elementary_imsets_array_rownames <- rownames(elementary_imsets_array)
  known_CIS_vec <- replace_abc_to_num(known_CIS_abc_vec)
  known_CIS_semi_elementary_imsets_array <- convert_CIS_vec_into_semi_elementary_imsets_array(known_CIS_vec, n2, elementary_imsets_array_rownames)
  unknown_CIS <- replace_abc_to_num(unknown_CIS_abc_txt)
  unknown_CIS_semi_elementary_imset_array <- convert_CIS_vec_into_semi_elementary_imsets_array(unknown_CIS, n2, elementary_imsets_array_rownames)
  res_test_by_BHLS <- test_by_BHLS(n2, known_CIS_semi_elementary_imsets_array, unknown_CIS_semi_elementary_imset_array, elementary_imsets_array, tol = 1e-8)
  return(res_test_by_BHLS)
}
test_by_TSTS_extra <- function(n2, known_CIS_semi_elementary_imsets_array, unknown_CIS_semi_elementary_imset_array, elementary_imsets_array, method = 1, tol = 1e-8, known_CIS_vec = NULL){
  ## elementary imsets
  elementary_imsets_array_rownames <- rownames(elementary_imsets_array)
  elementary_imsets_array_colnames <- colnames(elementary_imsets_array)
  ## known and unknown abc names
  if(is.null(known_CIS_vec)){
    known_CIS_vec <- colnames(known_CIS_semi_elementary_imsets_array)
  }
  known_CIS_vec <- replace_num_to_abc(known_CIS_vec)
  unknown_CIS_txt <- replace_num_to_abc(colnames(unknown_CIS_semi_elementary_imset_array))
  ## extra conditional independence statements
  extra_CIS_elementary_imsets_array <- get_extra_CIS_elementary_imsets_array(colnames(unknown_CIS_semi_elementary_imset_array), elementary_imsets_array, method)
  if(is.null(extra_CIS_elementary_imsets_array)){
    extra_CIS_elementary_imsets_array_colnames <- NULL
  } else{
    extra_CIS_elementary_imsets_array_colnames <- colnames(extra_CIS_elementary_imsets_array)
  }
  ## Amat xvec = bvec ?
  if(is.null(extra_CIS_elementary_imsets_array)){
   Amat <- cbind(known_CIS_semi_elementary_imsets_array, -known_CIS_semi_elementary_imsets_array)
  } else{
   Amat <- cbind(known_CIS_semi_elementary_imsets_array, extra_CIS_elementary_imsets_array, -known_CIS_semi_elementary_imsets_array, -extra_CIS_elementary_imsets_array)
  }
  bvec <- as.matrix(unknown_CIS_semi_elementary_imset_array)
  bvec_length <- dim(bvec)[1]
  for(i in (1:bvec_length)){
    if(bvec[i] < 0){
      bvec[i] <- -bvec[i]
      Amat[i,] <- -Amat[i,]
    }
  }
  ## two-phase simplex method
  ## min (0 | 1)^T (xvec^T | zvec^T)^T
  ## s.t. AAmat (xvec^T | zvec^T)^T = bvec, xvec >= 0, zvec >= 0
  nrow.Amat <- nrow(Amat)
  ncol.Amat <- ncol(Amat)
  cvec <- c( rep(0, ncol.Amat), rep(1, nrow.Amat) )
  Amat <- cbind(Amat, diag(nrow.Amat))
  res <- lp(direction = "min", objective.in=cvec, const.mat=Amat, const.dir=rep("=", nrow.Amat), const.rhs=bvec)
  opt <- res$objval
  sol <- as.matrix(res$solution)
  if(is.null(extra_CIS_elementary_imsets_array)){
    rownames(sol) <- (c(known_CIS_vec,
                        paste("-", known_CIS_vec, sep=""),
                        paste("dummy", (1:nrow.Amat), sep="")
                        )
                      )
  } else{
    rownames(sol) <- (c(known_CIS_vec,
                        paste("{", replace_num_to_abc(extra_CIS_elementary_imsets_array_colnames), "}", sep=""),
                        paste("-", known_CIS_vec, sep=""),
                        paste("-{", replace_num_to_abc(extra_CIS_elementary_imsets_array_colnames), "}", sep=""),
                        paste("dummy", (1:nrow.Amat), sep="")
                        )
                      )
  }
  colnames(sol) <- (unknown_CIS_txt)
  solved <- FALSE
  if(abs(opt) < tol){
    solved <- TRUE
  }
  return(list(solved=solved, sol=sol, res=res))
}
test_by_TSTS_extra_abc <- function(n, known_CIS_abc_vec, unknown_CIS_abc_txt, method = 1, tol = 1e-8){ # ex) known_CIS_abc_vec = c(<a,b|c>,<a,c|>,...)
  n1 <- n-1
  n2 <- (2**n)-1
  elementary_imsets_array <- gen_elementary_imsets(n1, n2)
  elementary_imsets_array_rownames <- rownames(elementary_imsets_array)
  known_CIS_vec <- replace_abc_to_num(known_CIS_abc_vec)
  known_CIS_semi_elementary_imsets_array <- convert_CIS_vec_into_semi_elementary_imsets_array(known_CIS_vec, n2, elementary_imsets_array_rownames)
  unknown_CIS <- replace_abc_to_num(unknown_CIS_abc_txt)
  unknown_CIS_semi_elementary_imset_array <- convert_CIS_vec_into_semi_elementary_imsets_array(unknown_CIS, n2, elementary_imsets_array_rownames)
  res_test_by_TSTS_extra <- test_by_TSTS_extra(n2, known_CIS_semi_elementary_imsets_array, unknown_CIS_semi_elementary_imset_array, elementary_imsets_array, method = method, tol = 1e-8)
  return(res_test_by_TSTS_extra)
}
get_indispensable_vec <- function(unknown_CIS_txt, elementary_imsets_array, method=1){
  # method = 0( ABCD \setminus(AC \cup BC) )
  # method = 1( ABCD \setminus(ACD \cup BCD) )
  nrow.elementary_imsets_array <- nrow(elementary_imsets_array)
  elementary_imsets_array_rownames <- rownames(elementary_imsets_array)
  dispensable_vec <- c(1)
  AC_vec <- div_set_to_vec(div_CIS_into_four_parts(unknown_CIS_txt)[1])
  BC_vec <- div_set_to_vec(div_CIS_into_four_parts(unknown_CIS_txt)[2])
  ABC_vec <- div_set_to_vec(div_CIS_into_four_parts(unknown_CIS_txt)[4])
  if(method == 0){
    for(i in (2:nrow.elementary_imsets_array)){
      temp_vec <- div_set_to_vec(elementary_imsets_array_rownames[i])
      if(check_inclusion_vecs(temp_vec, AC_vec) || check_inclusion_vecs(temp_vec, BC_vec)){
        dispensable_vec <- c(dispensable_vec, i)
      }
    }
  } else if(method == 1){
    D_vec <- as.character(((0:(n-1)))[is.na(match(as.character((0:(n-1))), ABC_vec))])
    ACD_vec <- c(AC_vec, D_vec)
    BCD_vec <- c(BC_vec, D_vec)
    for(i in (2:nrow.elementary_imsets_array)){
      temp_vec <- div_set_to_vec(elementary_imsets_array_rownames[i])
      if(check_inclusion_vecs(temp_vec, ACD_vec) || check_inclusion_vecs(temp_vec, BCD_vec)){
        dispensable_vec <- c(dispensable_vec, i)
      }
    }
  } else{
     ## ??
  }
  indispensable_vec <- (1:nrow.elementary_imsets_array)[-dispensable_vec]
  return(indispensable_vec)
}
test_by_TSTS_restriction <- function(n2, known_CIS_semi_elementary_imsets_array, unknown_CIS_semi_elementary_imset_array, elementary_imsets_array, method = method, tol = 1e-8, known_CIS_vec = NULL, unknown_CIS_txt = NULL){
  ## elementary imsets
  elementary_imsets_array_rownames <- rownames(elementary_imsets_array)
  elementary_imsets_array_colnames <- colnames(elementary_imsets_array)
  ## known and unknown CISs
  if(is.null(known_CIS_vec)){
    known_CIS_vec <- colnames(known_CIS_semi_elementary_imsets_array)
  }
  if(is.null(unknown_CIS_txt)){
    unknown_CIS_txt <- colnames(unknown_CIS_semi_elementary_imset_array)
  }
  # indispensable vec
  indispensable_vec <- get_indispensable_vec(unknown_CIS_txt, elementary_imsets_array, method)
  colnames_known_CIS_semi_elementary_imsets_array <- colnames(known_CIS_semi_elementary_imsets_array)
  known_CIS_semi_elementary_imsets_array <- as.matrix(known_CIS_semi_elementary_imsets_array[indispensable_vec,])
  # print(known_CIS_semi_elementary_imsets_array)
  colnames(known_CIS_semi_elementary_imsets_array) <- colnames_known_CIS_semi_elementary_imsets_array
  colnames_unknown_CIS_semi_elementary_imset_array_fst_col <- colnames(unknown_CIS_semi_elementary_imset_array)[1]
  unknown_CIS_semi_elementary_imset_array <- as.matrix(unknown_CIS_semi_elementary_imset_array[indispensable_vec,])
  colnames(unknown_CIS_semi_elementary_imset_array) <- colnames_unknown_CIS_semi_elementary_imset_array_fst_col
  elementary_imsets_array <- elementary_imsets_array[indispensable_vec,]
  ## known and unknown abc names
  known_CIS_abc_vec <- replace_num_to_abc(known_CIS_vec)
  unknown_CIS_abc_txt <- replace_num_to_abc(unknown_CIS_txt)
  ## the case for x with p(x;ABC)>0
  ## Amat xvec = bvec ?
  Amat <- cbind(known_CIS_semi_elementary_imsets_array, -known_CIS_semi_elementary_imsets_array)
  bvec <- as.matrix(unknown_CIS_semi_elementary_imset_array) # \delta_{ABC}
  ## two-phase simplex method
  ## min (0 | 1)^T (xvec^T | zvec^T)^T
  ## s.t. Amat (xvec^T | zvec^T)^T = bvec, xvec >= 0, zvec >= 0
  nrow.Amat <- nrow(Amat)
  ncol.Amat <- ncol(Amat)
  cvec <- c( rep(0, ncol.Amat), rep(1, nrow.Amat) )
  Amat <- cbind(Amat, diag(nrow.Amat))
  res <- lp(direction = "min", objective.in=cvec, const.mat=Amat, const.dir=rep("=", nrow.Amat), const.rhs=bvec)
  opt <- res$objval
  sol <- as.matrix(res$solution)
  rownames(sol) <- (c(known_CIS_abc_vec, paste("-", known_CIS_abc_vec, sep=""), paste("dummy", (1:nrow.Amat), sep="")))
  colnames(sol) <- (unknown_CIS_abc_txt)
  solved <- FALSE
  if(abs(opt) < tol){
    solved <- TRUE
  }
  if(method == 0){ # the case for x with p(x;ABC)=0
    known_CIS_semi_elementary_cut_imsets_array <- convert_CIS_vec_into_semi_elementary_cut_imsets_array(known_CIS_vec, n2, elementary_imsets_array_rownames)
    known_CIS_semi_elementary_cut_imsets_array <- as.matrix(known_CIS_semi_elementary_cut_imsets_array[indispensable_vec,])
    colnames(known_CIS_semi_elementary_cut_imsets_array) <- colnames_known_CIS_semi_elementary_imsets_array
    Amat <- known_CIS_semi_elementary_cut_imsets_array
    nrow.Amat <- nrow(Amat)
    ncol.Amat <- ncol(Amat)
    cvec <- c( rep(0, ncol.Amat), rep(1, nrow.Amat) )
    Amat <- cbind(Amat, diag(nrow.Amat))
    bvec <- as.matrix(unknown_CIS_semi_elementary_imset_array) # \delta_{ABC}
    res.zero <- lp(direction = "min", objective.in=cvec, const.mat=Amat, const.dir=rep(">=", nrow.Amat), const.rhs=bvec)
    opt.zero <- res.zero$objval
    sol.zero <- as.matrix(res.zero$solution)
    rownames(sol.zero) <- (c(known_CIS_abc_vec, paste("dummy", (1:nrow.Amat), sep="")))
    colnames(sol.zero) <- (unknown_CIS_abc_txt)
    solved.zero <- FALSE
    if(abs(opt.zero) < tol){
      solved.zero <- TRUE
    }
    return(list(solved=all(solved, solved.zero), solved.positive=solved, sol.positive=sol, res.positive=res, solved.zero=solved.zero, sol.zero=sol.zero, res.zero=res.zero, indispensable_rows=(replace_num_to_abc(elementary_imsets_array_rownames))[indispensable_vec], indispensable_vec=indispensable_vec))
  } else{
    return(list(solved=solved, sol=sol, res=res, indispensable_rows=(replace_num_to_abc(elementary_imsets_array_rownames))[indispensable_vec], indispensable_vec=indispensable_vec))
  }
}
test_by_TSTS_restriction_abc <- function(n, known_CIS_abc_vec, unknown_CIS_abc_txt, method = 1, tol = 1e-8){ # ex) known_CIS_abc_vec = c(<a,b|c>,<a,c|>,...)
  n1 <- n-1
  n2 <- (2**n)-1
  elementary_imsets_array <- gen_elementary_imsets(n1, n2)
  elementary_imsets_array_rownames <- rownames(elementary_imsets_array)
  known_CIS_vec <- replace_abc_to_num(known_CIS_abc_vec)
  known_CIS_semi_elementary_imsets_array <- convert_CIS_vec_into_semi_elementary_imsets_array(known_CIS_vec, n2, elementary_imsets_array_rownames)
  unknown_CIS_txt <- replace_abc_to_num(unknown_CIS_abc_txt)
  unknown_CIS_semi_elementary_imset_array <- convert_CIS_vec_into_semi_elementary_imsets_array(unknown_CIS_txt, n2, elementary_imsets_array_rownames)
  res_test_by_TSTS_restriction <- test_by_TSTS_restriction(n2, known_CIS_semi_elementary_imsets_array, unknown_CIS_semi_elementary_imset_array, elementary_imsets_array, method = method, tol = 1e-8)
  return(res_test_by_TSTS_restriction)
}
get_active_CIS_mat <- function(CIS_mat){
  active_CIS_mat <- (CIS_mat!=0)
  active_CIS_mat_colnames <- colnames(CIS_mat)
  active_CIS_mat_rownames <- rownames(CIS_mat)[active_CIS_mat]
  active_CIS_mat <- as.matrix((CIS_mat)[active_CIS_mat])
  colnames(active_CIS_mat) <- active_CIS_mat_colnames
  rownames(active_CIS_mat) <- active_CIS_mat_rownames
  return(active_CIS_mat)
}
##
get_tex_format_from_four_parts <- function(CIS_vec, m=1, tol = 1e-8){ # m: power.
  if(CIS_vec[3]==""){
    CIS_vec[3] <- paste("\\emptyset", sep="")
  }
  if(abs(m-1) < (tol*100)){
    CIS_tex <- paste("p(", CIS_vec, ")", sep="")
  } else{
    CIS_tex <- paste("p(", CIS_vec, ")^{", m, "}", sep="")
  }
  CIS_tex <- paste("{\\frac{", CIS_tex[4], CIS_tex[3], "}{", CIS_tex[1], CIS_tex[2], "}}", sep="")
  return(CIS_tex)
}
get_int_coef_CIS_mat <- function(CIS_mat, m=0, m_max=103, tol = 1e-8){
  if(m <= 0){
    for(i in (1:m_max)){
      remainder <- abs((i*CIS_mat) %% 1)
      gap <- abs(1-remainder)
      if( (max(apply(cbind(remainder, gap), 1, min))) < (tol*100) ){
        m <- i
        break
      }
    }
  }
  if(m == 0){
   m=1
  }
  CIS_mat <- m * CIS_mat
  return(list(mat=CIS_mat, m=m))
}
get_tex_format_from_CIS_mat <- function(CIS_mat, m=0, m_max=103, tol = 1e-8){ # CIS_mat: 1 column, rownames, colnames. m: power.
  CIS_mat_rownames <- replace_abc_to_num(rownames(CIS_mat))
  CIS_mat_colnames <- replace_abc_to_num(colnames(CIS_mat))
  left_hand_side_known_tex_vec <- c()
  left_hand_side_extra_tex_vec <- c()
  right_hand_side_known_tex_vec <- c()
  right_hand_side_extra_tex_vec <- c()
  CIS_mat <- get_int_coef_CIS_mat(CIS_mat, m, m_max, tol)$mat
  unknown_CIS_txt <- div_CIS_into_four_parts(CIS_mat_colnames)
  unknown_CIS_tex <- replace_num_to_abc(get_tex_format_from_four_parts(unknown_CIS_txt, m))
  CIS_mat_length <- length(CIS_mat)
  for(i in (1:CIS_mat_length)){
    CIS_txt <- CIS_mat_rownames[i]
    if(is.na(charmatch("-", CIS_txt))){
      if(is.na(charmatch("{", CIS_txt))){
        CIS_vec <- div_CIS_into_four_parts(CIS_txt)
        CIS_tex <- replace_num_to_abc(get_tex_format_from_four_parts(CIS_vec, m*CIS_mat[i]))
        left_hand_side_known_tex_vec <- c(left_hand_side_known_tex_vec, CIS_tex)
      } else{
        CIS_txt <- gsub("\\{", "", CIS_txt)
        CIS_txt <- gsub("\\}", "", CIS_txt)
        CIS_vec <- div_CIS_into_four_parts(CIS_txt)
        CIS_tex <- replace_num_to_abc(get_tex_format_from_four_parts(CIS_vec, m*CIS_mat[i]))
        left_hand_side_extra_tex_vec <- c(left_hand_side_extra_tex_vec, CIS_tex)
      }
    } else{
      CIS_txt <- gsub("-", "", CIS_txt)
      if(is.na(charmatch("{", CIS_txt))){
        CIS_vec <- div_CIS_into_four_parts(CIS_txt)
        CIS_tex <- replace_num_to_abc(get_tex_format_from_four_parts(CIS_vec, m*CIS_mat[i]))
        right_hand_side_known_tex_vec <- c(right_hand_side_known_tex_vec, CIS_tex)
      } else{
        CIS_txt <- gsub("\\{", "", CIS_txt)
        CIS_txt <- gsub("\\}", "", CIS_txt)
        CIS_vec <- div_CIS_into_four_parts(CIS_txt)
        CIS_tex <- replace_num_to_abc(get_tex_format_from_four_parts(CIS_vec, m*CIS_mat[i]))
        right_hand_side_extra_tex_vec <- c(right_hand_side_extra_tex_vec, CIS_tex)
      }
    }
  }
  left_hand_side_known_tex <- paste(left_hand_side_known_tex_vec, collapse="")
  left_hand_side_known_tex <- paste(left_hand_side_known_tex, "\n", sep="")
  if(length(left_hand_side_extra_tex_vec)==0){
    left_hand_side_extra_tex <- ""
  } else{
    left_hand_side_extra_tex <- paste(left_hand_side_extra_tex_vec, collapse="")
    left_hand_side_extra_tex <- paste("\\left\\{", left_hand_side_extra_tex, "\\right\\}", "\n", sep="")
  }
  right_hand_side_known_tex <- ""
  right_hand_side_extra_tex <- ""
  if(length(right_hand_side_known_tex_vec)==0){
    right_hand_side_known_tex <- ""
  } else{
    right_hand_side_known_tex <- paste(right_hand_side_known_tex_vec, collapse="")
    right_hand_side_known_tex <- paste(right_hand_side_known_tex, "\n", sep="")
  }
  if(length(right_hand_side_extra_tex_vec)==0){
    right_hand_side_extra_tex <- ""
  } else{
    right_hand_side_extra_tex <- paste(right_hand_side_extra_tex_vec, collapse="")
    right_hand_side_extra_tex <- paste("\\left\\{", right_hand_side_extra_tex, "\\right\\}", "\n", sep="")
  }
  CIS_tex <- ""
  CIS_tex <- paste(CIS_tex, "\\begin{eqnarray*}\n\\lefteqn{\n", sep="")
  CIS_tex <- paste(CIS_tex, left_hand_side_known_tex, sep="")
  CIS_tex <- paste(CIS_tex, left_hand_side_extra_tex, sep="")
  CIS_tex <- paste(CIS_tex, "} & & \\\\ & = & \n", sep="")
  CIS_tex <- paste(CIS_tex, "\\left[", unknown_CIS_tex, "\\right]\n", sep="")
  CIS_tex <- paste(CIS_tex, right_hand_side_known_tex, sep="")
  CIS_tex <- paste(CIS_tex, right_hand_side_extra_tex, sep="")
  CIS_tex <- paste(CIS_tex, "\\end{eqnarray*}\n", sep="")
  return(CIS_tex)
}
convert_tri <- function(CIS_txt){
  CIS_txt <- gsub("-", "", CIS_txt)
  CIS_txt <- gsub("\\{", "", CIS_txt)
  CIS_txt <- gsub("\\}", "", CIS_txt)
  CIS_txt <- gsub("<", "{", CIS_txt)
  CIS_txt <- gsub(">", "}", CIS_txt)
  CIS_txt <- gsub(",", "}{", CIS_txt)
  CIS_txt <- gsub("\\|", "}{", CIS_txt)
  return(CIS_txt)
}
not_one <- function(num){
  txt <- ""
  num <- round(num)
  if(num != 1){
    txt <- paste(num)
  }
  return(txt)
}
get_imset_relation_tex_format_from_CIS_mat <- function(CIS_mat, m=0, m_max=103, tol = 1e-8){ # CIS_mat: 1 column, rownames, colnames. m: power.
  CIS_mat <- get_active_CIS_mat(CIS_mat)
  CIS_mat_rownames <- replace_num_to_abc(rownames(CIS_mat))
  CIS_mat_colnames <- replace_num_to_abc(colnames(CIS_mat))
  CIS_mat_list <- get_int_coef_CIS_mat(CIS_mat, m, m_max, tol)
  CIS_mat <- CIS_mat_list$mat
  m <- CIS_mat_list$m
  CIS_mat_length <- length(CIS_mat)
  left_hand_side_known_tex_vec <- c()
  left_hand_side_extra_tex_vec <- c()
  right_hand_side_known_tex_vec <- c()
  right_hand_side_extra_tex_vec <- c()
  unknown_CIS_tex <- paste(not_one(m), "u_{\\tri", convert_tri(CIS_mat_colnames[1]), "}", sep="", collapse="")
  for(i in (1:CIS_mat_length)){
    if(is.na(charmatch("-", CIS_mat_rownames[i]))){
      if(is.na(charmatch("{", CIS_mat_rownames[i]))){
        left_hand_side_known_tex_vec <- c(left_hand_side_known_tex_vec, paste(not_one(CIS_mat[i]), "u_{\\tri", convert_tri(CIS_mat_rownames[i]), "}", sep="", collapse=""))
      } else{
        left_hand_side_extra_tex_vec <- c(left_hand_side_extra_tex_vec, paste(not_one(CIS_mat[i]), "u_{\\tri", convert_tri(CIS_mat_rownames[i]), "}", sep="", collapse=""))
      }
    } else{
      if(is.na(charmatch("-{", CIS_mat_rownames[i]))){
        right_hand_side_known_tex_vec <- c(right_hand_side_known_tex_vec, paste(not_one(CIS_mat[i]), "u_{\\tri", convert_tri(CIS_mat_rownames[i]), "}", sep="", collapse=""))
      } else{
        right_hand_side_extra_tex_vec <- c(right_hand_side_extra_tex_vec, paste(not_one(CIS_mat[i]), "u_{\\tri", convert_tri(CIS_mat_rownames[i]), "}", sep="", collapse=""))
      }
    }
  }
  left_hand_side_known_tex <- paste(left_hand_side_known_tex_vec, collapse=" + ")
  left_hand_side_known_tex <- paste(left_hand_side_known_tex, "\n", sep="")
  if(length(left_hand_side_extra_tex_vec)==0){
    left_hand_side_extra_tex <- ""
  } else{
    left_hand_side_extra_tex <- paste(left_hand_side_extra_tex_vec, collapse=" + ")
    left_hand_side_extra_tex <- paste("+ \\left\\{", left_hand_side_extra_tex, "\\right\\}", "\n", sep="")
  }
  right_hand_side_known_tex <- ""
  right_hand_side_extra_tex <- ""
  if(length(right_hand_side_known_tex_vec)==0){
    right_hand_side_known_tex <- ""
  } else{
    right_hand_side_known_tex <- paste(right_hand_side_known_tex_vec, collapse=" + ")
    right_hand_side_known_tex <- paste(" + ", right_hand_side_known_tex, "\n", sep="")
  }
  if(length(right_hand_side_extra_tex_vec)==0){
    right_hand_side_extra_tex <- ""
  } else{
    right_hand_side_extra_tex <- paste(right_hand_side_extra_tex_vec, collapse=" + ")
    right_hand_side_extra_tex <- paste(" + \\left\\{", right_hand_side_extra_tex, "\\right\\}", "\n", sep="")
  }
  CIS_tex <- ""
  CIS_tex <- paste(CIS_tex, "% \\usepackage{bm}\n", sep="")
  CIS_tex <- paste(CIS_tex, "% \\newcommand{\\Bu}{\\bm{u}}\n", sep="")
  CIS_tex <- paste(CIS_tex, "% \\newcommand{\\tri}[3]{\\langle #1, #2 \\, | \\, #3 \\rangle}\n", sep="")
  CIS_tex <- paste(CIS_tex, "\\begin{eqnarray*}\n\\lefteqn{\n", sep="")
  CIS_tex <- paste(CIS_tex, left_hand_side_known_tex, sep="")
  CIS_tex <- paste(CIS_tex, left_hand_side_extra_tex, sep="")
  CIS_tex <- paste(CIS_tex, "} & & \\\\ & = & \n", sep="")
  CIS_tex <- paste(CIS_tex, "\\left[", unknown_CIS_tex, "\\right]\n", sep="")
  CIS_tex <- paste(CIS_tex, right_hand_side_known_tex, sep="")
  CIS_tex <- paste(CIS_tex, right_hand_side_extra_tex, sep="")
  CIS_tex <- paste(CIS_tex, "\\end{eqnarray*}\n", sep="")
  return(CIS_tex)
}
get_indep_relation_tex_format_from_CIS_mat <- function(CIS_mat, m=0, m_max=103, tol = 1e-8){ # CIS_mat: 1 column, rownames, colnames. m: power.
  CIS_mat <- get_active_CIS_mat(CIS_mat)
  CIS_mat_rownames <- replace_num_to_abc(rownames(CIS_mat))
  CIS_mat_colnames <- replace_num_to_abc(colnames(CIS_mat))
  CIS_mat_list <- get_int_coef_CIS_mat(CIS_mat, m, m_max, tol)
  CIS_mat <- CIS_mat_list$mat
  m <- CIS_mat_list$m
  CIS_mat_length <- length(CIS_mat)
  left_hand_side_known_tex_vec <- c()
  left_hand_side_extra_tex_vec <- c()
  right_hand_side_known_tex_vec <- c()
  right_hand_side_extra_tex_vec <- c()
  unknown_CIS_tex <- paste("{\\indtris", convert_tri(CIS_mat_colnames[1]), "}", sep="", collapse="")
  for(i in (1:CIS_mat_length)){
    if(is.na(charmatch("-", CIS_mat_rownames[i]))){
      if(is.na(charmatch("{", CIS_mat_rownames[i]))){
        left_hand_side_known_tex_vec <- c(left_hand_side_known_tex_vec, paste("{\\indtris", convert_tri(CIS_mat_rownames[i]), "}", sep="", collapse=""))
      } else{
        left_hand_side_extra_tex_vec <- c(left_hand_side_extra_tex_vec, paste("{\\indtris", convert_tri(CIS_mat_rownames[i]), "}", sep="", collapse=""))
      }
    } else{
      if(is.na(charmatch("-{", CIS_mat_rownames[i]))){
        right_hand_side_known_tex_vec <- c(right_hand_side_known_tex_vec, paste("{\\indtris", convert_tri(CIS_mat_rownames[i]), "}", sep="", collapse=""))
      } else{
        right_hand_side_extra_tex_vec <- c(right_hand_side_extra_tex_vec, paste("{\\indtris", convert_tri(CIS_mat_rownames[i]), "}", sep="", collapse=""))
      }
    }
  }
  left_hand_side_known_tex <- paste(left_hand_side_known_tex_vec, collapse=" ,\\, ")
  left_hand_side_known_tex <- paste(left_hand_side_known_tex, "\n", sep="")
  if(length(left_hand_side_extra_tex_vec)==0){
    left_hand_side_extra_tex <- ""
  } else{
    left_hand_side_extra_tex <- paste(left_hand_side_extra_tex_vec, collapse=" ,\\, ")
    left_hand_side_extra_tex <- paste(",\\, ", left_hand_side_extra_tex, "\n", sep="") # paste(",\\, \\left\\{", left_hand_side_extra_tex, "\\right\\}", "\n", sep="")
  }
  right_hand_side_known_tex <- ""
  right_hand_side_extra_tex <- ""
  if(length(right_hand_side_known_tex_vec)==0){
    right_hand_side_known_tex <- ""
  } else{
    right_hand_side_known_tex <- paste(right_hand_side_known_tex_vec, collapse=" ,\\, ")
    right_hand_side_known_tex <- paste(" ,\\, ", right_hand_side_known_tex, "\n", sep="")
  }
  if(length(right_hand_side_extra_tex_vec)==0){
    right_hand_side_extra_tex <- ""
  } else{
    right_hand_side_extra_tex <- paste(right_hand_side_extra_tex_vec, collapse=" ,\\, ")
    right_hand_side_extra_tex <- paste(" ,\\, ", right_hand_side_extra_tex, "\n", sep="") # paste(" ,\\, \\left\\{", right_hand_side_extra_tex, "\\right\\}", "\n", sep="")
  }
  CIS_tex <- ""
  CIS_tex <- paste(CIS_tex, "% \\newcommand{\\indep}{\\mathop{\\perp\\!\\!\\!\\perp}}\n", sep="")
  CIS_tex <- paste(CIS_tex, "% \\newcommand{\\indtris}[3]{#1 \\indep #2 \\, | \\, #3  }\n", sep="")
  CIS_tex <- paste(CIS_tex, "\\begin{eqnarray*}\n\\lefteqn{\n", sep="")
  CIS_tex <- paste(CIS_tex, left_hand_side_known_tex, sep="")
  CIS_tex <- paste(CIS_tex, left_hand_side_extra_tex, sep="")
  CIS_tex <- paste(CIS_tex, "} & & \\\\ & \\Leftrightarrow & \n", sep="")
  CIS_tex <- paste(CIS_tex, unknown_CIS_tex, "\n", sep="") # "\\left[", unknown_CIS_tex, "\\right]\n", sep="")
  CIS_tex <- paste(CIS_tex, right_hand_side_known_tex, sep="")
  CIS_tex <- paste(CIS_tex, right_hand_side_extra_tex, sep="")
  CIS_tex <- paste(CIS_tex, "\\end{eqnarray*}\n", sep="")
  return(CIS_tex)
}
CIS.test.extra <- function(n, known, unknown, method, prncat=FALSE){
  n <- n
  known_CIS_abc_vec <- known
  unknown_CIS_abc_txt <- unknown
  method <- method
  ####################
  ## method = 0({EFG} \subset ({AC} \cup {BC}))
  ## method = 1({EFG} \subset ({ACD} \cup {BCD}))
  ## method = 2({EF} \subset ({ACD} \cup {BCD}))
  if(method == 0){
    method_txt <- "method = 0: Extra_CIS(for<A,B|C>) = {<E,F|G> | {EFG} subset ({AC} cup {BC})}\n"
  } else if(method == 1){
    method_txt <- "method = 1: Extra_CIS(for<A,B|C>) = {<E,F|G> | {EFG} subset ({ACD} cup {BCD})}\n"
  } else if(method == 2){
    method_txt <- "method = 2: Extra_CIS(for<A,B|C>) = {<E,F|G> | {EF} subset ({ACD} cup {BCD})}\n"
  }
  ####################
  n1 <- n-1
  n2 <- (2**n)-1
  elementary_imsets_array <- gen_elementary_imsets(n1, n2)
  elementary_imsets_array_colnames <- colnames(elementary_imsets_array)
  elementary_imsets_array_rownames <- rownames(elementary_imsets_array)
  nrow.elementary_imsets_array <- nrow(elementary_imsets_array)
  ####################
  # known_CIS_abc_vec <- replace_num_to_abc(known_CIS_vec)
  known_CIS_vec <- replace_abc_to_num(known_CIS_abc_vec)
  known_CIS_vec_org <- known_CIS_vec
  known_CIS_vec_names <- known_CIS_vec
  known_CIS_semi_elementary_imsets_array <- convert_CIS_vec_into_semi_elementary_imsets_array(known_CIS_vec, n2, elementary_imsets_array_rownames)
  #
  unknown_CIS_txt <- replace_abc_to_num(unknown_CIS_abc_txt)
  unknown_CIS_semi_elementary_imset_array <- convert_CIS_vec_into_semi_elementary_imsets_array(unknown_CIS_txt, n2, elementary_imsets_array_rownames)
  unknown_CIS_semi_elementary_imset_array_one_col <- as.matrix(unknown_CIS_semi_elementary_imset_array[,1])
  colnames(unknown_CIS_semi_elementary_imset_array_one_col) <- colnames(unknown_CIS_semi_elementary_imset_array)[1]
  #
  elementary_CIS_vec <- colnames(elementary_imsets_array)
  elementary_CIS_abc_vec <- replace_num_to_abc(elementary_CIS_vec)
  catalyst_CIS_vec <- elementary_CIS_vec[is.na(match(elementary_CIS_vec, known_CIS_vec))]
  catalyst_CIS_list <- list()
  res.test <- list()
  res.unknown.test <- res.test
  solved <- FALSE
  for(j in (1:nrow.elementary_imsets_array)){
    res.BHLS_test <- test_by_BHLS(n2, known_CIS_semi_elementary_imsets_array, unknown_CIS_semi_elementary_imset_array_one_col, elementary_imsets_array, known_CIS_vec=known_CIS_vec_names)
    res.test <- c(res.BHLS_test, method="BHLS", iter=j)
    res.unknown.test <- res.test
    if(res.test$solved==TRUE){
      solved <- TRUE
      break
    }
    res.TSTS_extra_test <- test_by_TSTS_extra(n2, known_CIS_semi_elementary_imsets_array, unknown_CIS_semi_elementary_imset_array_one_col, elementary_imsets_array, method = method, known_CIS_vec=known_CIS_vec_names)
    res.test <- c(res.TSTS_extra_test, method="TSTS_extra", iter=j)
    res.unknown.test <- res.test
    if(res.test$solved==TRUE){
      solved <- TRUE
      break
    }
    #
    len.catalyst_CIS_vec <- length(catalyst_CIS_vec)
    solved.catalyst_CIS_list <- list(CIS=c(), names=c())
    catalyst_CIS_semi_elementary_imsets_array <- convert_CIS_vec_into_semi_elementary_imsets_array(catalyst_CIS_vec, n2, elementary_imsets_array_rownames)
    changed_flag <- FALSE
    for(k in (1:len.catalyst_CIS_vec)){
      catalyst_CIS_semi_elementary_imsets_array_one_col <- as.matrix(catalyst_CIS_semi_elementary_imsets_array[,k])
      colnames(catalyst_CIS_semi_elementary_imsets_array_one_col) <- colnames(catalyst_CIS_semi_elementary_imsets_array)[k]
      res.BHLS_test <- test_by_BHLS(n2, known_CIS_semi_elementary_imsets_array, catalyst_CIS_semi_elementary_imsets_array_one_col, elementary_imsets_array)
      res.test <- c(res.BHLS_test, method="BHLS", iter=j)
      if(res.test$solved==TRUE){
        changed_flag <- TRUE
        catalyst_CIS_list <- c( catalyst_CIS_list, list(CIS=c(replace_num_to_abc(catalyst_CIS_vec[k]), method="BHLS", iter=j, res.test)) )
        solved.catalyst_CIS_list$CIS <- c(solved.catalyst_CIS_list$CIS, catalyst_CIS_vec[k])
        solved.catalyst_CIS_list$names <- c(solved.catalyst_CIS_list$names, paste(catalyst_CIS_vec[k], "[BHLS ", j, "th]", sep="", collapse=""))
        if(prncat==TRUE){
          print("BHLS")
          print(j)
          print(replace_num_to_abc(catalyst_CIS_vec[k]))
        }
        next
      }
      res.TSTS_extra_test <- test_by_TSTS_extra(n2, known_CIS_semi_elementary_imsets_array, catalyst_CIS_semi_elementary_imsets_array_one_col, elementary_imsets_array, method = method)
      res.test <- c(res.TSTS_extra_test, method="TSTS_extra", iter=j)
      if(res.test$solved==TRUE){
        changed_flag <- TRUE
        catalyst_CIS_list <- c( catalyst_CIS_list, list(CIS=c(replace_num_to_abc(catalyst_CIS_vec[k]), method="TSTS_extra", iter=j, res.test)) )
        solved.catalyst_CIS_list$CIS <- c(solved.catalyst_CIS_list$CIS, catalyst_CIS_vec[k])
        solved.catalyst_CIS_list$names <- c(solved.catalyst_CIS_list$names, paste(catalyst_CIS_vec[k], "[TSTS_extra ", j, "th]", sep="", collapse=""))
        if(prncat==TRUE){
          print("TSTS_extra")
          print(j)
          print(replace_num_to_abc(catalyst_CIS_vec[k]))
        }
        next
      }
    }
    ## renew known_CIS_semi_elementary_imsets_array here -->
    if(!changed_flag){
      break
    }
    known_CIS_vec <- c(known_CIS_vec, solved.catalyst_CIS_list$CIS)
    known_CIS_vec_names <- c(known_CIS_vec_names, solved.catalyst_CIS_list$names)
    known_CIS_semi_elementary_imsets_array <- convert_CIS_vec_into_semi_elementary_imsets_array(known_CIS_vec, n2, elementary_imsets_array_rownames)
    catalyst_CIS_vec <- elementary_CIS_vec[is.na(match(elementary_CIS_vec, known_CIS_vec))]
  }
  #
  return(list(main=res.unknown.test, catalyst=catalyst_CIS_list))
}
CIS.test.restriction <- function(n, known, unknown, method, prncat=FALSE){
  n <- n
  known_CIS_abc_vec <- known
  unknown_CIS_abc_txt <- unknown
  method <- method
  ####################
  ## method = 0( ABCD \setminus (AC \cup BC) )
  ## method = 1( ABCD \setminus (ACD \cup BCD) )
  if(method == 0){
    method_txt <- "method = 0: K = ABCD setminus (AC cup BC)\n"
  } else if(method == 1){
    method_txt <- "method = 1: K = ABCD setminus (ACD cup BCD)\n"
  }
  ####################
  n1 <- n-1
  n2 <- (2**n)-1
  elementary_imsets_array <- gen_elementary_imsets(n1, n2)
  elementary_imsets_array_colnames <- colnames(elementary_imsets_array)
  elementary_imsets_array_rownames <- rownames(elementary_imsets_array)
  nrow.elementary_imsets_array <- nrow(elementary_imsets_array)
  ####################
  # known_CIS_abc_vec <- replace_num_to_abc(known_CIS_vec)
  known_CIS_vec <- replace_abc_to_num(known_CIS_abc_vec)
  known_CIS_vec_org <- known_CIS_vec
  known_CIS_vec_names <- known_CIS_vec
  known_CIS_semi_elementary_imsets_array <- convert_CIS_vec_into_semi_elementary_imsets_array(known_CIS_vec, n2, elementary_imsets_array_rownames)
  #
  unknown_CIS_txt <- replace_abc_to_num(unknown_CIS_abc_txt)
  unknown_CIS_semi_elementary_imset_array <- convert_CIS_vec_into_semi_elementary_imsets_array(unknown_CIS_txt, n2, elementary_imsets_array_rownames)
  unknown_CIS_semi_elementary_imset_array_one_col <- as.matrix(unknown_CIS_semi_elementary_imset_array[,1])
  colnames(unknown_CIS_semi_elementary_imset_array_one_col) <- colnames(unknown_CIS_semi_elementary_imset_array)[1]
  #
  elementary_CIS_vec <- colnames(elementary_imsets_array)
  elementary_CIS_abc_vec <- replace_num_to_abc(elementary_CIS_vec)
  catalyst_CIS_vec <- elementary_CIS_vec[is.na(match(elementary_CIS_vec, known_CIS_vec))]
  catalyst_CIS_list <- list()
  res.test <- list()
  res.unknown.test <- res.test
  solved <- FALSE
  for(j in (1:nrow.elementary_imsets_array)){
    res.BHLS_test <- test_by_BHLS(n2, known_CIS_semi_elementary_imsets_array, unknown_CIS_semi_elementary_imset_array_one_col, elementary_imsets_array, known_CIS_vec=known_CIS_vec_names)
    res.test <- c(res.BHLS_test, method="BHLS", iter=j)
    res.unknown.test <- res.test
    if(res.test$solved==TRUE){
      solved <- TRUE
      break
    }
    res.TSTS_extra_test <- test_by_TSTS_extra(n2, known_CIS_semi_elementary_imsets_array, unknown_CIS_semi_elementary_imset_array_one_col, elementary_imsets_array, method = method, known_CIS_vec=known_CIS_vec_names)
    res.test <- c(res.TSTS_extra_test, method="TSTS_extra", iter=j)
    res.unknown.test <- res.test
    if(res.test$solved==TRUE){
      solved <- TRUE
      break
    }
    #
    len.catalyst_CIS_vec <- length(catalyst_CIS_vec)
    solved.catalyst_CIS_list <- list(CIS=c(), names=c())
    catalyst_CIS_semi_elementary_imsets_array <- convert_CIS_vec_into_semi_elementary_imsets_array(catalyst_CIS_vec, n2, elementary_imsets_array_rownames)
    changed_flag <- FALSE
    for(k in (1:len.catalyst_CIS_vec)){
      catalyst_CIS_semi_elementary_imsets_array_one_col <- as.matrix(catalyst_CIS_semi_elementary_imsets_array[,k])
      colnames(catalyst_CIS_semi_elementary_imsets_array_one_col) <- colnames(catalyst_CIS_semi_elementary_imsets_array)[k]
      res.BHLS_test <- test_by_BHLS(n2, known_CIS_semi_elementary_imsets_array, catalyst_CIS_semi_elementary_imsets_array_one_col, elementary_imsets_array)
      res.test <- c(res.BHLS_test, method="BHLS", iter=j)
      if(res.test$solved==TRUE){
        changed_flag <- TRUE
        catalyst_CIS_list <- c( catalyst_CIS_list, list(CIS=c(replace_num_to_abc(catalyst_CIS_vec[k]), method="BHLS", iter=j, res.test)) )
        solved.catalyst_CIS_list$CIS <- c(solved.catalyst_CIS_list$CIS, catalyst_CIS_vec[k])
        solved.catalyst_CIS_list$names <- c(solved.catalyst_CIS_list$names, paste(catalyst_CIS_vec[k], "[BHLS ", j, "th]", sep="", collapse=""))
        if(prncat==TRUE){
          print("BHLS")
          print(j)
          print(replace_num_to_abc(catalyst_CIS_vec[k]))
        }
        next
      }
      res.TSTS_extra_test <- test_by_TSTS_extra(n2, known_CIS_semi_elementary_imsets_array, catalyst_CIS_semi_elementary_imsets_array_one_col, elementary_imsets_array, method = method)
      res.test <- c(res.TSTS_extra_test, method="TSTS_extra", iter=j)
      if(res.test$solved==TRUE){
        changed_flag <- TRUE
        catalyst_CIS_list <- c( catalyst_CIS_list, list(CIS=c(replace_num_to_abc(catalyst_CIS_vec[k]), method="TSTS_extra", iter=j, res.test)) )
        solved.catalyst_CIS_list$CIS <- c(solved.catalyst_CIS_list$CIS, catalyst_CIS_vec[k])
        solved.catalyst_CIS_list$names <- c(solved.catalyst_CIS_list$names, paste(catalyst_CIS_vec[k], "[TSTS_extra ", j, "th]", sep="", collapse=""))
        if(prncat==TRUE){
          print("TSTS_extra")
          print(j)
          print(replace_num_to_abc(catalyst_CIS_vec[k]))
        }
        next
      }
    }
    ## renew known_CIS_semi_elementary_imsets_array here -->
    if(!changed_flag){
      break
    }
    known_CIS_vec <- c(known_CIS_vec, solved.catalyst_CIS_list$CIS)
    known_CIS_vec_names <- c(known_CIS_vec_names, solved.catalyst_CIS_list$names)
    known_CIS_semi_elementary_imsets_array <- convert_CIS_vec_into_semi_elementary_imsets_array(known_CIS_vec, n2, elementary_imsets_array_rownames)
    catalyst_CIS_vec <- elementary_CIS_vec[is.na(match(elementary_CIS_vec, known_CIS_vec))]
  }
  #
  return(list(main=res.unknown.test, catalyst=catalyst_CIS_list))
}
CIS_for4ti2 <- function(n, known){
  n <- n
  known_CIS_abc_vec <- known
  known_CIS_vec <- replace_abc_to_num(known_CIS_abc_vec)
  elementary_imsets_array <- gen_elementary_imsets(n1, n2)
  elementary_imsets_array_colnames <- colnames(elementary_imsets_array)
  for4ti2_mat <- elementary_imsets_array[, !is.na(match(elementary_imsets_array_colnames, known_CIS_vec))]
  systime_txt <- Sys.time()
  systime_txt <- gsub(":", "", systime_txt)
  systime_txt <- gsub("-", "", systime_txt)
  systime_txt <- gsub("[[:space:]]", "", systime_txt)
  filename <- paste("./for4ti2_", n, "var", dim(for4ti2_mat)[1], "row", dim(for4ti2_mat)[2], "col_", systime_txt, ".mat", sep="")
  out <- file(filename, "w")
  writeLines(paste(dim(for4ti2_mat)[1], " ", dim(for4ti2_mat)[2], "\n", sep="", collapse=""), out, sep="")
  write.table(for4ti2_mat, out, sep=" ", row.names = FALSE, col.names = FALSE)
  close(out)
}
##
CIS_test <- function(n, known_CIS_abc_vec_org, known_CIS_abc_vec, unknown_CIS_abc_vec, method = "BHLS", positive=TRUE, tol = 1e-8, printout=TRUE){
  # method = BHLS, TSTS
  # positive = TRUE, FALSE
  known_CIS_abc_vec_txt <- paste(known_CIS_abc_vec_org, collapse=", ", sep="")
  n1 <- n-1
  n2 <- (2**n)-1
  elementary_imsets_array <- gen_elementary_imsets(n1, n2)
  elementary_imsets_array_rownames <- rownames(elementary_imsets_array)
  known_CIS_vec <- replace_abc_to_num(known_CIS_abc_vec)
  known_CIS_semi_elementary_imsets_array <- convert_CIS_vec_into_semi_elementary_imsets_array(known_CIS_vec, n2, elementary_imsets_array_rownames)
  unknown_CIS_vec <- replace_abc_to_num(unknown_CIS_abc_vec)
  unknown_CIS_semi_elementary_imsets_array <- convert_CIS_vec_into_semi_elementary_imsets_array(unknown_CIS_vec, n2, elementary_imsets_array_rownames)
  len.unknown_CIS_vec <- length(unknown_CIS_abc_vec)
  res.solved.vec <- c()
  res.solved.positive.vec <- c()
  res.solved.zero.vec <- c()
  method.txt <- ""
  res.txt <- ""
  if(method=="BHLS"){
    method.txt <- "[BHLS]"
    for(i in (1:len.unknown_CIS_vec)){
      colnames_unknown_CIS_semi_elementary_imset_array_ith_col <- colnames(unknown_CIS_semi_elementary_imsets_array)[i]
      unknown_CIS_semi_elementary_imset_array <- as.matrix(unknown_CIS_semi_elementary_imsets_array[,i])
      colnames(unknown_CIS_semi_elementary_imset_array) <- colnames_unknown_CIS_semi_elementary_imset_array_ith_col
      res_test_by_BHLS <- test_by_BHLS(n2, known_CIS_semi_elementary_imsets_array, unknown_CIS_semi_elementary_imset_array, elementary_imsets_array, tol)
      res.solved.vec <- c(res.solved.vec, res_test_by_BHLS$solved)
      res.txt <- paste(res.solved.vec[i], method.txt, ": ", known_CIS_abc_vec_txt, " => ", unknown_CIS_abc_vec[i], collapse="", sep="")
      if(printout){
        print(res.txt)
      }
    }
  } else if(method=="TSTS"){
    if(positive==TRUE){
      method.txt <- "[TSTS(positive)]"
      for(i in (1:len.unknown_CIS_vec)){
        colnames_unknown_CIS_semi_elementary_imset_array_ith_col <- colnames(unknown_CIS_semi_elementary_imsets_array)[i]
        unknown_CIS_semi_elementary_imset_array <- as.matrix(unknown_CIS_semi_elementary_imsets_array[,i])
        colnames(unknown_CIS_semi_elementary_imset_array) <- colnames_unknown_CIS_semi_elementary_imset_array_ith_col
        res_test_by_TSTS <- test_by_TSTS_restriction(n2, known_CIS_semi_elementary_imsets_array, unknown_CIS_semi_elementary_imset_array, elementary_imsets_array, method=1, tol)
        res.solved.vec <- c(res.solved.vec, res_test_by_TSTS$solved)
        res.txt <- paste(res.solved.vec[i], method.txt, ": ", known_CIS_abc_vec_txt, " => ", unknown_CIS_abc_vec[i], collapse="", sep="")
        if(printout){
          print(res.txt)
        }
      }
    } else{
      method.txt <- "[TSTS(not necessarily positive)]"
      for(i in (1:len.unknown_CIS_vec)){
        colnames_unknown_CIS_semi_elementary_imset_array_ith_col <- colnames(unknown_CIS_semi_elementary_imsets_array)[i]
        unknown_CIS_semi_elementary_imset_array <- as.matrix(unknown_CIS_semi_elementary_imsets_array[,i])
        colnames(unknown_CIS_semi_elementary_imset_array) <- colnames_unknown_CIS_semi_elementary_imset_array_ith_col
        res_test_by_TSTS <- test_by_TSTS_restriction(n2, known_CIS_semi_elementary_imsets_array, unknown_CIS_semi_elementary_imset_array, elementary_imsets_array, method=0, tol)
        res.solved.positive.vec <- c(res.solved.positive.vec, res_test_by_TSTS$solved.positive)
        res.solved.zero.vec <- c(res.solved.zero.vec, res_test_by_TSTS$solved.zero)
        res.solved.vec <- c(res.solved.vec, res_test_by_TSTS$solved)
        res.txt <- paste(res.solved.vec[i], "(positive=", res.solved.positive.vec[i], ", zero=", res.solved.zero.vec[i], ")", method.txt, ": ", known_CIS_abc_vec_txt, " => ", unknown_CIS_abc_vec[i], collapse="", sep="")
        if(printout){
          print(res.txt)
        }
      }
    }
  } else if(method=="MIXED"){
    for(i in (1:len.unknown_CIS_vec)){
      colnames_unknown_CIS_semi_elementary_imset_array_ith_col <- colnames(unknown_CIS_semi_elementary_imsets_array)[i]
      unknown_CIS_semi_elementary_imset_array <- as.matrix(unknown_CIS_semi_elementary_imsets_array[,i])
      colnames(unknown_CIS_semi_elementary_imset_array) <- colnames_unknown_CIS_semi_elementary_imset_array_ith_col
      res_test_by_BHLS <- test_by_BHLS(n2, known_CIS_semi_elementary_imsets_array, unknown_CIS_semi_elementary_imset_array, elementary_imsets_array, tol)
      if(res_test_by_BHLS$solved){
        method.txt <- "[MIXED:BHLS]"
        res.solved.vec <- c(res.solved.vec, res_test_by_BHLS$solved)
        res.txt <- paste(res.solved.vec[i], method.txt, ": ", known_CIS_abc_vec_txt, " => ", unknown_CIS_abc_vec[i], collapse="", sep="")
      } else{
        if(positive==TRUE){
          method.txt <- "[MIXED:TSTS(positive)]"
          res_test_by_TSTS <- test_by_TSTS_restriction(n2, known_CIS_semi_elementary_imsets_array, unknown_CIS_semi_elementary_imset_array, elementary_imsets_array, method=1, tol)
          res.solved.vec <- c(res.solved.vec, res_test_by_TSTS$solved)
          if(res_test_by_TSTS$solved){
            res.txt <- paste(res.solved.vec[i], method.txt, ": ", known_CIS_abc_vec_txt, " => ", unknown_CIS_abc_vec[i], collapse="", sep="")
          } else{
            res.txt <- paste(res.solved.vec[i], "[MIXED:BHLS,TSTS(positive)]", ": ", known_CIS_abc_vec_txt, " => ", unknown_CIS_abc_vec[i], collapse="", sep="")
          }
        } else{
          method.txt <- "[MIXED:TSTS(not necessarily positive)]"
          res_test_by_TSTS <- test_by_TSTS_restriction(n2, known_CIS_semi_elementary_imsets_array, unknown_CIS_semi_elementary_imset_array, elementary_imsets_array, method=0, tol)
          res.solved.positive.vec <- c(res.solved.positive.vec, res_test_by_TSTS$solved.positive)
          res.solved.zero.vec <- c(res.solved.zero.vec, res_test_by_TSTS$solved.zero)
          res.solved.vec <- c(res.solved.vec, res_test_by_TSTS$solved)
          if(res_test_by_TSTS$solved){
            res.txt <- paste(res.solved.vec[i], "(positive=", res.solved.positive.vec[i], ", zero=", res.solved.zero.vec[i], ")", method.txt, ": ", known_CIS_abc_vec_txt, " => ", unknown_CIS_abc_vec[i], collapse="", sep="")
          } else{
            res.txt <- paste(res.solved.vec[i], "[MIXED:BHLS,TSTS(not necessarily positive)]", ": ", known_CIS_abc_vec_txt, " => ", unknown_CIS_abc_vec[i], collapse="", sep="")
          }
        }
      }
      if(printout){
        print(res.txt)
      }
    }
  }
  else{
    print("ERROR!! -- INVALID METHOD (TK)--")
  }
  return(list(solved=res.solved.vec, solved.positive=res.solved.positive.vec, solved.zero=res.solved.zero.vec))
}
CIS_closure <- function(n, known_CIS_abc_vec, method = "BHLS", positive=TRUE, tol = 1e-8, loop_max=100, printout=TRUE){
  # method = BHLS, TSTS, MIXED
  # positive = TRUE, FALSE
  known_CIS_abc_vec_org <- known_CIS_abc_vec
  known_CIS_abc_vec_txt <- paste(known_CIS_abc_vec, collapse=", ", sep="")
  n1 <- n-1
  n2 <- (2**n)-1
  elementary_imsets_array <- gen_elementary_imsets(n1, n2)
  elementary_imsets_array_rownames <- rownames(elementary_imsets_array)
  # ncol.elementary_imsets_array <- ncol(elementary_imsets_array)
  known_CIS_vec <- replace_abc_to_num(known_CIS_abc_vec)
  known_CIS_semi_elementary_imsets_array <- convert_CIS_vec_into_semi_elementary_imsets_array(known_CIS_vec, n2, elementary_imsets_array_rownames)
  method.txt <- paste("[", method, "]", collapse="", sep="")
  mixed.flag <- 0

  unsolved_elementary_imsets_array <- elementary_imsets_array
  # closure
  for(j in (1:loop_max)){
    ncol.unsolved_elementary_imsets_array <- ncol(unsolved_elementary_imsets_array)
    solved_vec <- c()
    for(i in (1:ncol.unsolved_elementary_imsets_array)){
      colnames_unknown_CIS_semi_elementary_imset_array_ith_col <- colnames(unsolved_elementary_imsets_array)[i]
      unknown_CIS_semi_elementary_imset_array <- as.matrix(unsolved_elementary_imsets_array[,i])
      colnames(unknown_CIS_semi_elementary_imset_array) <- colnames_unknown_CIS_semi_elementary_imset_array_ith_col
      if(method=="BHLS"){
        res <- test_by_BHLS(n2, known_CIS_semi_elementary_imsets_array, unknown_CIS_semi_elementary_imset_array, elementary_imsets_array, tol = 1e-8)
      } else if(method=="TSTS"){
        if(positive==TRUE){
          res <- test_by_TSTS_restriction(n2, known_CIS_semi_elementary_imsets_array, unknown_CIS_semi_elementary_imset_array, elementary_imsets_array, method=1, tol = 1e-8)
        } else{
          res <- test_by_TSTS_restriction(n2, known_CIS_semi_elementary_imsets_array, unknown_CIS_semi_elementary_imset_array, elementary_imsets_array, method=0, tol = 1e-8)
        }
      } else if(method=="MIXED"){
        if(j%%2==0){
          res <- test_by_BHLS(n2, known_CIS_semi_elementary_imsets_array, unknown_CIS_semi_elementary_imset_array, elementary_imsets_array, tol = 1e-8)
        } else{
          if(positive==TRUE){
            res <- test_by_TSTS_restriction(n2, known_CIS_semi_elementary_imsets_array, unknown_CIS_semi_elementary_imset_array, elementary_imsets_array, method=1, tol = 1e-8)
          } else{
            res <- test_by_TSTS_restriction(n2, known_CIS_semi_elementary_imsets_array, unknown_CIS_semi_elementary_imset_array, elementary_imsets_array, method=0, tol = 1e-8)
          }
        }
      } else{
        print("ERROR!! -- INVALID METHOD [CIS_closed_test] (TK)--")
      }
      if(res$solved){
        if(!(any(colnames(known_CIS_semi_elementary_imsets_array)==colnames(unknown_CIS_semi_elementary_imset_array)))){
          known_CIS_semi_elementary_imsets_array <- cbind(known_CIS_semi_elementary_imsets_array, unknown_CIS_semi_elementary_imset_array)
        }
        solved_vec <- c(solved_vec, i)
      }
    }
    if(method=="MIXED"){
      if(length(solved_vec)==0){
        mixed.flag <- mixed.flag + 1
      } else{
        mixed.flag <- 0
        colnames_unsolved_elementary_imsets_array <- colnames(unsolved_elementary_imsets_array)[-solved_vec]
        unsolved_elementary_imsets_array <- as.matrix(unsolved_elementary_imsets_array[,-solved_vec])
        colnames(unsolved_elementary_imsets_array) <- colnames_unsolved_elementary_imsets_array
      }
      if(mixed.flag>=2){
        break
      }
    } else{
      if(length(solved_vec)==0){
        break
      } else{
        colnames_unsolved_elementary_imsets_array <- colnames(unsolved_elementary_imsets_array)[-solved_vec]
        unsolved_elementary_imsets_array <- as.matrix(unsolved_elementary_imsets_array[,-solved_vec])
        colnames(unsolved_elementary_imsets_array) <- colnames_unsolved_elementary_imsets_array
      }
    }
    if(ncol(unsolved_elementary_imsets_array)==0){
      break
    }
    if(j==loop_max){
      print("WARNING!! -- NOT CONVERGED (TK)--")
    }
  }
  # test
  closure_of_known_CIS_abc_vec <- replace_num_to_abc(colnames(known_CIS_semi_elementary_imsets_array))
  return(closure_of_known_CIS_abc_vec)
}
CIS_closed_test <- function(n, known_CIS_abc_vec, unknown_CIS_abc_vec, method = "BHLS", positive=TRUE, tol = 1e-8, loop_max=100, printout=TRUE){
  closure_of_known_CIS_abc_vec <- CIS_closure(n, known_CIS_abc_vec, method, positive, tol, loop_max, printout)
  res <- CIS_test(n, known_CIS_abc_vec, closure_of_known_CIS_abc_vec, unknown_CIS_abc_vec, method=method, positive=positive, tol=tol, printout=TRUE)
  return(c(res, list(closure=c(closure_of_known_CIS_abc_vec))))
}
CIS_five_types_of_closed_tests <- function(n, known_CIS_abc_vec, unknown_CIS_abc_vec, tol = 1e-8, loop_max=100, printout=TRUE){
  res.CIS.BHLS <- CIS_closed_test(n, known_CIS_abc_vec, unknown_CIS_abc_vec, method = "BHLS", positive=TRUE, tol=tol, loop_max=loop_max, printout=printout)
  res.CIS.TSTS.positive <- CIS_closed_test(n, known_CIS_abc_vec, unknown_CIS_abc_vec, method = "TSTS", positive=TRUE, tol=tol, loop_max=loop_max, printout=printout)
  res.CIS.TSTS.zero <- CIS_closed_test(n, known_CIS_abc_vec, unknown_CIS_abc_vec, method = "TSTS", positive=FALSE, tol=tol, loop_max=loop_max, printout=printout)
  res.CIS.MIXED.positive <- CIS_closed_test(n, known_CIS_abc_vec, unknown_CIS_abc_vec, method = "MIXED", positive=TRUE, tol=tol, loop_max=loop_max, printout=printout)
  res.CIS.MIXED.zero <- CIS_closed_test(n, known_CIS_abc_vec, unknown_CIS_abc_vec, method = "MIXED", positive=FALSE, tol=tol, loop_max=loop_max, printout=printout)
  return(c(res.BHLS=list(res.CIS.BHLS), res.TSTS.positive=list(res.CIS.TSTS.positive), res.TSTS.zero=list(res.CIS.TSTS.zero), res.MIXED.positive=list(res.CIS.MIXED.positive), res.MIXED.zero=list(res.CIS.MIXED.zero)))
}
CIS_five_types_of_closures <- function(n, known_CIS_abc_vec, tol = 1e-8, loop_max=100, printout=TRUE){
  res.CIS.BHLS <- CIS_closure(n, known_CIS_abc_vec, method = "BHLS", positive=TRUE, tol=tol, loop_max=loop_max, printout=printout)
  res.CIS.TSTS.positive <- CIS_closure(n, known_CIS_abc_vec, method = "TSTS", positive=TRUE, tol=tol, loop_max=loop_max, printout=printout)
  res.CIS.TSTS.zero <- CIS_closure(n, known_CIS_abc_vec, method = "TSTS", positive=FALSE, tol=tol, loop_max=loop_max, printout=printout)
  res.CIS.MIXED.positive <- CIS_closure(n, known_CIS_abc_vec, method = "MIXED", positive=TRUE, tol=tol, loop_max=loop_max, printout=printout)
  res.CIS.MIXED.zero <- CIS_closure(n, known_CIS_abc_vec, method = "MIXED", positive=FALSE, tol=tol, loop_max=loop_max, printout=printout)
  return(c(org=list(known_CIS_abc_vec), closure.BHLS=list(res.CIS.BHLS), closure.TSTS.positive=list(res.CIS.TSTS.positive), closure.TSTS.zero=list(res.CIS.TSTS.zero), closure.MIXED.positive=list(res.CIS.MIXED.positive), closure.MIXED.zero=list(res.CIS.MIXED.zero)))
}
#### subfunc
####################
