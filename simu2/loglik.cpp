//#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <Rmath.h>
#include <iterator>
#include <iostream>
#include <vector>
#include <algorithm>

#include "commonf.h"

//using namespace Rcpp;
using namespace arma;


// Obtain environment containing function
//Rcpp::Environment package_env("package:mnormt"); 

// Make function callable from C++
//Rcpp::Function bipmvnrm = package_env["pmnorm"];    

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export()]]
NumericVector Lcomp_adjusted_1j(NumericVector l, NumericVector r, NumericVector cov, 
                                NumericVector alp, NumericVector brks){
  
  int n = l.size();
  NumericVector res(n);
  
  for(int i=0; i<n; i++){
    
    NumericVector lev = alp*exp(cov[i]);
    res[i] = ppwc_double(r[i], brks, lev, 1, 0) - ppwc_double(l[i], brks, lev, 1, 0);
  }
  
  return res;
}

// [[Rcpp::export()]]
List get_indata_s2_jjp(NumericVector phi1_j, NumericVector phi1_jp,
                       NumericVector lt_j, NumericVector rt_j,
                       NumericVector lt_jp, NumericVector rt_jp,
                       NumericMatrix Xmat_j, NumericMatrix Xmat_jp, 
                       NumericVector cuts_j, NumericVector cuts_jp){
  
  int n = lt_j.size();
  
  int XSize_j = Xmat_j.ncol();
  int XSize_jp = Xmat_jp.ncol();
  int cutSize_j = cuts_j.size();
  int cutSize_jp = cuts_jp.size();
  
  NumericVector levels_j = exp(phi1_j[Range(0,cutSize_j)]);
  NumericVector levels_jp = exp(phi1_jp[Range(0,cutSize_jp)]);
  
  NumericVector eta_j = phi1_j[Range(cutSize_j+1, cutSize_j+XSize_j+1)];
  NumericVector eta_jp = phi1_jp[Range(cutSize_jp+1, cutSize_jp+XSize_jp+1)];
  
  NumericVector pi_marg_j(n), survlt_j(n), survrt_j(n), invPhisurvlt_j(n), invPhisurvrt_j(n);
  NumericVector pi_marg_jp(n), survlt_jp(n), survrt_jp(n), invPhisurvlt_jp(n), invPhisurvrt_jp(n);
  
  for(int i=0; i<n; i++){
    
    survlt_j[i] = ppwc_double(lt_j[i], cuts_j, levels_j, 0, 0);
    survrt_j[i] = ppwc_double(rt_j[i], cuts_j, levels_j, 0, 0);
    survlt_jp[i] = ppwc_double(lt_jp[i], cuts_jp, levels_jp, 0, 0);
    survrt_jp[i] = ppwc_double(rt_jp[i], cuts_jp, levels_jp, 0, 0);
    invPhisurvlt_j[i] = R::qnorm(survlt_j[i], 0.0, 1.0, 1, 0);
    invPhisurvrt_j[i] = R::qnorm(survrt_j[i], 0.0, 1.0, 1, 0);
    invPhisurvlt_jp[i] = R::qnorm(survlt_jp[i], 0.0, 1.0, 1, 0);
    invPhisurvrt_jp[i] = R::qnorm(survrt_jp[i], 0.0, 1.0, 1, 0);
    
    double tmpcov_j = eta_j[0] + sum(eta_j[Range(1,XSize_j)]*Xmat_j(i,_));
    double tmpcov_jp = eta_jp[0] + sum(eta_jp[Range(1,XSize_jp)]*Xmat_jp(i,_));
    pi_marg_j[i] = exp(tmpcov_j)*pow((1+exp(tmpcov_j)), -1);
    pi_marg_jp[i] = exp(tmpcov_jp)*pow((1+exp(tmpcov_jp)), -1);
    
  }
  
  return List::create(Named("pi_marg_j")=pi_marg_j,
                      Named("pi_marg_jp")=pi_marg_jp,
                      Named("survlt_j")=survlt_j,
                      Named("survrt_j")=survrt_j,
                      Named("survlt_jp")=survlt_jp,
                      Named("survrt_jp")=survrt_jp,
                      Named("invPhisurvlt_j")=invPhisurvlt_j,
                      Named("invPhisurvrt_j")=invPhisurvrt_j,
                      Named("invPhisurvlt_jp")=invPhisurvlt_jp,
                      Named("invPhisurvrt_jp")=invPhisurvrt_jp);
}


// [[Rcpp::export]]
NumericVector newevalLogCLPair_jjp(double gamma, NumericVector d2_jjp, 
                                   NumericVector xi_j, NumericVector xi_jp,
                                   NumericVector pi_marg_j, NumericVector pi_marg_jp,
                                   NumericVector survlt_j, NumericVector survrt_j,
                                   NumericVector survlt_jp, NumericVector survrt_jp){
  
  int n = pi_marg_j.size();
  
  //double gamma = phi2_jjp(0); // , rho = atan(phi2_jjp(1))*2/M_PI;
  
  //NumericMatrix sigma(2,2);
  
  //Rcout << "sigma is " << sigma << std::endl;
  
  //sigma(0,1) = rho;
  //sigma(1,0) = rho;
  
  //Rcout << "sigma is " << sigma << std::endl;
  
  //sigma.fill_diag(1.0);
  
  NumericVector loglikPair(n);
  
  for(int i=0; i<n; i++){
    
    double tmpCL = 0.0;
    
    //NumericMatrix mat1(1, 2), mat2(1, 2), mat3(1, 2), mat4(1, 2);
    
    //Rcout << " mat1 is " << mat1 << std::endl;
    
    //mat1(0,0) = invPhisurvlt_j[i];
    //mat1(0,1) = invPhisurvlt_jp[i];
    //mat2(0,0) = invPhisurvlt_j[i];
    //mat2(0,1) = invPhisurvrt_jp[i];
    //mat3(0,0) = invPhisurvrt_j[i];
    //mat3(0,1) = invPhisurvlt_jp[i];
    //mat4(0,0) = invPhisurvrt_j[i];
    //mat4(0,1) = invPhisurvrt_jp[i];
    
    //Rcout << " mat1 is " << mat1 << std::endl;
    
    //Rcout << " cholmat1 is " << arma::chol(as<arma::mat>(sigma)) << std::endl;
    
    
    //NumericVector res1, res2, res3, res4;
    //arma::rowvec mean(2);
    //res1 = wrap(dmvnrm(as<arma::mat>(mat1), mean.fill(0.0), as<arma::mat>(sigma), 0));
    //res2 = as<NumericVector>(wrap(dmvnrm(as<arma::mat>(mat2), mean.fill(0.0), as<arma::mat>(sigma), 0)));
    //res3 = as<NumericVector>(wrap(dmvnrm(as<arma::mat>(mat3), mean.fill(0.0), as<arma::mat>(sigma), 0)));
    //res4 = as<NumericVector>(wrap(dmvnrm(as<arma::mat>(mat4), mean.fill(0.0), as<arma::mat>(sigma), 0)));
    
    // Obtain environment containing function
    //Rcpp::Environment package_env("package:mnormt"); 
    
    // Make function callable from C++
    //Rcpp::Function bipmvnrm = package_env["pmnorm"];       
    
    //res1 = bipmvnrm(mat1, mean.fill(0.0), sigma);
    //res2 = bipmvnrm(mat2, mean.fill(0.0), sigma);
    //res3 = bipmvnrm(mat3, mean.fill(0.0), sigma);
    //res4 = bipmvnrm(mat4, mean.fill(0.0), sigma);
    
    //res1 = ifelse(is_na(res1), 0.0, res1);
    //res2 = ifelse(is_na(res1), 0.0, res2);
    //res3 = ifelse(is_na(res3), 0.0, res3);
    //res4 = ifelse(is_na(res4), 0.0, res4);
    
    NumericMatrix w_jkjpkp(2, 2);
    w_jkjpkp(1,1) = pi11_jkjpkp(pi_marg_j[i], pi_marg_jp[i], gamma);
    w_jkjpkp(0,1) =  pi_marg_jp[i] - w_jkjpkp(1,1); 
    w_jkjpkp(1,0) = pi_marg_j[i] - w_jkjpkp(1,1); 
    w_jkjpkp(0,0) = 1 - pi_marg_j[i] - pi_marg_jp[i] + w_jkjpkp(1,1); 
    
    
    if(xi_j[i]){
      
      if(xi_jp[i]){
        
        tmpCL += sum(log(w_jkjpkp(1,1)) + log(d2_jjp[i]));
      }else{
        
        tmpCL += sum(log( w_jkjpkp(1,1)*d2_jjp[i] 
                        + w_jkjpkp(1,0)*(survlt_j[i] - survrt_j[i])));
      }
      
    }else{
      
      if(xi_jp[i]){
        
        tmpCL += log(w_jkjpkp(1,1)*d2_jjp[i] 
                       + w_jkjpkp(0,1)*(survlt_jp[i] - survrt_jp[i]));
      }else{
        
        tmpCL += log( w_jkjpkp(1,1)*d2_jjp[i]
                        + w_jkjpkp(1,0)*(survlt_j[i] - survrt_j[i])
                        + w_jkjpkp(0,1)*(survlt_jp[i] - survrt_jp[i])
                        + 1 - pi_marg_j[i] - pi_marg_jp[i] + w_jkjpkp(1,1) );
      }
    }
    
    
    loglikPair[i] = tmpCL;
    
    
    //Rcout << " res1 is " << res1 << std::endl;
    //Rcout << " res2 is " << res2 << std::endl;
    //Rcout << " res3 is " << res3 << std::endl;
    //Rcout << " res4 is " << res4 << std::endl;
    
    //Rcout << " eta_jk is " << eta_jk << std::endl;
    //Rcout << " eta_jpkp is " << eta_jpkp << std::endl;
    
  }
  
  return loglikPair;
}





