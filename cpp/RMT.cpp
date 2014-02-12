#ifndef RMT_VARIOUS
#define RMT_VARIOUS


#include <itpp/itbase.h>
#include <itpp/stat/misc_stat.h>
#include <cpp/itpp_ext_math.cpp> 
#include "RMT.h"
namespace RMT{ // {{{ Headers
  itpp::mat RandomGOEDeltaOne(int const );
  itpp::Mat<std::complex<double> > RandomCUE(int const );
  void FlatSpectrumGSE(itpp::Mat<std::complex<double> >&, itpp::Vec<double>& );
} // }}}
namespace RMT{ // {{{ Implementation
// unfolding things {{{
  double unfolding_function(double const e){ //{{{
    return (asin(e)+e*sqrt(1-pow(e,2)))/itpp::pi;

  } // }}}
  bool unfoldcircular(itpp::Vec<double>& eigenvalues){ // {{{
    // aca se supone que entran de -1 a 1
    double e;
    bool is_everything_ok=true;
    for (int i=0; i< eigenvalues.size(); i++){
      e=eigenvalues(i);
      if(std::abs(e)<1){
        e=eigenvalues.size()*unfolding_function(e);
      }
      else{
        e=eigenvalues.size()*0.5*itpp::sign(e);
#ifndef NOWARNDRMT
        // std::cerr  << "Warning, unfoldcircular got wierd value "<< "Ignore if comming from safe routine "<<e <<std::endl ;
#endif /* NOWARNRMT */
        is_everything_ok=false;
      }
      eigenvalues(i)=e;
    }
    return is_everything_ok;
  }
  // }}}
// }}}
// Others (e.g. antisymmetric, isometries) {{{
  itpp::mat RandomAntisymmetricDeltaOne(int const dim){ //{{{
    itpp::mat temp(dim, dim);
    temp=itpp::randn(dim,dim);
    return sqrt(1/2.)*(temp-temp.hermitian_transpose());
  }
  // }}}
  itpp::Mat<std::complex<double> > RandomIsometryToSingleQubit(const int dim){ //{{{
    itpp::Mat<std::complex<double> > tmp=itpp::randn_c(1, dim);
    double norm=itpp::norm(tmp,2);
    return tmp/=norm ;
  };
  // }}}
// }}}
// GSE Collection {{{
  itpp::cmat RandomGSEDeltaOne(int const dim){ //{{{
    // Los eigenvalores de esta matriz estan entre +- Sqrt[8 * nsize]
    if (dim %2 != 0){
      std::cerr << "RandomGSEDeltaOne::dimension es impar, eso es invalido dim="<<dim
        << std::endl;
      abort();


    }
    itpp::mat Ha, Hb ; int d2=dim/2;
    itpp::cmat H(dim, dim);
    std::complex<double> I(0.,1.);
    Ha=RandomGOEDeltaOne(d2); Hb=RandomAntisymmetricDeltaOne(d2);
    H.set_submatrix(0,0,Ha-I*Hb);  H.set_submatrix(d2,d2,Ha+I*Hb);
    Ha=RandomAntisymmetricDeltaOne(d2); Hb=RandomAntisymmetricDeltaOne(d2);
    H.set_submatrix(0,d2,-Ha-I*Hb); H.set_submatrix(d2,0,Ha-I*Hb);
    return H;
  } // }}}
  itpp::cmat FlatGSE(int d){ // {{{
//     itpp::Mat<std::complex<double> >& U, itpp::Vec<double>& eigenvalues){
    itpp::cmat U(d,d);
    itpp::vec eig(d);
    FlatSpectrumGSE(U,eig);
    return U*diag(eig)*U.hermitian_transpose();
  }
  // }}}
  void FlatSpectrumGSE(itpp::Mat<std::complex<double> >& U, itpp::Vec<double>& eigenvalues){ // {{{
    int dim=eigenvalues.size();
    itpp::Mat<std::complex<double> > temp(dim,dim);
    double pi=3.1415926535897932384626433832795;
    bool exito=false;
    for( ; !exito ; ){
      temp=RandomGSEDeltaOne(dim);
      itpp::eig_sym(temp, eigenvalues, U);
      eigenvalues*=sqrt(1/(8.*dim));
      exito=unfoldcircular(eigenvalues);
    }
  }
  // }}}
  void FlatSpectrumGSE(itpp::Vec<double>& eigenvalues){ // {{{
    itpp::Mat<std::complex<double> > U;
    FlatSpectrumGSE(U,eigenvalues);
  } // }}}
  itpp::vec FlatSpectrumGSE(int d){ // {{{
    itpp::vec eigenvalues(d);
    FlatSpectrumGSE(eigenvalues);
    return eigenvalues; 
  }
  // }}}
// }}}
// GOE Collection {{{
  itpp::mat RandomGOEDeltaOne(int const dim){ //{{{
    itpp::mat temp(dim, dim);
    temp=itpp::randn(dim,dim);
    return sqrt(1/2.)*(temp+temp.hermitian_transpose());
  }
  // }}}
  itpp::Vec<double> FlatSpectrumGOE(int const dim, double const percen_out){ //{{{
    int extra_per_side=itpp::round_i(std::ceil(0.5*percen_out*dim));
    int dim_large=itpp::round_i(std::ceil(dim+2*extra_per_side));
    itpp::mat temp(dim_large,dim_large);
    itpp::Vec<double> eigenvalues(dim);
    itpp::Vec<double> eigenvalues_enlarged(dim_large);
    // 		eigenvalues.set_size(dim_enlarged);
    temp=RandomGOEDeltaOne(dim_large);
    itpp::eig_sym(temp, eigenvalues_enlarged);
    eigenvalues_enlarged /= sqrt(4*dim_large);
    if (eigenvalues_enlarged(extra_per_side)<-1 ||
        eigenvalues_enlarged(dim+extra_per_side-1)> 1){
      std::cout<<"UUps, we need to a bigger part of the "
        <<" spectrum in routine FlatSpectrumGUE";
      exit(1);
    }
    unfoldcircular(eigenvalues_enlarged);
    eigenvalues=eigenvalues_enlarged.get(extra_per_side,
        dim+extra_per_side-1);
    return eigenvalues;
  }
  // }}}
  itpp::mat RandomGOE(int const dim, std::string normalization="sigma_offdiag=1", double const percentage_away=0.1){ //{{{
    if (normalization=="sigma_offdiag=1"){
      return RandomGOEDeltaOne(dim);
    } else {
      std::cerr  << "Illegal normalization RandomGOE" << percentage_away;
      exit(1);
    }	
  }
  // }}}
  // }}}
// PUE Collection {{{
  itpp::vec RandomPUEspectrum(int const dim, double eigen){ //{{{
    itpp::vec temp(dim);
    temp=2*eigen*(itpp::randu(dim)-0.5);
    return temp;
  }
  // }}}
  itpp::cmat RandomPUE(int const dim, double eigen){ //{{{
    itpp::vec lambda=RandomPUEspectrum(dim, eigen);
    itpp::cmat U=RandomCUE(dim);
    return itppextmath::UDiagUdagger(U, lambda);
//     temp=itpp::randu(dim)-0.5;
//     return temp;
  }
  itpp::vec RandomPUEspectrum_G(int const dim, double eigen){ //{{{
    itpp::vec temp(dim);
    temp=eigen*(itpp::randn(dim));//Gaussian distribution with sigma = eigen
    return temp;
  }
  // }}}
  itpp::cmat RandomPUE_G(int const dim, double eigen){ //{{{
    itpp::vec lambda=RandomPUEspectrum_G(dim, eigen);
    itpp::cmat U=RandomCUE(dim);
    return itppextmath::UDiagUdagger(U, lambda);
  }
  // }}}
  // }}}
// GUE, CUE Collection {{{
  itpp::Mat<std::complex<double> > RandomGUEDeltaOne(int const dim){ //{{{
    itpp::Mat<std::complex<double> > temp(dim, dim);
    temp=itpp::randn_c(dim,dim);
    return sqrt(1/2.)*(temp+temp.hermitian_transpose());
  }
  // }}}
  void FlatSpectrumGUE(itpp::Mat<std::complex<double> >& U, itpp::Vec<double>& eigenvalues){ // {{{
    int dim=eigenvalues.size();
    if(dim != U.rows() || dim != U.cols()){
      std::cerr << "error, dimensionenes no apropiadas "
        << " en FlatSpectrumGUE" <<std::endl  ;
      exit(1);
    }
    itpp::Mat<std::complex<double> > temp(dim,dim);
    bool exito=false;
    for( ; !exito ; ){
      temp=RandomGUEDeltaOne(dim);
      itpp::eig_sym(temp, eigenvalues, U);
      eigenvalues/=sqrt(4*dim);
      exito=unfoldcircular(eigenvalues);
    }
  }
  // }}}
  void FlatSpectrumGUE(itpp::Vec<double>& eigenvalues){ // {{{
    itpp::Mat<std::complex<double> > U;
    FlatSpectrumGUE(U,eigenvalues);
  }
  // }}}
  itpp::Vec<double> FlatSpectrumGUE(int const dim, double const percen_out=0.2){ //{{{
    int extra_per_side=itpp::round_i(std::ceil(0.5*percen_out*dim));
    int dim_large=itpp::round_i(std::ceil(dim+2*extra_per_side));
    itpp::Mat<std::complex<double> > temp(dim_large,dim_large);
    itpp::Vec<double> eigenvalues(dim);
    itpp::Vec<double> eigenvalues_enlarged(dim_large);
    // 		eigenvalues.set_size(dim_enlarged);
    temp=RandomGUEDeltaOne(dim_large);
    itpp::eig_sym(temp, eigenvalues_enlarged);
    eigenvalues_enlarged /= sqrt(4*dim_large);
    if (eigenvalues_enlarged(extra_per_side)<-1 ||
        eigenvalues_enlarged(dim+extra_per_side-1)> 1){
      std::cout<<"UUps, we need to a bigger part of the "
        <<" spectrum in routine FlatSpectrumGUE";
      exit(1);
    }
    unfoldcircular(eigenvalues_enlarged);
    eigenvalues=eigenvalues_enlarged.get(extra_per_side,
        dim+extra_per_side-1);
    return eigenvalues;
  }
  // }}}
  itpp::Mat<std::complex<double> > RandomGUE(int const dim, std::string normalization="sigma_offdiag=1", double const percentage_away=0.1){ //{{{
    if (normalization=="sigma_offdiag=1"){
      return RandomGUEDeltaOne(dim);
    }	
    else if (normalization=="unfolded mean_level_spacing=1"){
      itpp::Mat<std::complex<double> > U(dim, dim), tmp(dim,dim);
      itpp::Vec<std::complex<double> > vec1(dim);
      itpp::Vec<double> eigenvalues(dim);
      FlatSpectrumGUE(U, eigenvalues);
      for (int i=0; i<dim; i++){
        vec1=itpp::elem_mult(conj(U.get_col(i)), to_cvec(eigenvalues));
        for (int j=i; j<dim; j++){
          tmp(i,j)=vec1*U.get_col(j);
          if (i<j){tmp(j,i)=conj(tmp(i,j));}
        }
      }
      return tmp;
      // std::cout  << eigenvalues << 	std::endl ;
      // ya dentro de temp tenemos a una matriz que no esta 
      // unfolded. tengo que encontrar los eigenvectores,
      // luego los eigenvalores, eso diagonalizarlos y chan 
    } else {
      std::cerr  << "Illegal normalization RandomGUE" << percentage_away;
      exit(1);
    }
    // Aca poner un factor de normalizacion opcional
  }
  // }}}
  itpp::Mat<std::complex<double> > RandomCUE(int const dim){ //{{{
    itpp::Mat<std::complex<double> > H;
    H = RandomGUEDeltaOne(dim);

    itpp::Mat<std::complex<double> > U(dim,dim);
    itpp::Vec<double> eigenvalues(dim);
    eig_sym(H, eigenvalues, U);
//     FlatSpectrumGUE(U, eigenvalues);
    //     return exp( 2*itpp::pi*std::complex<double>(0.,1.)*  itpp:randu())*U;
    return exp(2.*itpp::pi*std::complex<double>(0.,1.) *itpp::randu() )*U;
  }
  // }}}
  // }}}
// Checking things {{{
  void PrintElementsGUE(){ // {{{
    int dim=4;
#ifdef ASK
    std::clog  << "Inserte dim "; cin >> dim;
#endif
    std::clog  << "dim="<<dim;
    itpp::cmat H;
    for (int n_ensemble=0;n_ensemble<100;n_ensemble++){
      H=RandomGUE(dim);
      for(int i=0;i<dim;i++){
        std::cout  << real(H(i,i)) << std::endl ;
        for(int j=i+1;j<dim;j++){
          std::cout <<H(j,i)<<std::endl ;
        }
      }

    }
  }
  // }}}
  void CheckElementsGUE(){ // {{{
    int dim=40, ElementsEnsemble=1000;
#ifdef ASK
    std::clog  << "Inserte dim "; cin >> dim; std::cout <<std::endl ;
    std::clog  << "Inserte ElementsEnsemble "; cin >> ElementsEnsemble;std::cout <<std::endl ;
#endif
    std::clog  << "dim="<<dim<<" ElementsEnsemble="<<ElementsEnsemble<<std::endl ;
    itpp::cmat H;
    double DiagonalAverage=0, DiagonalSigma=0,OffDiagonalSigma=0;
    std::complex<double>  OffDiagonalAverage=0;
    for (int iEnsemble=0; iEnsemble<ElementsEnsemble;iEnsemble++){
      H=RandomGUE(dim);
      for(int i=0;i<dim;i++){
        DiagonalAverage+=real(H(i,i));
        DiagonalSigma+=pow(abs(H(i,i)),2);
        for(int j=i+1;j<dim;j++){
          // 				OffDiagonalAverage+=H(i,j);
          OffDiagonalAverage+= H(i,j) ;
          // 				OffDiagonalSigma+=pow(abs(real(H(i,j))),2)
          // 				  +complex(0,1)*pow(abs(imag(H(i,j))),2);
          OffDiagonalSigma+=pow(abs(H(i,j)),2);
        }
      }
    }
    DiagonalAverage/=(dim*ElementsEnsemble);
    DiagonalSigma/=(dim*ElementsEnsemble);
    OffDiagonalAverage/=(0.5*dim*(dim-1)*ElementsEnsemble);
    OffDiagonalSigma/=(0.5*dim*(dim-1)*ElementsEnsemble);
    std::clog  << "Aca se muestra que la matriz cumple con las relaciones de"<<std::endl ;
    std::clog  << "\\bar{V_ij} = 0 " <<std::endl ;
    std::clog  << "\\bar{V_ij V_kl} = \\delta_il \\delta_jk " <<std::endl ;
    std::clog  << "\\bar{|V_ij|^2} = 1 " <<std::endl ;
    std::clog  << "tando para elementos diagonales como no diagonales. " <<std::endl ;
    std::clog  << "Como sigue, en orden ,deben ser 0, 1, 0 , 0, 1 " <<std::endl ;
    std::cout  << DiagonalAverage <<std::endl ;
    std::cout  << DiagonalSigma   <<std::endl ;
    std::cout  << OffDiagonalAverage.real() <<std::endl ;
    std::cout  << OffDiagonalAverage.imag() <<std::endl ;
    std::cout  << OffDiagonalSigma<<std::endl ;


  }
  // }}}
  void CheckUnfolded(){ // {{{

    std::clog  << "Esta rutina toca que tenga la lindea donde digo ... hist.update(eigenvalues);... "
      << "en el loop. Ahora la comento porque se enchocha haciendo la optimizacion -O2"
      << std::endl ;
    int dim=40, ElementsEnsemble=1000;
#ifdef ASK
    std::clog  << "Inserte dim "; cin >> dim; std::cout <<std::endl ;
    std::clog  << "Inserte ElementsEnsemble "; cin >> ElementsEnsemble;std::cout <<std::endl ;
#endif
    std::clog  << "dim="<<dim<<" ElementsEnsemble="<<ElementsEnsemble<<std::endl ;
    itpp::cmat H;
    itpp::Vec<double> eigenvalues(dim);
    itpp::Histogram<double> hist(-0.5*dim+0.5, 0.5*dim-.5, dim);
    for (int i=0; i<ElementsEnsemble; i++){
      FlatSpectrumGUE(eigenvalues);
      //  		hist.update(eigenvalues);
    }
    std::cout  << hist.get_bins() - ElementsEnsemble <<std::endl ;
    std::cout  << "los elementos deben ser muy chicos con respecto a "
      <<ElementsEnsemble<<std::endl ;
  } // }}}
// }}}
} // }}}
#endif // RMT_VARIOUS

