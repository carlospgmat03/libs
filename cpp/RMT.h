#ifndef RMT_H
#define RMT_H
#include <itpp/itbase.h>
#include <itpp/stat/histogram.h>

#include <cmath>
#include <stdio.h>


namespace RMT{
	bool unfoldcircular(itpp::Vec<double>& eigenvalues);
	itpp::Mat<std::complex<double> > RandomGUEDeltaOne(int const dim);
}

#endif
