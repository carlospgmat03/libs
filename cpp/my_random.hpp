#ifndef RANDOM_GAUSSIAN
#define RANDOM_GAUSSIAN

#include <iostream>
#include <boost/random.hpp>
#include <complex>

/** \brief Holds functions to get random variables
 *
 * Here we have functions to genereat random gaussian varibles.
 * We also want to test doxygen
 * It relyies on boost library (in particular random), see
 * <A HREF="http://www.boost.org/"> http://www.boost.org/ </A>
 * or <A HREF="http://boost.org/libs/random/random-generators.html"> 
 * http://www.boost.org/ libs/random/random-generators.html</A>.
 *
 * The Random Number generator is set to lagged_fibonacci44497
 * Other valid options are minstd_rand or ecuyer1988 or mt19937 
 *
 * \attention 
 * You have to define "wp", example:
 * \code typedef long double wp; \endcode
 *
 * \attention 
 * you must also have to have available the boost libraries as
 * \code #include <boost/random.hpp> \endcode
 *
 * \author Carlos Pineda
 */
namespace random_variables{
	/** \brief The kind of random number generator to use
	 *
	 * This random number generator can be chosen differently. 
	 * for example choose minstd_rand or ecuyer1988 or mt19937
	 * instead of lagged_fibonacci44497. Upto this point the change
	 * must be done here in the source code.
         * see <A HREF="http://boost.org/libs/random/random-generators.html"> 
         * http://www.boost.org/ libs/random/random-generators.html</A>.
	 *
	 */
	typedef boost::lagged_fibonacci44497 base_generator_type;

	/** \brief contains vital information
	 *
	 * It contains properly the random number generator. A method that
	 * can be accessed is seed. To initialize one can for example use 
	 *   \code global_random_number_generator.seed(42); \endcode
	 *   or
	 *   \code global_random_number_generator(static_cast<unsigned> (std::time(0))) \endcode
	 *
	 */
	static base_generator_type global_random_number_generator(static_cast<unsigned> (std::time(0)));

	/** Function to get gaussian real numbers with a given mean and given standard deviation
	 */
	wp SampleNormal (wp sigma = 1 , wp mean = 0)
	{
    		//select Gaussian probability distribution
		boost::normal_distribution<wp> norm_dist(mean, sigma);
		// bind random number generator to distribution, forming a function
		boost::variate_generator< base_generator_type&, boost::normal_distribution<wp> >  normal_sampler(global_random_number_generator, norm_dist);
		// sample from the distribution
    		return normal_sampler();
	}
	/** Function to get gaussian complex numbers with a given mean and given standard deviation
	 *
	 */
	std::complex<wp> ComplexSampleNormal (  wp sigma = 1  , std::complex<wp> mean = 0 )
	{
		std::complex<wp> number ( SampleNormal ( sigma, std::real(mean)) ,SampleNormal ( sigma, std::imag(mean))  );
		return number;
	}
}
#endif // RANDOM_GAUSSIAN
