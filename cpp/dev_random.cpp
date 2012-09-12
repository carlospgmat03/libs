/**
 * \file random.h
 * \author Guy Rutenberg <guyrutenberg@gmail.com>
 * \version 1.0
 * \date 2007
 *
 * from http://www.guyrutenberg.com/2007/09/15/random-a-random-number-generator-class/
 *
 * Header file for Random.
 *
 * License: The MIT License (included in the file).
*/
 
/*****************************************************************************
The MIT License
 
Copyright (c) 2007 Guy Rutenberg
 
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
 
The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.
 
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*****************************************************************************/
 
#ifndef __RANDOM_H__
#define __RNADOM_H__
 
#define DEVRANDOMERROR 1;
 
#include <iostream>
#include <fstream>
using namespace std;
 
/**
 * \brief Interface for strong integer random number generators.
 *
 * Random number generators provided by this class are base on
 * /dev/random and /dev/urandom.
*/
class Random {
private:
	ifstream m_random;
	ifstream m_urandom;
 
public:
	inline Random();
	inline unsigned int secure();
	inline unsigned int strong();
	inline double secure_real();
	inline double strong_real();
	inline unsigned int secure_range(unsigned int x);
	inline unsigned int strong_range(unsigned int x);
};
 
inline Random::Random()
{
	m_random.open("/dev/random",ios::in|ios::binary);
	m_urandom.open("/dev/urandom",ios::in|ios::binary);
	if (!m_random || !m_urandom)
		throw DEVRANDOMERROR;
}
 
/**
 * \brief Generates secure, cryptography strong random numbers.
 *
 * Uses /dev/random to generate cryptography strong random numbers.
 * \note This function may take a long time to return if system entropy levels
 * are low as /dev/random will block.
 * \return randomly generated unsigned int.
*/
inline unsigned int Random::secure()
{
	unsigned int num = 0;
	m_random.read((char*)&num, sizeof(unsigned int));
	return num;
}
 
/**
 * \brief Generates strong pseudo random numbers.
 *
 * Uses /dev/urandom to generate strong random numbers. Unlike
 * Random::secure() this function won't block if entropy levels are low but
 * instead it will generate random numbers with lower randomness level.
*/
inline unsigned int Random::strong()
{
	unsigned int num = 0;
	m_urandom.read((char*)&num, sizeof(unsigned int));
	return num;
}
 
/**
 * \brief Generates secure pseudo random real number.
 *
 * Uses /dev/random to generate secure random real numbers between 0 and 1.
*/
inline double Random::secure_real()
{
	return (double)secure()/(unsigned int)(-1);
}
 
/**
 * \brief Generates strong pseudo random real number.
 *
 * Uses /dev/urandom to generate strong random real numbers between 0 and 1.
*/
inline double Random::strong_real()
{
	return (double)strong()/(unsigned int)(-1);
}
 
/**
 * \brief Generates secure pseudo random number in a limited range.
 *
 * Uses /dev/random to generate secure random real numbers between 0 and 1.
 * \param x [unsigned int] random numbers will be generated between 0 and x-1 (including).
*/
inline unsigned int Random::secure_range(unsigned int x)
{
	// we do mod x to handle the extreme case secure_real() returns 1
	return (unsigned int)(secure_real() * x)%x;
}
 
/**
 * \brief Generates strong pseudo random number in a limited range.
 *
 * Uses /dev/random to generate secure random real numbers between 0 and 1.
 * \param x [unsigned int] random numbers will be generated between 0 and x-1 (including).
*/
inline unsigned int Random::strong_range(unsigned int x)
{
	// we do mod x to handle the extreme case secure_real() returns 1
	return (unsigned int)(strong_real() * x)%x;
}
 
#endif /*__RANDOM_H__*/
