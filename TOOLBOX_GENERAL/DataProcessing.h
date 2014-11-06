/*
 *  DataProcessing.h
 *  
 *
 *  Created by Jonas Funke on 8/30/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */



#ifndef DATAPROCESSING_H
#define DATAPROCESSING_H
#include <vector>

class DataProcessing
{
public:
	static void nlf(std::vector< double > * , std::vector< double > *  , std::vector<int> * , int ); //vectzor in, vector out, window sizes, p
	static void nlf_fret(std::vector< double > * , std::vector< double > * , std::vector< double > *,  std::vector< double > *  , std::vector<int> * , int ); //vectzor green, vector fret, donor out, fret out, window sizes, p
	static void nlf_alex(std::vector< double > * , std::vector< double > * ,std::vector< double > *, std::vector< double > *, std::vector< double > *,  std::vector< double > *  , std::vector<int> * , int ); //vectzor donor,,vector acceptor,  vector fret, donor out, acceptor out ,fret out, window sizes, p
	static int max(int ,int);
	static int min(int ,int);
	
	static double fit_bleaching(std::vector< double > * , std::vector< double > *  ); //trace, parameters_out
	static double fit_2step(std::vector< double > * , std::vector< double > *);//trace, parameters_out

	static int fit_traces(std::vector<double> * , std::vector<double> *  ,std::vector<double> * , std::vector<double> * ); //d, a, f, parameter.... returns which option was found

	static void studentT(std::vector< double > * x, std::vector< double > * x_out, int win_size );

    static void find_peaks1d(std::vector< double > * x, int width, double h_min, int absolute, std::vector< double > * peaks);
private:
	static double theta(double x);
	static double my_var(std::vector< double > * , double , int, int );
	static double my_avg(std::vector< double > * , int , int );
	static double chi2_bleaching(std::vector< double > * ,double , int , double); //, x, hight a, t0, bg
	static double chi2_2step(std::vector< double > * , int , int, double, double, double); //x, t0, t1, a, b, bg
	static double quad_diff(std::vector<double> * , double , int , int );

		
	
};

#endif
