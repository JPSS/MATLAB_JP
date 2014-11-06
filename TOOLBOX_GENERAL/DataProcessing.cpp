/*
 *  dataprocessing.cpp
 *  
 *
 *  Created by Jonas Funke on 8/30/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "math.h"
#include "stdlib.h"
#include "mex.h"
#include <vector>
#include "DataProcessing.h"
using namespace std;


//PRIVATE METHODS
int DataProcessing::max(int a, int b){
	if(b>a) return b;
	else return a;
}

int DataProcessing::min(int a, int b){
	if(b<a) return b;
	else return a;
}
 

 


double DataProcessing::theta(double x){
	if(x < 0) {return 0.0;}
	else{ 
		if(x > 0) {return 1.0;}
		else {return 0.5;}
	}
}

double DataProcessing::chi2_bleaching(vector<double> * x, double a, int t0, double bg){ //heigth, t, bg
	double chi2 = 0.0;
	for(int i=0; i < x->size(); i++){
		chi2 += pow(x->at(i)-a*theta(t0-i)-bg*theta(i-t0) , 2);
	}
	return chi2;
}

double DataProcessing::chi2_2step(std::vector< double > * x, int t0, int t1, double I_1, double I_2, double I_bg){ //x, t0, t1, a, b
	return quad_diff(x, I_1, 0, t0) + quad_diff(x, I_2, t0+1, t1) + quad_diff(x, I_bg, t1+1, x->size());
	


}



double DataProcessing::fit_2step(std::vector< double > * x , std::vector< double > * p){//trace, parameters_out
	double I_1, I_1_min, I_2, I_2_min, I_bg, I_bg_min, chi2, chi2_min;
	int t0, t0_min, t1, t1_min;
	chi2_min = 1e20;
	
	for(int t0=0; t0 < x->size()-2; t0++){  // t0 in [0, size(x)-2]
		I_1 = my_avg(x, 0, t0);		//mean to the left
		for(int t1=t0+1; t1 < x->size()-1; t1++){ //t1 in [t0, size(x)-1]
			I_2 = my_avg(x, t0, t1); //mean in the middle
			I_bg = my_avg(x, t1+1, x->size()-1);
			chi2 = chi2_2step(x, t0, t1, I_1, I_2, I_bg);
			if (chi2 < chi2_min) {
				I_1_min = I_1;
				I_2_min = I_2;
				I_bg_min = I_bg;
				t0_min = t0;
				t1_min = t1;
				chi2_min = chi2;
			}
		}
	}
	p->clear();
	p->push_back(t0_min);
	p->push_back(t1_min);
	p->push_back(I_1_min);
	p->push_back(I_2_min);
	p->push_back(I_bg_min);
	return chi2_min;

}

//calculates the average of evctor x from start to stop
double DataProcessing::my_avg(vector< double > * x, int start, int stop){
	if(start < 0 ){ mexPrintf("Error in avg: start index is smaller than 0\n");}
	if(start >= x->size() ){ mexPrintf("Error in avg: start index is larger than size\n");}
	if(stop < 0 ){ mexPrintf("Error in avg: stop index is smaller than 0\n");}
	if(start >= x->size() ){ mexPrintf("Error in avg: stop index is lagrer than size\n");} 
	double avg = 0.0; int n=0;
	for(int i=max(start,0); i <= min(x->size()-1, stop); i++){
		avg += x->at(i); n++;
	}
	avg = avg / n;
	return avg;
}

//calculates the variance of x from start to stop
double DataProcessing::my_var(vector< double > * x, double avg, int start, int stop){
	if(start < 0 ){ mexPrintf("Error in var: start index is smaller than 0\n");}
	if(start >= x->size() ){ mexPrintf("Error in var: start index is larger than size\n");}
	if(stop < 0 ){ mexPrintf("Error in var: stop index is smaller than 0\n");}
	if(start >= x->size() ){ mexPrintf("Error in var: stop index is lagrer than size\n");} 
	
	double sum2 = 0.0;
	double sum3 = 0.0;
	int n=0;
	for(int i=max(0, start); i <= min(x->size()-1, stop); i++){
		sum2 += pow(x->at(i) - avg , 2);
		sum3 += (x->at(i) - avg);
		n++;
	}
	return (sum2 - pow(sum3,2)/n) /(n - 1);
}

double DataProcessing::quad_diff(vector<double> * x, double avg, int start, int stop){
	if(start < 0 ){ mexPrintf("Error in var: start index is smaller than 0\n");}
	if(start >= x->size() ){ mexPrintf("Error in var: start index is larger than size\n");}
	if(stop < 0 ){ mexPrintf("Error in var: stop index is smaller than 0\n");}
	if(start >= x->size() ){ mexPrintf("Error in var: stop index is lagrer than size\n");} 	
	int n=0;
	double qd = 0.0;
	for(int i=max(0, start); i <= min(x->size()-1, stop); i++){
		qd += pow(avg - x->at(i), 2);
		n++;
	}
	return qd;
}


//END PRIVATE METHODS


//PUBLIC METHODS
double DataProcessing::fit_bleaching(std::vector< double > * x, std::vector< double > * parameters ){ //trace, parameters
	double a, bg, t, a_min, bg_min, t_min, chi2, chi2_min;
	//initial 
	a_min = my_avg(x, 0, 0);
	bg_min = my_avg(x, 0, x->size()-1);
	t_min = 0;
	chi2_min = chi2_bleaching(x, a_min, 0, bg_min);
	
	for (int i=1; i < x->size()-1; i++) {
		a = my_avg(x, 0, i-1);
		bg = my_avg(x, i+1, x->size()-1);
		chi2 = chi2_bleaching(x, a, i, bg);
		if (chi2 < chi2_min) {
			a_min = a;
			t_min = i;
			bg_min = bg;
			chi2_min = chi2;
		}
		//mexPrintf("%i %e %e %e\n",i, a, bg, chi2);
	}
	
	parameters->clear();
	parameters->push_back(a_min);
	parameters->push_back(t_min);
	parameters->push_back(bg_min);
	return chi2_min;
}

void DataProcessing::nlf(vector< double > * x, vector< double > * y , vector<int> * win_size, int p){
	int windows = win_size->size();
	double norm, varf, varb;
	
	//make sure y has same length as x
	while(y->size() > x->size()){
		y->pop_back();
	}
	while(y->size() < x->size() ){
		y->push_back(0.0);
	}
	
	//coefficient/weights
	double f[windows];
	double b[windows];
	double avg_f[windows];
	double avg_b[windows];
	double var_f[windows];
	double var_b[windows];
	
	for(int i=0; i < x->size(); i++){ //loop over data
		//init weights
		for(int j=0; j <  windows; j++){ 
			f[j] = 0; 
			b[j] = 0; 
			avg_b[j] =0; 
			avg_f[j] = 0;
			var_b[j] =0; 
			var_f[j] = 0;
		} 
		
		//determine f ,b weights
		for(int j=0; j <  windows; j++){
			if(i - win_size->at(j) >= 0){ //only use this weight if left index is >= 0 otherwise it will be 0
				avg_f[j] = my_avg(x, i-win_size->at(j), i);  //rethink that ===> boundaries??????
				var_f[j] = my_var(x, avg_f[j], i-win_size->at(j), i);	
				f[j] = pow(var_f[j],-1.*p);
				//mexPrintf("Forward: %e\t%e\t%e\n", pow(1./varf, p), 1./pow(varf,p), f[j]);
			}
			
			if(i+win_size->at(j) < x->size()){
				avg_b[j] = my_avg( x, i, i+win_size->at(j) );  //rethink that ===> boundaries??????
				var_b[j] = my_var(x, avg_b[j], i, i+win_size->at(j));
				b[j] = pow(var_b[j],-1.*p);
				//mexPrintf("i=%i win_size=%i Back: avg_b=%.4f varb=%.4f b=%.4f\n", i, win_size->at(j), avg_b[j], varb, b[j]);
			}
		}
		//normalize weights
		norm = 0.0;
		for(int j=0; j <  windows; j++){
			norm += f[j]+b[j];
		}
		
		/*for(int j=0; j < windows; j++){
			norm = 0.0;			
			for(int k=0; k < windows; k++){
				norm += pow(var_f[j]/var_f[k],p)+ pow(var_f[j]/var_b[k],p);
			}
			f[j] = 1./norm;
			norm = 0.0;			
			for(int k=0; k < windows; k++){
				norm += pow(var_b[j]/var_f[k],p)+ pow(var_b[j]/var_b[k],p);
			}
			b[j] = 1./norm;
			
		}*/
		//mexPrintf("norm =%e\t", norm);
		y->at(i) = 0.0;	
		for(int j=0; j <  windows; j++){
			y->at(i) += (f[j]*avg_f[j] + b[j] * avg_b[j])/norm;  //f[j]*avg_f[j] + b[j] * avg_b[j];
		}
		
		
		
	} // end i 
	
	
}



void DataProcessing::nlf_fret(vector< double > * x, vector<double> * fret,  vector< double > * donor_out, vector< double > * fret_out , vector<int> * win_size, int p){
	int windows = win_size->size();
	double norm, varf, varb, dout,  fout;
	int N = min(x->size(), fret->size());
	
	donor_out->clear();
	fret_out->clear();
	
	//coefficient/weights
	double f[windows];
	double b[windows];
	double avg_f_x[windows];
	double avg_b_x[windows];
	double avg_f_fret[windows];
	double avg_b_fret[windows];
	
	
	for(int i=0; i < N; i++){ //loop over data
		//init weights
		for(int j=0; j <  windows; j++){ 
			f[j] = 0; 
			b[j] = 0; 
			avg_b_x[j] =0; 
			avg_f_x[j] = 0;
			avg_b_fret[j] =0; 
			avg_f_fret[j] = 0;
		} 
		
		//determine f ,b weights
		for(int j=0; j <  windows; j++){
			if(i - win_size->at(j) >= 0){ //only use this weight if left index is >= 0 otherwise it will be 0
				avg_f_x[j] = my_avg(x, i-win_size->at(j), i);  
				avg_f_fret[j] = my_avg(fret, i-win_size->at(j), i);  

				varf = my_var(x, avg_f_x[j], i-win_size->at(j), i) + my_var(fret, avg_f_fret[j], i-win_size->at(j), i);	
				f[j] = pow(varf,-1.*p);
			}
			
			if(i+win_size->at(j) < N){
				avg_b_x[j] = my_avg( x, i, i+win_size->at(j) );  //rethink that ===> boundaries??????
				avg_b_fret[j] = my_avg( fret, i, i+win_size->at(j) );
				varb = my_var(x, avg_b_x[j], i, i+win_size->at(j)) + my_var(fret, avg_b_fret[j], i, i+win_size->at(j));
				b[j] = pow(varb,-1.*p);
			}
		}
		//normalize weights
		norm = 0.0;
		for(int j=0; j <  windows; j++){
			norm += f[j]+b[j];
		}
		//calc output
		dout  = 0.0;	
		fout = 0.0;	
		for(int j=0; j <  windows; j++){
			dout += (f[j]*avg_f_x[j] + b[j] * avg_b_x[j])/norm;
			fout += (f[j]*avg_f_fret[j] + b[j] * avg_b_fret[j])/norm;
		}
		donor_out->push_back(dout);
		fret_out->push_back(fout);
		
		
		
	} // end i 
	
	//pad data to fill up to original size
	double last = donor_out->back();
	while(donor_out->size() < x->size()){
		donor_out->push_back(last);
	}
	last = fret_out->back();
	while(fret_out->size() < fret->size()){
		fret_out->push_back(last);
	}
	
}





void DataProcessing::nlf_alex(vector< double > * donor, vector<double> * acceptor ,vector<double> * fret,  vector< double > * donor_out, vector<double> * acceptor_out, vector< double > * fret_out , vector<int> * win_size, int p){
	int windows = win_size->size();
	double norm, varf, varb, dout,aout, fout;
	int N = min(acceptor->size(), min(donor->size(), fret->size()));
	
	donor_out->clear();
	acceptor_out->clear();
	fret_out->clear();
	
	//coefficient/weights
	double f[windows];
	double b[windows];
	double avg_f_donor[windows];
	double avg_b_donor[windows];
	double avg_f_acceptor[windows];
	double avg_b_acceptor[windows];
	double avg_f_fret[windows];
	double avg_b_fret[windows];
	
	
	for(int i=0; i < N; i++){ //loop over data
		//init weights
		for(int j=0; j <  windows; j++){ 
			f[j] = 0; 
			b[j] = 0; 
			avg_b_donor[j] =0; 
			avg_f_donor[j] = 0;
			avg_b_acceptor[j] =0; 
			avg_f_acceptor[j] = 0;
			avg_b_fret[j] =0; 
			avg_f_fret[j] = 0;
		} 
		
		//determine f ,b weights
		for(int j=0; j <  windows; j++){
			if(i - win_size->at(j) >= 0){ //only use this weight if left index is >= 0 otherwise it will be 0
				avg_f_donor[j] = my_avg(donor, i-win_size->at(j), i);
				avg_f_acceptor[j] = my_avg(acceptor, i-win_size->at(j), i);  
				avg_f_fret[j] = my_avg(fret, i-win_size->at(j), i);  

				varf = my_var(donor, avg_f_donor[j], i-win_size->at(j), i) + my_var(acceptor, avg_f_acceptor[j], i-win_size->at(j), i)+ my_var(fret, avg_f_fret[j], i-win_size->at(j), i);	
				f[j] = pow(varf,-1.*p);
			}
			
			if(i+win_size->at(j) < N){
				avg_b_donor[j] = my_avg( donor, i, i+win_size->at(j) );  //rethink that ===> boundaries??????
				avg_b_acceptor[j] = my_avg( acceptor, i, i+win_size->at(j) ); 
				avg_b_fret[j] = my_avg( fret, i, i+win_size->at(j) );
				varb = my_var(donor, avg_b_donor[j], i, i+win_size->at(j)) + my_var(acceptor, avg_b_acceptor[j], i, i+win_size->at(j)) + my_var(fret, avg_b_fret[j], i, i+win_size->at(j));
				b[j] = pow(varb,-1.*p);
			}
		}
		//normalize weights
		norm = 0.0;
		for(int j=0; j <  windows; j++){
			norm += f[j]+b[j];
		}
		//calc output
		dout  = 0.0;
		aout = 0.0;	
		fout = 0.0;	
		for(int j=0; j <  windows; j++){
			dout += (f[j]*avg_f_donor[j] + b[j] * avg_b_donor[j])/norm;
			aout += (f[j]*avg_f_acceptor[j] + b[j] * avg_b_acceptor[j])/norm;
			fout += (f[j]*avg_f_fret[j] + b[j] * avg_b_fret[j])/norm;
		}
		donor_out->push_back(dout);
		acceptor_out->push_back(aout);
		fret_out->push_back(fout);
		
		
		
	} // end i 
	

	//pad data to fill up to original size
	double last = donor_out->back();
	while(donor_out->size() < donor->size()){
		donor_out->push_back(last);
	}
	last = acceptor_out->back();
	while(acceptor_out->size() < acceptor->size()){
		acceptor_out->push_back(last);
	}
	last = fret_out->back();
	while(fret_out->size() < fret->size()){
		fret_out->push_back(last);
	}
	
}







int DataProcessing::fit_traces(std::vector<double> * d , std::vector<double> * a   ,std::vector<double> * f , std::vector<double> * parameters  ){
	int ta_min, td_min;
	double E, E_min, E_a, E_d,E_f, avg_l_a, avg_l_d, avg_l_f, avg_m_d, avg_r_a, avg_r_f, avg_r_d;
	int N = min(min(a->size(), d->size()), f->size());
	int found = -1;
	int found_min = found;
	
	//initial value 
	avg_l_a = my_avg(a, 0, 0);
	avg_r_a = my_avg(a, 0, a->size()-1);
	E_a = quad_diff(a,avg_l_a, 0, 0 ) + quad_diff(a, avg_r_a , 0, a->size()-1);
	//donor		
	avg_l_d = my_avg(d, 0, 0);
	avg_r_d = my_avg(d, 0, d->size()-1);
	E_d = quad_diff(d,avg_l_d, 0, 0 ) + quad_diff(d, avg_r_d , 0, d->size()-1); 
	//fret
	avg_l_f = my_avg(f, 0, 0);
	avg_r_f = my_avg(f, 0, f->size()-1);
	E_f = quad_diff(f ,avg_l_f, 0, 0 ) + quad_diff(f, avg_r_f , 0, f->size()-1);
	//overall energy
	if(pow(avg_l_a-avg_r_a,2) > 0 && pow(avg_l_d-avg_r_d,2) > 0 && pow(avg_l_f-avg_r_f,2)){
		//E_min = E_a/pow(avg_l_a-avg_r_a,2) + E_d/pow(avg_l_d-avg_r_d,2) + E_f/pow(avg_l_f-avg_r_f,2);		
		E_min = E_a+E_f+E_d;
	}
	else {
		E_min = 1e15;
		mexPrintf("WARNING: Initialization of E_min not possible...E_min=%e\n", E_min);
	}

	
	
	
	
	for (int ta=0; ta < N-1 ; ta++) {
		for (int td=0; td < N-1 ; td++) {
			if(td > ta){ //acceptor bleaches first
				//acceptor
				avg_l_a = my_avg(a, 0, ta);
				avg_r_a = my_avg(a, ta+1, a->size()-1);
				E_a = quad_diff(a,avg_l_a, 0, ta ) + quad_diff(a, avg_r_a , ta+1, a->size()-1);
				//fret
				avg_l_f = my_avg(f, 0, ta);
				avg_r_f = my_avg(f, ta+1, f->size()-1);
				E_d = quad_diff(f ,avg_l_f, 0, ta ) + quad_diff(f, avg_r_f , ta+1, f->size()-1);
				//donor
				avg_l_d = my_avg(d, 0, ta);
				avg_m_d = my_avg(d, ta+1, td);
				avg_r_d = my_avg(d, td+1, d->size()-1);
				E_f = quad_diff(d ,avg_l_d, 0, ta ) + quad_diff(d, avg_m_d , ta+1, td) + quad_diff(d, avg_r_d , td+1, d->size()-1);
				if(pow(avg_l_a-avg_r_a,2) > 0 && pow(avg_l_f-avg_r_f,2) > 0 && pow(avg_l_d-avg_r_d,2)){
					//E = E_a/pow(avg_l_a-avg_r_a,2) + E_f/pow(avg_l_f-avg_r_f,2) + E_d/pow(avg_l_d-avg_r_d,2);
				}
				else {
					mexPrintf("WARNING: Differences are zero.\n");
				}
				E = E_a +E_d + E_f;
				found = 0;
				
			}
			else{ //donor bleaches first
				//acceptor		
				avg_l_a = my_avg(a, 0, ta);
				avg_r_a = my_avg(a, ta+1, a->size()-1);
				E_a = quad_diff(a,avg_l_a, 0, ta ) + quad_diff(a, avg_r_a , ta+1, a->size()-1);
				//donor		
				avg_l_d = my_avg(d, 0, td);
				avg_r_d = my_avg(d, td+1, d->size()-1);
				E_d = quad_diff(d,avg_l_d, 0, td ) + quad_diff(d, avg_r_d , td+1, d->size()-1); 
				//fret
				avg_l_f = my_avg(f, 0, td);
				avg_r_f = my_avg(f, td+1, f->size()-1);
				E_f = quad_diff(f ,avg_l_f, 0, td ) + quad_diff(f, avg_r_f , td+1, f->size()-1);
				//overall energy
				if(pow(avg_l_a-avg_r_a,2) > 0 && pow(avg_l_d-avg_r_d,2) > 0 && pow(avg_l_f-avg_r_f,2)){
					//E = E_a/pow(avg_l_a-avg_r_a,2) + E_d/pow(avg_l_d-avg_r_d,2) + E_f/pow(avg_l_f-avg_r_f,2);	
				}
				else {
					mexPrintf("WARNING: Differences are zero.\n");
				}
				E = E_a +E_d + E_f;
				found = 1;
			}
			//mexPrintf("%i\t%e\n", found, E);
			if (E < E_min) {
				E_min = E;
				found_min = found;
				ta_min = ta;
				td_min = td;
				//mexPrintf("%e\t", E );
				//mexPrintf("%i\t%i\t%i\n", td, ta, found);
			}
			
			
		}
	}
	
	
	parameters->clear();
	parameters->push_back(E_min);
	parameters->push_back(td_min);
	parameters->push_back(ta_min);
	if (found_min == 0) { //acceptor bleaches first
		parameters->push_back(ta_min); //t_f
		parameters->push_back(my_avg(d, 0, ta_min)); //avg_l_d
		parameters->push_back(my_avg(d, ta_min, td_min )); //avg_m_d
		parameters->push_back(my_avg(d, td_min, d->size()-1)); //avg_r_d
		
		parameters->push_back(my_avg(a, 0, ta_min)); //avg_l_a
		parameters->push_back(my_avg(a, ta_min, a->size()-1)); //avg_r_a
		
		parameters->push_back(my_avg(f, 0, ta_min)); //avg_l_f
		parameters->push_back(my_avg(f, ta_min, f->size()-1)); //avg_r_f
	}
	if (found_min == 1) {//donor bleaches first
		parameters->push_back(td_min); //t_f
		parameters->push_back(my_avg(d, 0, td_min)); //avg_l_d
		parameters->push_back(my_avg(d, td_min, d->size()-1)); //avg_r_d
		
		parameters->push_back(my_avg(a, 0, ta_min)); //avg_l_a
		parameters->push_back(my_avg(a, ta_min, a->size()-1)); //avg_r_a
		
		parameters->push_back(my_avg(f, 0, td_min)); //avg_l_f
		parameters->push_back(my_avg(f, td_min, f->size()-1)); //avg_r_f
		
	}
	
	/*mexPrintf("%i\t", found_min);
for (int i=0; i<parameters->size(); i++) {
	mexPrintf("%.2f\t", parameters->at(i));
}*/
	
	return found_min;
	
	
}


 void DataProcessing::find_peaks1d(std::vector< double > * x, int width, double h_min, int absolute, std::vector< double > * peaks){
    double local_max; 
    double abs_min;
    double local_avg;
    
    if (x->size() > 0){
        abs_min = x->at(0);
    }
    
            
    int j_max;
    peaks->clear();
    
    for(int i=0; i < x->size(); i++){
        if(abs_min > x->at(i)){
                abs_min = x->at(i);
        } 
    }
    
    //loop through data and find candidates
     for(int i=0; i < x->size(); i++){
        local_max = abs_min;
        j_max = -1;

        if (absolute==1) {
            local_avg = 0;//my_avg(x, max(0,i-4*width), min(x->size()-1, i+4*width)); //local average with respect to 4 * width
        }
        else {
            local_avg = my_avg(x, max(0,i-2*width), min(x->size()-1, i+2*width)); //local average with respect to 2 * width
        }       
        for(int j=max(0,i-width); j <= min(x->size()-1, i+width) ;j++){
            if(x->at(j) > local_max){
                local_max = x->at(j);
                j_max = j;
            }
        }
        
        
        if ( (local_max == x->at(i)) && (j_max == i) && (local_max-local_avg >= h_min)){
            peaks->push_back(i);
        }
        
     }
         
     
 }








void DataProcessing::studentT(std::vector< double > * x, std::vector< double > * x_out, int N ){
	double t, avg_l, avg_r, var_l, var_r;
	x_out->clear();
	if (x->size() >= 2*N-1 && N>1) {
		//first couple of values
		for (int i=0; i < N-1; i++) {
			x_out->push_back(0.);
		}
		//compute measure
		for (int i=N-1; i <= x->size()-N; i++) {
			avg_l = my_avg(x, i-N+1, i);
			var_l = my_var(x, avg_l, i-N+1,i);
			
			avg_r = my_avg(x, i, i+N-1);
			var_r = my_var(x, avg_r,i, i+N-1);

			t = (avg_r-avg_l);// / sqrt(var_l/N + var_r/N);
			
			x_out->push_back(t);
		}
		for (int i=x->size()-N+1; i < x->size(); i++) {
			x_out->push_back(0.);
		}
	}
	else {
		for (int i=0; i < x->size() ; i++) {
			x_out->push_back(0.);
		}
	}

	
	
	
	
	
}




