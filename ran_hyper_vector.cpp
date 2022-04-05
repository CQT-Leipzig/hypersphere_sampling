#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>

#include <cmath>
#include <vector>
#include <algorithm>
#include <random>

#include "rangen.h"
#include "drand48.h"

#ifndef M_PI
#define M_PI 3.14159265358979324
#endif

//std::mt19937 gen;
std::mt19937_64 gen;

//inline double random_num(void){return my_drand48();}
inline double random_num(void){return static_cast<double>(gen())/gen.max();}
//inline double random_num(void){return rangen.RAN01();}

// a struct allowing the sorting of triples of numbers based on the
// value of the first number
struct triple{
    double S,a,b;
    bool operator<(const triple& other) const { return (S < other.S); }
};

// a struct allowing the sorting of pairs of numbers based on the absolute value of the first number
// this is needed for the in-situ version
// the inteded usage S = a*a + b*b will be modified to S = sign(a)*(a*a + b*b)
// to allow the reconstruction of a from S and b
struct my_pair{
    double S,b;
    bool operator<(const my_pair& other) const { return ( fabs(S) < fabs(other.S) ); } // note the fabs()
};


/**************************************** Mueller's method using Box-Mueller ***************************************/
// This method employs a multivariate normal distribution - initially, each coordinate is a Gaussian random number -
// which is spherically symmetric

void rand_vect_mueller_odd(int n,double* v){
    double r,p,S(0.);
    for(int k=0;k<n-1;k+=2){
        r = sqrt( -2*log(1.0-random_num()) );
        p = random_num();
        v[k] = cos(2*M_PI*p) ;
        v[k+1] = p<0.5 ?  r*sqrt(1.0 -v[k]*v[k]) 
                       : -r*sqrt(1.0 -v[k]*v[k]) ; // = r*sin(2*M_PI*p) ;
        v[k] *= r;
        S += v[k]  *v[k];
        S += v[k+1]*v[k+1];
    }
    r = sqrt( -2*log(1.0-random_num()) );
    p = 2*M_PI*random_num();
    v[n-1] = r * cos(p) ;
    S = sqrt( S + v[n-1]*v[n-1] );
    for(int k=0;k<n;++k)
        v[k] /= S;
}

void rand_vect_mueller_even(int n,double* v){
    double r,p,S(0.);
    for(int k=0;k<n;k+=2){
        r = sqrt( -2*log(1.0-random_num()) );
        p = random_num();
        v[k] = cos(2*M_PI*p) ;
        v[k+1] = p<0.5 ?  r*sqrt(1.0 -v[k]*v[k]) 
                       : -r*sqrt(1.0 -v[k]*v[k]) ;// = r*sin(2*M_PI*p) ;
        v[k] *= r;
        S += v[k]  *v[k];
        S += v[k+1]*v[k+1];
    }
    S = sqrt(S);
    for(int k=0;k<n;++k)
        v[k] /= S;
}

void rand_vect_mueller(int n,double* v){ return n&1 ? rand_vect_mueller_odd(n,v) : rand_vect_mueller_even(n,v) ;}

/*****************************************Basic********************************************************************/

void rand_vect_basic_even(int n,double* v){
    static std::vector< triple > twoDpts;
    if( twoDpts.size() != n/2 )
        twoDpts.resize(n/2 , triple() );
    for(int i=0;i<n/2;++i){
        do{
            twoDpts[i].a = 2*random_num()-1.0;
            twoDpts[i].b = 2*random_num()-1.0;
            twoDpts[i].S = twoDpts[i].a*twoDpts[i].a + twoDpts[i].b*twoDpts[i].b ;
        }
        while( twoDpts[i].S > 1.0 );
    }
    std::sort( twoDpts.begin() , twoDpts.end() );
    double t;
    double s = 1.0 / twoDpts[n/2-1].S ;
    double S_I,S_II(0.0);
    for(int i=0;i<n/2;++i){
        S_I  = S_II ;
        S_II = twoDpts[i].S ;
        t = sqrt( (1.0 - S_I/S_II )*s ) ;
        v[2*i]   = twoDpts[i].a * t ;
        v[2*i+1] = twoDpts[i].b * t ;
    }
    return;
}

void rand_vect_basic_odd(int n,double* v){
    ++n; // starting as if n is even
    static std::vector< triple > twoDpts;
    if( twoDpts.size() != n/2 )
        twoDpts.resize( n/2 , triple() );
    for(int i=0;i<n/2;++i){
        do{
            twoDpts[i].a = 2*random_num()-1.0;
            twoDpts[i].b = 2*random_num()-1.0;
            twoDpts[i].S = twoDpts[i].a*twoDpts[i].a + twoDpts[i].b*twoDpts[i].b ;
        }
        while( twoDpts[i].S > 1.0 );
    }
    std::sort( twoDpts.begin() , twoDpts.end() );
    double t, s( 1.0 / twoDpts[n/2-1].S );
//  in contrast two the even n version we throw away one coordinate ( twoDpts[0].a * sqrt( s ) )
//  and modify s such that the final vector will have length 1    
    s = s / ( 1.0 - s * twoDpts[0].a * twoDpts[0].a );
    double S_I(0.0), S_II(twoDpts[0].S);
    v[0] = twoDpts[0].b * sqrt( s ) ;
    for(int i=1;i<n/2;++i){
        S_I  = S_II ;
        S_II = twoDpts[i].S ;
        t = sqrt( (1.0 - S_I/S_II )*s ) ;
        v[2*i-1] = twoDpts[i].a * t ;
        v[2*i]   = twoDpts[i].b * t ;
    }
    return;
}

// function that is supposed to be called:
void rand_vect_basic(int n,double* v){ return n&1 ? rand_vect_basic_odd(n,v) : rand_vect_basic_even(n,v) ;}

/*****************************************with BucketSort*********************************************************/

const int BUCKET_SIZE=16;

void rand_vect_bs_even(int n,double* v){
    static std::vector< std::vector<triple>  > twoDpts;
    int nb = 1+n/BUCKET_SIZE; //
    twoDpts.resize( nb , std::vector<triple>(0) );
    triple tr;
    for(int i=0;i<n/2;++i){
        do{
            tr.a = 2*random_num()-1.0;
            tr.b = 2*random_num()-1.0;
            tr.S = tr.a*tr.a + tr.b*tr.b ;
        }
        while( tr.S > 1.0 );
        twoDpts[tr.S*nb].push_back( tr );
    }
    for (int bi = 0; bi < nb; ++bi)
        sort(twoDpts[bi].begin(),twoDpts[bi].end());
    
    int last_bucket = nb-1;
    while(!twoDpts[last_bucket ].size())
        --last_bucket;
    double t;
    double S_I,S_II(0.0);
    double s = 1.0 / twoDpts[last_bucket].rbegin()->S ;
    int i=0;
    for (int bi = 0; bi <=last_bucket; ++bi){
        for (int j = 0; j<twoDpts[bi].size(); ++j){
            S_I  = S_II ;
            S_II = twoDpts[bi][j].S ;
            t = sqrt( (1.0 - S_I/S_II )*s ) ;
            v[i++] = twoDpts[bi][j].a * t ;
            v[i++] = twoDpts[bi][j].b * t ;
        }
        twoDpts[bi].resize(0);
    }
    return;
}

void rand_vect_bs_odd(int n,double* v){
    ++n;
    static std::vector< std::vector<triple>  > twoDpts;
    int nb = 1+n/BUCKET_SIZE; //
    twoDpts.resize( nb , std::vector<triple>(0) );
    triple tr;
    for(int i=0;i<n/2;++i){
        do{
            tr.a = 2*random_num()-1.0;
            tr.b = 2*random_num()-1.0;
            tr.S = tr.a*tr.a + tr.b*tr.b ;
        }
        while( tr.S > 1.0 );
        twoDpts[tr.S*nb].push_back( tr );
    }
    for (int bi = 0; bi < nb; ++bi)
        sort(twoDpts[bi].begin(),twoDpts[bi].end());
    
    int first_bucket = 0;
    while(!twoDpts[first_bucket ].size()) 
        ++first_bucket;
    int last_bucket = nb-1;
    while( !twoDpts[last_bucket].size() )
        --last_bucket;
    double t;
    double s = 1.0 / twoDpts[last_bucket].rbegin()->S ;
    s = s / ( 1.0 - s * twoDpts[first_bucket].begin()->a * twoDpts[first_bucket].begin()->a );
    double S_I,S_II;
    v[0] = twoDpts[first_bucket].begin()->b * sqrt( s ) ;
    S_II = twoDpts[first_bucket].begin()->S ;
    int i=1;
    int j=1;
    for (int bi = first_bucket; bi <=last_bucket; ++bi){
        for ( ;j<twoDpts[bi].size(); ++j){
            S_I  = S_II ;
            S_II = twoDpts[bi][j].S ;
            t = sqrt( (1.0 - S_I/S_II )*s ) ;
            v[i++] = twoDpts[bi][j].a * t ;
            v[i++] = twoDpts[bi][j].b * t ;
        }
        j=0;
        twoDpts[bi].resize(0);
    }
    return;
}

void rand_vect_bs(int n,double* a){ return n&1 ? rand_vect_bs_odd(n,a) : rand_vect_bs_even(n,a) ;}

/*****************************************inSitu*******************************************************************/

// as a prove of concept a in-situ version for even n:

void rand_vect_in_situ_even(int n,double* v){
    my_pair* twoDpts( reinterpret_cast<my_pair*>( v ) );
    for(int i=0;i<n/2;++i){
        do{
            twoDpts[i].S = 2*random_num()-1.0;
            twoDpts[i].b = 2*random_num()-1.0;
            // store squared radius and sign of a in S:
            twoDpts[i].S = (twoDpts[i].S>0) ?   twoDpts[i].S*twoDpts[i].S + twoDpts[i].b*twoDpts[i].b 
                                            : -(twoDpts[i].S*twoDpts[i].S + twoDpts[i].b*twoDpts[i].b) ;
        }
        while( fabs(twoDpts[i].S) > 1.0 );
    }
    std::sort( twoDpts , twoDpts + n/2 );
    double t;
    double s = 1.0 / fabs(twoDpts[n/2-1].S) ;
    double S_II(0.0);
    for(int i=0;i<n/2;++i){
        t  = S_II ;
        S_II =fabs( twoDpts[i].S) ;
        t = sqrt( (1.0 - t/S_II )*s ) ;
        twoDpts[i].S = twoDpts[i].S>0 ?  sqrt(  twoDpts[i].S - twoDpts[i].b*twoDpts[i].b ) * t
                                      : -sqrt( -twoDpts[i].S - twoDpts[i].b*twoDpts[i].b ) * t;
        twoDpts[i].b = twoDpts[i].b * t ;
    }
    return;
}

// We don't see an elegant way to solve the technical difficulties related to an insitu function for odd n
// mainly sorting (n+1)/2 pairs on an array of size n

/******************************************************************************************************************/

int main(void)
{
    double *v;
    w = new double[1000000];

    int D = 2;
	time_t t_start,t_curr;
    unsigned int N;
    int dt(10);
    
    while(D<1e7)
    {
        std::cout<<D;
/********************** Mueller **********************/
        time(&t_start);
        N = 0;
        do{
            rand_vect_mueller( D , w );
            ++N;
            time(&t_curr);
        }while( N<10 || difftime(t_curr,t_start)<dt );
        std::cout<<" "<<difftime(t_curr,t_start)*1.0/N;

/*********************** Basic ***********************/
        rand_vect_basic( D , w ); // call once to allocate memory
        time(&t_start);
        N = 0;
        do{
            rand_vect_basic( D , w );
            ++N;
            time(&t_curr);
        }while( N<10 || difftime(t_curr,t_start)<dt );
        std::cout<<" "<<difftime(t_curr,t_start)*1.0/N;

/********************* BucketSort ********************/
        rand_vect_bs( D , w ); // call once to allocate memory
        time(&t_start);
        N = 0;
        do{
            rand_vect_bs( D , w );
            ++N;
            time(&t_curr);
        }while( N<10 || difftime(t_curr,t_start)<dt );
        std::cout<<" "<<difftime(t_curr,t_start)*1.0/N;
        
/********************* InSitu (only even) ***********/
        if(!(D&1) )
        {
            time(&t_start);
            N = 0;
            do{
                rand_vect_in_situ_even( D , w );
                ++N;
                time(&t_curr);
            }while( N<10 || difftime(t_curr,t_start)<dt );
            std::cout<<" "<<difftime(t_curr,t_start)*1.0/N;
        }
        else
            std::cout<<" -1"; // to preserve output format

/****************************************************/
        std::cout<<std::endl;
/********************* increase D *******************/
        if(D&1)
            D = 2*floor(D*0.80901699435) ;
        else
            D = D+1;
    }
    return 0;
}

