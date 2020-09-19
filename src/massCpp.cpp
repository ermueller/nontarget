#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <numeric>
#include <vector>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R.h>

using namespace std;

struct dat{
       vector<double> mass;
       vector<double> retent;
       bool operator() (size_t i, size_t j) { return (mass[i]<mass[j]);}
       };

extern "C"
{
void mass(  double *mass,
            double *retent,
            int *a,
            double *masstol,
            double *massfrac,
            double *rttollow,
            double *rttolup,
            int *n_isos,
            double *isomat1,
            int *isomat3,
            double *maxmass1,
            int *isomat4,
            int *entry,
            int *ppm2,
            int *getit1,
            int *getit2,
            int *getit4,
            int *getit5,
            int *getit6
        )
    {
    
    cout << "Start" << endl;
    
    //Initialize variables
    size_t i=0, k=0, j=0, do_at=0, l=0, upcount=0, lowcount=0, rtup=0, rtlow=0, howmany=0;
    int entry2=*entry;
    double uptol, lowtol, thismasslow, thismasslow2, thismassup, thismassup2;
    

    //Generate index vectors:
    vector<double> datmass;
    vector<int> index;
    vector<int> index2;
    vector<int> index_test;
    
    // Read in the variables from R
    SEXP getit1bstore;
    PROTECT(getit1bstore = NEW_INTEGER(*a+1));
    int *getit1b;
    getit1b = INTEGER_POINTER(getit1bstore);
    for(k=0;k<(size_t)(*a+1);k++){*(getit1b+k) = 0;}

    SEXP getit2bstore;
    PROTECT(getit2bstore = NEW_INTEGER(*a+1));
    int *getit2b;
    getit2b = INTEGER_POINTER(getit2bstore);
    for(k=0;k<(size_t)(*a+1);k++){*(getit2b+k) = 0;}

    SEXP getit4bstore;
    PROTECT(getit4bstore = NEW_INTEGER(*a+1));
    int *getit4b;
    getit4b = INTEGER_POINTER(getit4bstore);
    for(k=0;k<(size_t)(*a+1);k++){*(getit4b+k) = 0;}

    SEXP getit5bstore;
    PROTECT(getit5bstore = NEW_INTEGER(*a+1));
    int *getit5b;
    getit5b = INTEGER_POINTER(getit5bstore);
    for(k=0;k<(size_t)(*a+1);k++){*(getit5b+k) = 0;}

    SEXP getit6bstore;
    PROTECT(getit6bstore = NEW_INTEGER(*a+1));
    int *getit6b;
    getit6b = INTEGER_POINTER(getit6bstore);
    for(k=0;k<(size_t)(*a+1);k++){*(getit6b+k) = 0;}

    // Read data into dat1
    dat dat1;
    
    // m/z
    for(i=0; i<(unsigned)*a; i++){
        dat1.mass.push_back(mass[i]);
    }
    
    // Retention time
    for(i=0; i<(unsigned)*a; i++){
        dat1.retent.push_back(retent[i]);
    }
    
    bool first_iter, rt_change;
    
    // Search loop:
    for(i = 0;i < (unsigned)*a; i++){
        
        // Find out where in the loop we are or if RT changed from the last peak
        first_iter = (i == 0);
        rt_change = (dat1.retent[i]!=dat1.retent[i-1]);
        
        // If this is the first iteration or the RT was changed:
        // Build new index vector with all peaks within the rt tolerance
        if(first_iter || rt_change){ 
            
            datmass = dat1.mass;
            
            //Get upper and lower tolerance
            uptol  = dat1.retent[i] + *rttolup;
            lowtol = dat1.retent[i] + *rttollow;
            
            //Get indices rtlow,rtup of the retetion times that are within the specified tolerance
            while(((rtlow + 1) < (unsigned)*a) && (dat1.retent[rtlow] < lowtol)) {
                rtlow++;
            }
            
            rtup = rtlow;
            
            while(((rtup + 1) < (unsigned)*a) && (dat1.retent[rtup + 1] <= uptol)){
                rtup++;
            }
            
            // Fill index with the numbers from rtlow to rtup
            index.clear();

            for(j = rtlow; j <= rtup; j++){
                    index.push_back(j);
            }
            index_test = index;
            

            // Sort the index by using a lambda (internal functions in C++, basically)
            sort(index.begin(), index.end(), [&](int i1, int i2) { return datmass[i1] < datmass[i2]; });
            
            // Alternative old way: Sort index vector by mass (see "dat.operator()")
            // This is slow because of the struct thing, I think. Lambdas are the fastest you could probably do here.
            
            // sort(index.begin(), index.end(), dat1);

            // Get the index of the mass in dat1 with the correct retention time
            do_at = 0;
            while((((do_at + 1) < index.size()) && (dat1.mass[index[do_at]] <= dat1.mass[i]))){
                do_at++;
            }

        }else{
            
            // Else, index is still valid from the last iteration
            
            // Get the index of the mass in dat1 with the correct retention time
            if(dat1.mass[i]<dat1.mass[i-1]){
                do_at=0;
            }
            
            while ((((do_at+1) < index.size()) && (dat1.mass[index[do_at]] <= dat1.mass[i]))) {do_at++;};
        }
        
        
        // How many peaks are in the specified RT window?
        howmany = index.size(); 
        
        // If there are peaks in the vector
        if(howmany>0){
            if(*ppm2==1){
                
                //Get mass tolerance (if supplied in ppm)
                thismasslow=(dat1.mass[i]-(dat1.mass[i]**masstol/1E6));
                thismasslow2=(dat1.mass[i]-(((dat1.mass[i]**masstol/1E6))**massfrac));
                thismassup=(dat1.mass[i]+(dat1.mass[i]**masstol/1E6));
                thismassup2=(dat1.mass[i]+(((dat1.mass[i]**masstol/1E6))**massfrac));
            }else{
                
                //Get mass tolerance (if supplied as absolute)
                thismasslow=(dat1.mass[i]-*masstol);
                thismasslow2=(dat1.mass[i]-(*masstol**massfrac));
                thismassup=(dat1.mass[i]+*masstol);
                thismassup2=(dat1.mass[i]+(*masstol**massfrac));
            }
            
            lowcount = do_at;
            
            // Go through all isotopes
            for(k = 0; k<(unsigned)*n_isos; k++){
                
                //Set the mass window
                while(((lowcount+1) < howmany) && (dat1.mass[index[lowcount]] < (thismasslow + isomat1[k]))){
                    lowcount++;
                }
                
                upcount = lowcount;
                
                while(((upcount+1) < howmany ) && (dat1.mass[index[upcount+1]] <= (thismassup + isomat1[k]))){
                    upcount++;
                }
                
                // Go through all peaks and add them to the output
                for(l=lowcount; l<=upcount; l++){
                    
                    // If the mass is within the isotope range, add it as a candidate
                    if( (dat1.mass[index[l]]<=(thismassup + isomat1[k])) && (dat1.mass[index[l]]>=(thismasslow + isomat1[k]))){
                        
                        // From?
                        if(*(getit2b+index[l])<(*entry+1)){ 
                            getit2[(index[l]**entry)+(*(getit2b+index[l]))]=(i+1);
                            *(getit2b+index[l]) = (*(getit2b+index[l])+1);
                        }
                        
                        // To?
                        if( *(getit4b+i)<(*entry+1) ){ 
                            getit4[i**entry + *(getit4b+i)]=(index[l]+1);
                            *(getit4b+i) = (*(getit4b+i)+1);
                        }
                        
                        // Which isotope?
                        if(*(getit1b+i)<(*entry+1)){ 
                            getit1[i**entry+*(getit1b+i)]=(k+1);
                            *(getit1b+i) = (*(getit1b+i)+1);
                        }
                        
                        // Which charge level?
                        if(*(getit6b+i)<(*entry+1)){ 
                            getit6[i**entry+*(getit6b+i)]=(isomat4[k]);
                            *(getit6b+i) = (*(getit6b+i)+1);
                        }
                        
                        // Large or small mass tolerance?
                        if(*(getit5b+i)<(*entry+1)){ 
                            if( (dat1.mass[index[l]]<=(thismassup2 + isomat1[k])) && (dat1.mass[index[l]]>=(thismasslow2 + isomat1[k]))){
                                getit5[i**entry+*(getit5b+i)]=1; // 1 = small
                                *(getit5b+i) = (*(getit5b+i)+1);
                            }else{
                                getit5[i**entry+*(getit5b+i)]=1; // 2 = large
                                *(getit5b+i) = (*(getit5b+i)+1);
                            };
                        }
                        
                        isomat3[k]=isomat3[k]+1;
                        
                    } // if within mass-window
                } // for l
            } // for k
        }  // if howmany>0
    } // for i
    //////////////////////////////////////////////////////////////////////////

    // check if entry has reached limit //////////////////////////////////////
    for(j=0;j<(size_t)(*a-1);j++){if(*(getit1b+j)>entry2){*entry=*(getit1b+j);};};
    for(j=0;j<(size_t)(*a-1);j++){if(*(getit2b+j)>entry2){*entry=*(getit2b+j);};};
    for(j=0;j<(size_t)(*a-1);j++){if(*(getit4b+j)>entry2){*entry=*(getit4b+j);};};
    for(j=0;j<(size_t)(*a-1);j++){if(*(getit5b+j)>entry2){*entry=*(getit5b+j);};};
    
    UNPROTECT(5);

    } // main
} // extern "C"









