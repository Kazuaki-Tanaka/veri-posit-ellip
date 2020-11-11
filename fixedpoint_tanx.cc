/*********************************************************************************************************
* This function calculates a rigorous enclosure of 
* the first positive zero of x = tan(x).
* The zero coincides with the first positive zero of the Bessel function with order 1.5.
* The computation relies on the function "allsol" packaged in the kv library.
*********************************************************************************************************
* witten 11/11/2020 K.Tanaka and T.Asai
*********************************************************************************************************/

#include <kv/allsol.hpp>
typedef kv::interval<double> itv;

struct Func {
    template <class T> T operator() (const T& x){
        return x - tan(x);
    }
};

typedef kv::interval<double> itv;

int main()
{
    std::cout.precision(17);
    itv PI;
    PI=kv::constants<itv>::pi();
    kv::allsol(Func(), itv(PI.lower(), 3.*PI.upper()/2.));
}
