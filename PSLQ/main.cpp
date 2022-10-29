#include <iostream>
#include <iomanip>

#include "precision/fprecision.h";
#include "arpec_pslq/pslq1.h" 

using std::cout;
using std::endl;

 

int main(int argc, char** argv) {
 
    int n = 8;
    int nr_digits = 20;
    int n_eps = (nr_digits < 700 ? 10 : 20) - nr_digits;

    float_precision eps = pow(float_precision(10.0), float_precision(n_eps*10));
    cout << eps.precision() << " prec " << endl;
    matrix<float_precision> x(n);
    x(0) = "0.7522";
    x(1) = "0.14783";
    x(2) = "2.4556";
    x(3) = "-5.92046";
    x(4) = "0.035801";
    x(5) = "1.00333";
    x(6) = "0.002789";
    x(7) = "-0.00277728";

    matrix<int_precision> rel(n);

    cout << "nr_digits = " << nr_digits << endl;
    cout << "debug_level = " << debug_level << endl;
    cout << "n = " << n << endl;
 
    cout << "Starting PSLQ" << endl;
    /* Perform Level-1 PSLQ. */
    int result = pslq1(x, rel, eps);

    /* Output recovered relation. */
    if (result == RESULT_RELATION_FOUND) {
        cout << "Relation found:" << endl;
        cout << std::fixed << std::setprecision(0);
        for (int i = 0; i < n; i++) {
            cout << std::setw(3) << i;
            cout << std::setw(24) << rel(i) << endl;
        }
    }
    else {
        cout << "Precision exhausted." << endl;
    }

    return 0;
      
}

