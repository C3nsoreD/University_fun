"""
    Analytic calculations for FEM lab 2

"""
import math
import csv 

## Data 
IXX = 100970
IYY = 355416
A = 650 
E = 200000
DEN_ = 7.8e-9
L = 1140  

case1 = lambda i: 0.393 * ((2*i) - 1)**2
case2 = lambda i: 0.393 * ((2*i) + 1)**2
case3 = lambda i: 0.098 * ((4*i) + 1)**2
case4 = lambda i: 0.393 * ((2*i) - 1)**2

equations = {
    'case1': [0.560, 3.506, 9.825],
    'case2': [2.454, 7.951, 16.562],
    'case3': [3.560, 9.815, 19.257],
    'case4': [1, 1, 1],
}

numerial_results = {
    'case2': [27.179, 43.101, 101.98, 168.53, 192.62, 385.91, 463.72, 480.29],
    'case1': [49.696, 118.27, 152.36, 323.22, 377.10, 406.98, 617.58, 767.08],
    'case3': [150.92, 170.69, 358.37, 446.64, 460.61, 637.67, 876.03, 979.77],
    'case4': [180.07, 201.77, 413.37, 458.64, 506.73, 707.13, 928.28, 1059.1],
}


_term1 = 1 / L**2
_term2 = math.sqrt((E*IYY)/(DEN_*A))
c = _term1 * _term2

def _Fi_c(moi):
    """
        Calculates the f/c, that depends on the axis of bending moment (YZ, XZ) and 
        the bending mode 
    """

    F_C = {}
    F = {}
    for case, _f in equations.items():
        # f_c = [] 
        f = []
        f = [i * c for i in _f]
        
        # for i in range(1, 4):
        #     if i < 3:
        #         f_c.append(_f[i])
        
        # print("length of F_c \n", (f_c))

        # f =  [c * f_c[i] for i in range(len(f_c)) if i < 2]
        # _v_3 = f_c[2](3) * c

        # f =  [f_c[2](3) ]
        # print("f \n", _v_3)
        F_C[case], F[case] = _f, f
        # f_values[c] = _values

    return F_C, F


def relative_error(analytic_res, n_results):
    _error = {}
    for k, v in analytic_res.items():
        _error[k] = [
            abs((n_results[k][i] - v[i]) / v[i]) for i in range(len(v))
        ]
    # print(_error)
    return _error


def write_csv(results, filename=None):
    field_names = [name for name in results.keys()]
    values = [v for v in results.values()]
    # print(values, end='\n')
    with open(filename, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(field_names)

        rows = zip(*[v for v in results.values()])
        csv_writer.writerows(rows)


if __name__ == "__main__":
    moi = [IXX, IYY]
    
    f_C_xx, f_xx = _Fi_c(moi)
    results_yy, f_yy = _Fi_c(moi)
    
    _error_xx = relative_error(f_xx, numerial_results)
    _error_yy = relative_error(f_yy, numerial_results)
    
    print("F results in xx\n", f_xx)
    print()
    print("F/C xx \n", f_C_xx)
    print("error in yy \n",_error_xx)
    # print("error in xx, \n",_error_yy)
    
    # write_csv(results_xx, "results_xx.csv")
    # write_csv(results_yy, "results_yy.csv")

    write_csv(_error_xx, "error_xx.csv")
    # write_csv(_error_yy, "error_yy.csv")
    write_csv(f_xx, "f_xx.csv")
    # write_csv(f_yy, "f_yy.csv")
    
    print()
    # Tests 
    print(c)

    