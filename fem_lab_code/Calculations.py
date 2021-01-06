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
DEN_ = 7.810e-9
L = 1140  

case1 = lambda i: 0.393 * math.pow((2*i - 1), 2)
case2 = lambda i: 0.098 * math.pow((4*i + 1), 2)
case3 = lambda i: 0.393 * math.pow((2*i + 1), 2)
case4 = lambda i: 0.393 * math.pow((2*i - 1), 2)

equations = {
    'case1': case1,
    'case2': case2,
    'case3': case3,
    'case4': case4,
}

numerial_results = {
    'case1': [27.179, 43.101, 101.98, 168.53, 192.62, 385.91, 463.72, 480.29],
    'case2': [49.696, 118.27, 152.36, 323.22, 377.10, 406.98, 617.58, 767.08],
    'case3': [150.92, 170.69, 358.37, 446.64, 460.61, 637.67, 876.03, 979.77],
    'case4': [180.07, 201.77, 413.37, 458.64, 506.73, 707.13, 928.28, 1059.1],
}


def c_analytic(moi):
    """ 
        Analytic equation thats dependant on moment of inertia
    
    """
    _term1 = math.sqrt((E*moi)/(DEN_*A))
    _term2 = 1 / L**2
    return (_term2 * _term1)


def _Fi_c(moi):
    """
        Calculates the f/c, that depends on the axis of bending moment (YZ, XZ) and 
        the bending mode 
    """

    results = {}
    f_values = {}
    for c, f in equations.items():
        # results[c] = [f(i) for i in range(1, 9)]
        _fi = []
        f_c = []
        for i in range(1, 9):
            if i in [1, 4, 7]:
                _fi.append(f(i) * c_analytic(moi[0]))
                f_c.append(f(i))
            else:
                _fi.append(f(i) * c_analytic(moi[1]))
                f_c.append(f(i))

        results[c], f_values[c] = f_c, _fi
        # f_values[c] = _values

    return results, f_values


def relative_error(analytic_res, n_results):
    _error = {}
    for k, v in analytic_res.items():
        _error[k] = [
            abs((n_results[k][i] - v[i]) / v[i]) for i in range(len(v))
        ]
    # print(_error)
    return _error


def write_csv(results, filename="results.csv"):
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
    
    results_xx, f_xx = _Fi_c(moi)

    results_yy, f_yy = _Fi_c(moi)
    
    _error_xx = relative_error(f_xx, numerial_results)
    _error_yy = relative_error(f_yy, numerial_results)
    
    write_csv(results_xx, "results_xx.csv")
    write_csv(results_yy, "results_yy.csv")

    write_csv(_error_xx, "error_xx.csv")
    write_csv(_error_yy, "error_yy.csv")
    write_csv(f_xx, "f_xx.csv")
    write_csv(f_yy, "f_yy.csv")
    
    # Tests 
    print(c_analytic(moi[0]))

    