import os
import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR   = os.path.dirname(SCRIPT_DIR)
VDIR       = os.path.join(ROOT_DIR, 'validation')


def convert_complex(txt):
    txt = txt.strip("()").replace(",", "+").replace("+-", "-") + "j"
    return complex(txt)


def calc_pct_diff(ref_list, comp_list, scalefactor):
    x1 = np.array(ref_list)
    x2 = np.array(comp_list)
    eps = 1.e-15
    x1[x1 < eps] = 0.0
    x2[x2 < eps] = 0.0
    pctdiff_arr = 100 * scalefactor * (abs(x2 - x1)) / x1
    np.nan_to_num(pctdiff_arr, nan=0.0)
    return pctdiff_arr


if __name__ == '__main__':
    with open(os.path.join(VDIR, 'TEST1.dat')) as f1:
        test1lines = f1.readlines()

    theta   = []
    res_ana = []
    res_m1  = []
    res_m2  = []

    for line in test1lines:
        thestr = line.split()
        theta.append(thestr[0])
        res_ana.append(convert_complex(thestr[1]))
        res_m1.append(convert_complex(thestr[2]))
        res_m2.append(convert_complex(thestr[3]))

    real_pctdiff_m1 = calc_pct_diff(
        [e.real for e in res_ana], [e.real for e in res_m1], 1)
    imag_pctdiff_m1 = calc_pct_diff(
        [e.imag for e in res_ana], [e.imag for e in res_m1], 1)

    print(imag_pctdiff_m1)

    fig, ax = plt.subplots()
    ax.plot(theta, real_pctdiff_m1, 'blue')
    ax.plot(theta, imag_pctdiff_m1, 'red')
    plt.show()
