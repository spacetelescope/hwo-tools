import sys
import pickle

import pytest
import numpy as np
from bokeh.models import ColumnDataSource

from main import do_recalculate_snr, load_initial
import catalog

from syotools.spectra.spec_defaults import pysyn_spectra_library


hwo = None
pivotwave = None
template_to_start_with = "Flat (AB)"

# to do: Allow nested values in the arguments (vary over values that are only valid/available
# given some other value)
def VaryMany(*args, **keys):
    inlist = []
    maxlen = 0
    # what is the longest list of parameters to vary over?
    for arg in args:
        #print(arg)
        if len(arg) > maxlen:
            maxlen = len(arg)
    # shuffle (and expand) all parameters to make random combinations
    expandargs = []
    for idx,arg in enumerate(args):
        #print(len(arg))
        # first, we can make repeats of the array as needed by resizing it to
        #  something larger
        arg = np.resize(arg,maxlen)
        #print(len(arg))
        # now we shuffle. We don't want to use randint() because every option must be
        #  picked at least once. This is an in-place operation.
        np.random.seed(idx)
        np.random.shuffle(arg)
        expandargs.append(arg)
    # now we have full lists of arguments and items, resized to the same length. Pair them up.
    for x in range(maxlen):
        calc = []
        for arg in expandargs:
            calc.append(arg[x])
        inlist.append(calc)

    return inlist

def check_relative_diff(actual, expected, rel_tol=0.1):
    """
    Simple function to check if two lists are approximately equal within a relative tolerance.
    Reports percentage differences for values that exceed the tolerance.

    Args:
        actual: List of actual values
        expected: List of expected values
        rel_tol: Relative tolerance (default: 0.1 or 10%)

    Returns:
        True if all values are within tolerance, False otherwise
    """
    if len(actual) != len(expected):
        print(f"Lists have different lengths: actual={len(actual)}, expected={len(expected)}")
        return False

    all_within_tolerance = True
    differences = []

    for i, (a, e) in enumerate(zip(actual, expected)):
        for j, (act, exp) in enumerate(zip(actual[i], expected[i])):
            if exp == 0:
                # Can't calculate relative difference if expected is zero
                if act != 0:
                    all_within_tolerance = False
                    differences.append((i, a, e, "inf"))
            else:
                rel_diff = abs(act - exp) / abs(exp)
                pct_diff = 100 * rel_diff

                if rel_diff > rel_tol:
                    all_within_tolerance = False
                    differences.append((i, j, act, exp, pct_diff))

    if not all_within_tolerance:
        print("\nValues exceeding relative tolerance:")
        for i, j, act, exp, pct in differences:
            if pct == "inf":
                print(f"  Index {i},{j}: actual={act}, expected={exp}, difference=infinite (division by zero)")
            else:
                print(f"  Index {i},{j}: actual={act:.6g}, expected={exp:.6g}, difference={pct:.2f}%")

    return all_within_tolerance


load_initial()

class item():
    def __init__(self, value=30):
        self.value = value

target_planet, target_star = catalog.load_catalog()

snr = item(10)
exptime = item(30.0)
diameter = item(7.0)
star = item("G2V star")
planet = item("Earth")
delta_mag = item(15)
distance = item(10)

snrs = np.linspace(0.1, 100.0, 30)
diameters = np.linspace(5, 15, 10)
stars = list(target_star)
planets = list(target_planet)
exptimes = np.linspace(1, 1000.0, 30)
distances = np.linspace(1.4, 100.0, 30)
delta_mags = np.linspace(10,30,20)

test_setups = VaryMany(snrs, diameters, stars, planets, exptimes, distances, delta_mags)

def do_comparisons(inputs):
    global snr
    global exptime
    global diameter
    global star
    global planet
    global delta_mag
    global distance

    snr = inputs[0]
    diameter = inputs[1]
    star = inputs[2]
    planet = inputs[3]
    exptime = inputs[4]
    distance = inputs[5]
    delta_mag = inputs[6]
    new_values = ColumnDataSource(data={"new_snr": [snr], "new_exp": [exptime], "new_diameter": [diameter], "new_star": [star], "new_distance": [distance], 
                  "new_planet": [planet], "delta_mag": [delta_mag], "observation": [True], "observatory": [True], "scene": [True]})
    obsdata = do_recalculate_snr(new_values)

    print(obsdata.data)

    return obsdata.data["exptime"], obsdata.data["FpFs"], obsdata.data["snr"]

def create_comparisons(reset):
    test_results = []
    for setup in test_setups:
        test_results.append((setup, do_comparisons(setup)))
    if reset:
        with open("test_coronspec.pickle", "wb") as picklefile:
            pickle.dump(test_results, picklefile)

'''
LOAD IT
'''
try:    
    with open("test_coronspec.pickle", "rb") as picklefile:
        test_results = pickle.load(picklefile)
except FileNotFoundError:
    create_comparisons(True)
    with open("test_coronspec.pickle", "rb") as picklefile:
        test_results = pickle.load(picklefile)

@pytest.mark.parametrize("testset", test_results)
def test_camera(testset):
    
    inputs = testset[0]
    expected = testset[1]
    print(inputs)
    exptime_out, FpFs, snr_out = do_comparisons(inputs)

    assert check_relative_diff(expected, (exptime_out, FpFs, snr_out))


if __name__ == "__main__":
    reset=False
    if (len(sys.argv) > 1) and (sys.argv[1] == "reset"):
        reset = True
    create_comparisons(reset)
