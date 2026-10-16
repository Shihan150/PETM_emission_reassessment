#!/usr/bin/env python3
"""Run the current paired PETM alkalinity-feedback simulations."""
import argparse
import os
from pathlib import Path

def run_scenario(scenario):
    from iloscar_backend import parameters, init_start, model_run
    os.chdir(Path(__file__).resolve().parent)
    if scenario == 'conventional':
        # run test exerpiment
        params = parameters(pcycle = 0, LOADFLAG=1, cinpflag = 0, cinp = 6000, RUNTYPE=1, FTYS = 1)
        params['exp_name'] = 'petm_test_no_rw'
        params['printflag'] = 1

        params['pref'] = 500
        params['pcsi'] = 1000
        params['rrain'] = 6.7
        params['fsh'] = 4
        params['coeff_flag'] = 2

        params['thc'] = 25
        params['cac'] = 0.02
        params['mgc'] = 0.03
        params['svresults'] = 1

        params['fsi0'] = 6
        params['finc0'] = 12

        params['LOADFLAG'] = 1
        params['svstart'] = 0

        params['nsi'] = 0.2
        params['ncc'] = 0.4
        params['sclim'] = 4
        params['thcweak'] = 0.1
        params['atlantic_extra_stagnation'] = 0.2
        params['r_deep_atlantic'] = 0.2

        params['initfile'] = 'petm_steady_fsed.dat'
        params['svfile'] = 'petm_steady_fsed.dat'

        params['RUN_TYPE'] = 1
        params['cinpflag'] = 2

        params['cinpfile'] = 'emi_rw.dat'

        params['cinp'] = 5000
        params['tcspan'] = 5000

        params['t0'] = 0
        params['tfinal'] = 2e5
        init_start(params)    
        model_run()
    elif scenario == 'modified':
        # run test exerpiment
        params = parameters(pcycle = 0, LOADFLAG=1, cinpflag = 0, cinp = 6000, RUNTYPE=1, FTYS = 1)
        params['exp_name'] = 'petm_test_rw'
        params['printflag'] = 1

        params['pref'] = 500
        params['pcsi'] = 1000
        params['rrain'] = 6.7
        params['fsh'] = 4
        params['coeff_flag'] = 2

        params['thc'] = 25
        params['cac'] = 0.02
        params['mgc'] = 0.03
        params['svresults'] = 1

        params['fsi0'] = 6
        params['finc0'] = 12
        params['fsil2'] = 10

        params['LOADFLAG'] = 1
        params['svstart'] = 0

        params['nsi'] = 0.2
        params['ncc'] = 0.4
        params['sclim'] = 4
        params['thcweak'] = 0.1
        params['atlantic_extra_stagnation'] = 0.2
        params['r_deep_atlantic'] = 0.2

        params['initfile'] = 'petm_steady_fsed.dat'

        params['RUN_TYPE'] = 1
        params['cinpflag'] = 2

        params['cinpfile'] = 'emi_rw.dat'

        params['t0'] = 0
        params['tfinal'] = 2e5
        init_start(params)    
        model_run()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("scenario", choices=['conventional', 'modified'])
    args = parser.parse_args()
    run_scenario(args.scenario)

if __name__ == "__main__":
    main()
