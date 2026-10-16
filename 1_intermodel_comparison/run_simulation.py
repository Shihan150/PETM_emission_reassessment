#!/usr/bin/env python3
"""Run finite-sediment comparison scenarios from the current experiment notebook."""
import argparse
import os
from pathlib import Path

def run_scenario(scenario):
    from iloscar_backend import parameters, init_start, model_run
    os.chdir(Path(__file__).resolve().parent)
    if scenario == 'gutjahr2017':
        # Gutjahr et al., 2017
        params = parameters(pcycle = 0, LOADFLAG=1, cinpflag = 0, cinp = 6000, RUNTYPE=1, FTYS = 1)
        params['exp_name'] = 'gutjahr2017_fsed'
        params['printflag'] = 1
        params['LOADFLAG'] = 0
        params['pref'] = 863
        params['pcsi'] = 863
        params['rrain'] = 6.7
        params['fsh'] = 1
        params['d13cin'] = 2
        params['d13cvc'] = -1.5
        params['thc'] = 25
        params['cac'] = 0.02
        params['mgc'] = 0.03
        params['svresults'] = 1

        params['fsi0'] = 7.5
        params['finc0'] = 17.5

        params['exp_rrain'] = 0
        params['initfile'] = 'gutjahr2017_fsed.dat'
        params['LOADFLAG'] = 1
        params['target'] = 'pH_d13c'
        params['target_file'] = 'data/gutjahr_ph.csv'
        params['nsi'] = 0.4
        params['ncc'] = 0.2

        params['target_file2'] = 'data/gutjahr_d13c.csv'
        params['tfinal'] = 2e5
        params['RUN_TYPE'] = 2
        params['toms_low'] = -0.1 
        params['toms_high'] = 1
        params['svresults'] = 1
        init_start(params)    
        model_run()
    elif scenario == 'zeebe2009':

        params = parameters(pcycle = 0, LOADFLAG=1, cinpflag = 0, cinp = 6000, RUNTYPE=1, FTYS = 1)
        params['exp_name'] = 'zeebe2009_fsed'
        params['printflag'] = 1

        params['pref'] = 834
        params['pcsi'] = 834
        params['rrain'] = 6.7
        params['fsh'] = 1
        params['d13cin'] = 2
        params['d13cvc'] = -1.5
        params['thc'] = 25
        params['cac'] = 0.02
        params['mgc'] = 0.03
        params['svresults'] = 1

        params['fsi0'] = 7.5
        params['finc0'] = 17.5

        params['initfile'] = 'gutjahr2017_fsed.dat'
        params['LOADFLAG'] = 1
        params['svstart'] = 0

        params['nsi'] = 0.4
        params['ncc'] = 0.2

        params['RUN_TYPE'] = 1
        params['cinpflag'] = 2

        params['cinpfile'] = 'Zeebe2009emi.dat'
        params['t0'] = 0
        params['tfinal'] = 2e5
        init_start(params)    
        model_run()
    elif scenario == 'this-study':
        # Gutjahr et al., 2017
        params = parameters(pcycle = 0, LOADFLAG=1, cinpflag = 0, cinp = 6000, RUNTYPE=1, FTYS = 1)
        params['exp_name'] = 'inv_rw'
        params['printflag'] = 1

        params['pref'] = 500
        params['pcsi'] = 1000
        params['rrain'] = 6.7
        params['fsh'] = 4.5
        params['coeff_flag'] = 2

        params['thc'] = 25
        params['cac'] = 0.02
        params['mgc'] = 0.03
        params['svresults'] = 1

        params['fsi0'] = 6
        params['finc0'] = 12

        params['fsil2'] = 5

        params['LOADFLAG'] = 1
        params['svstart'] = 0

        params['nsi'] = 0.2
        params['ncc'] = 0.4
        params['sclim'] = 4

        params['initfile'] = 'petm_steady_fsed.dat'

        params['exp_rrain'] = 0

        params['LOADFLAG'] = 1
        params['target'] = 'pH_d13c'
        params['target_file'] = 'data/gutjahr_ph.csv'

        params['target_file2'] = 'data/gutjahr_d13c.csv'
        params['tfinal'] = 2e5
        params['RUN_TYPE'] = 2
        params['toms_low'] = -0.1 
        params['toms_high'] = 1
        params['svresults'] = 1
        init_start(params)    
        model_run()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("scenario", choices=['gutjahr2017', 'zeebe2009', 'this-study'])
    args = parser.parse_args()
    run_scenario(args.scenario)

if __name__ == "__main__":
    main()
