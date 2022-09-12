#!/usr/bin/python3

import os, sys, shutil
from subprocess import Popen, PIPE
import numpy as np

def build_command(exe,params):
    cmd = exe
    for k,v in params.items():
        cmd += ' {} {}'.format(k,v)
    return cmd

def run_command(cmd, workdir, stdout, verbose):
  proc = Popen(cmd, cwd=workdir, stdout=PIPE, universal_newlines=True)
  if verbose:
    out = []
    while True:
      line = proc.stdout.readline()
      if not line:
        break
      else:
        if line is not None:
          sys.stdout.write(line)
          out.append(line)
    out = "".join(out)
  else:
    out, err = proc.communicate()
  open(os.path.join(workdir, stdout), "w").write(out)
  return out

# Make sure these executable are in place!
EXE = {
    'v1': '../rns.v1.1d/rns',
    'v2': '../rotstar/rotstar.x',
    'kepler': '../rotstar/kepler.x' 
}

class RNS:
    
    def __init__(self, version = 'v1', params = {}, eosfolder = '../rns.v1.1d/eos/', workdir="."):

        if version not in EXE.keys():
            raise ValueError('unknown code version {}'.format(version))
        
        self.exe = EXE[version]
        self.params = params
        self.workdir = workdir
        self.eosfolder = eosfolder

    def remove_none(self,d):
        return {k: v for k, v in d.items() if v is not None}

    def set_params_v1(self,
                      eos = 'tab',
                      N = 1,
                      eosfile = 'eosC',
                      ec = 2e15,
                      ecl = None,
                      a = 0.59,
                      M = None,
                      Mb = None,
                      Omg = None,
                      J = None,
                      task = 'model',
                      pr = 1,
                      relfact = None,
                      prconv = 0,
                      acc = 1e-5,
                      acc_q = 1e-4,
                      cf = 1.,
                      nomodels = None):
        """
        Set the options for rns v1
        (and reset the executable to this version)
        
        Note the tasks of v1 are:
        
        model : requires -e and -r
        gmass : requires -e and -m
        rmass : requires -e and -z
        omega : requires -e and -o
        jmoment : requires -e and -j
        static : requires -e
        kepler : requires -e
        eplot  : requires -e
        """
        params = {'-q': eos, # EOS type 'tab' : tabulated, 'poly' : analytic polytropic 
                  '-N': N, # polytropic index (P=K*e^(1+1/N))
                  '-f': '{}/{}'.format(self.eosfolder,eosfile), # EOS file 
                  '-e': ec, # central energy density in gr/cm^3
                  '-l': ecl,# central energy density of last model
                  '-r': a, #axes ratio
                  '-m': M, # mass in solar masses
                  '-z': Mb, # rest mass in solar masses
                  '-o': Omg, #angular velocity in 10^4 s^-1
                  '-j': J, # angular momentum in G*M_SUN^2/C
                  '-t': task, # task to be performed:
                  '-p': pr, # printing option  
                  '-c': cf, # relaxation factor
                  '-d': prconv, # if 0, do not monitor convergence 
                  '-a': acc, # accuracy in convergence
                  '-b': acc_q, # accuracy in fixing M, Omega, etc.
                  '-n': nomodels, # number of models, if sequence 
                }
        if eos == 'tab':
            params['-N'] = None
        elif eos == 'poly':
            params['-f'] = None
        else:
            raise ValueError('unknown type of EOS')        
        if task == 'model':
            params['-m'] = None
            params['-z'] = None
            params['-o'] = None
            params['-j'] = None
        elif task == 'gmass':
            params['-r'] = None
            params['-z'] = None
            params['-o'] = None
            params['-j'] = None
        elif task == 'rmass':
            params['-r'] = None
            params['-m'] = None
            params['-o'] = None
            params['-j'] = None
        elif task == 'omega':
            params['-r'] = None
            params['-m'] = None
            params['-z'] = None
            params['-j'] = None
        elif task == 'jmoment':
            params['-r'] = None
            params['-m'] = None
            params['-z'] = None
            params['-o'] = None
        elif task == 'static':
            params['-r'] = None
            params['-m'] = None
            params['-z'] = None
            params['-o'] = None
            params['-j'] = None
        elif task == 'kepler':
            params['-r'] = None
            params['-m'] = None
            params['-z'] = None
            params['-o'] = None
            params['-j'] = None
        elif task == 'eplot':
            params['-r'] = None
            params['-m'] = None
            params['-z'] = None
            params['-o'] = None
            params['-j'] = None
        else:
            raise ValueError('unknown taskS')
        self.params = self.remove_none(params)
        self.exe = EXE['v1'] 
        return

    def set_params_v2(self,
                      eos = 'poly',
                      N = 1,
                      K = 100,
                      eosfile = None,
                      ec = 1.1870e-03,
                      a = 8.5200e-01,
                      save = None, 
                      fname = None):
        """
        Set the options for rns v2/rotstar
        (and reset the executable to this version)
        """
        params = {'-q': eos, # EOS type 'tab' : tabulated, 'poly' : analytic polytropic 
                  '-N': N, # polytropic index (P=K*e^(1+1/N))
                  '-K': K, #polytropic constant (P=K*e^(1+1/N)) {100}
                  '-f': '{}/{}'.format(self.eosfolder,eosfile), # EOS file 
                  '-e': ec, # central energy density [CGS] tab EOS (gr/cm^3) [c=G=Msun=1] poly EOS 
                  '-r': a, #axes ratio
                  '-s': save, # save? 'yes'/'no'
                  '-O': fname, # output filename
                  }
        if eos == 'tab':
            params['-N'] = None
            params['-K'] = None
        elif eos == 'poly':
            params['-f'] = None
        else:
            raise ValueError('unknown type of EOS')        
        self.params = self.remove_none(params)
        self.exe = EXE['v2'] 
        return 

    def set_params_kepler(self,
                          eos = 'poly',
                          N = 1,
                          eosfile = 'eosC',
                          ec = 1e15):
        """
        Set the options for rns v2/kepler
        (and reset the executable to this version)
        """
        params = {'-q': eos, # EOS type 'tab' : tabulated, 'poly' : analytic polytropic 
                  '-N': N, # polytropic index (P=K*e^(1+1/N))
                  '-f': '{}/{}'.format(self.eosfolder,eosfile), # EOS file 
                  '-e': ec, # central energy density [CGS] tab EOS (gr/cm^3) [c=G=Msun=1] poly EOS 
                  }
        if eos == 'tab':
            params['-N'] = None
        elif eos == 'poly':
            params['-f'] = None
        else:
            raise ValueError('unknown type of EOS')
        self.params = self.remove_none(params)
        self.exe = EXE['kepler'] 
        return

    def run(self,stdout='rns.out',verbose=True):
        if not self.exe:
            raise ValueError('no executable')
        if not self.params:
            raise ValueError('empty parameters')
        cmd = build_command(self.exe,self.params)
        return run_command(cmd.split(), self.workdir, stdout, verbose)

    def parse_output_v1_1(self,out):
        """ Parse output created with option pr = 2
        """
        out = [s for s in out.splitlines() if s] 
        out = [s.strip() for s in out[2:-1]]
        return dict([("".join(s.split()[1:]),float(s.split()[0])) for s in out])

    def parse_output_v1_2(self,out,headerlines=5, keysline=3):
        """ Parse output created with option pr = 2
        FIXME: rns sometimes output values for more than 20 fields ...
        """
        out = [s for s in out.splitlines() if s] 
        keys = out[keysline].split() # usually len(keys) = 20
        out = list(np.concatenate([s.split() for s in out[headerlines:]]))
        out = [np.nan  if s == '---' else s for s in out]
        out = np.array([float(s) for s in out]).reshape((len(out)//len(keys),len(keys)))
        return dict([(k, out[:,i]) for i,k in enumerate(keys)])

    def parse_output_v1(self,out):
        if self.params['-p'] == 1:
            return self.parse_output_v1_1(out)
        elif self.params['-p'] == 2:
            return self.parse_output_v1_2(out)
        else:
            raise ValueError('unknown output')
            
    def parse_output_v2(self,out,headerlines=1,keysline=0):
        """ Parse rotstar output
        """
        if self.params['-q'] == 'tab':
            headerlines = 2
        out = [s for s in out.splitlines() if s] 
        keys = out[keysline].split()[1:] # rm #
        vals = [float(s) for s in out[headerlines].split()]
        return dict(zip(keys, vals))
    
    
if __name__ == "__main__":

    
    ns = RNS()
    ns.set_params_v1()
    ##ns.set_params_v1(pr=2) # this gives an error parsing 21 values for 20 fields...
    data = ns.parse_output_v1(ns.run())
    print(data)

    #ns.set_params_v2()
    ##ns.set_params_v2(eos = 'tab', eosfile = 'eosC', ec = 2e15, a = 5.68311e-01)
    #data = ns.parse_output_v2(ns.run())
    #print(data)

    
