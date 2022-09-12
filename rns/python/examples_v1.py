# Run the examples of v1:
#  fgrep "rns" ../rns.v1.1d/examples.test

import subprocess

rns = '../rns.v1.1d' # make sure rns is compiled there!

opt = {0:'-f eosC -t model -e 2e15 -r 0.59 -d 0',
       1:'-f eosC -t static -e 6e14 -l 2e15 -n 2 -p 2 -d 0',
       2:'-f eosC -t kepler -e 2e15 -d 0',
       3:'-f eosC -t gmass -e 1e15 -m 1.5 -d 0',
       4:'-f eosC -t rmass -e 1e15 -z 1.55 -d 0',
       5:'-f eosC -t omega -e 1e15 -o 0.5 -d 0',
       6:'-f eosC -t jmoment -e 1e15 -j 1.5 -d 0'}

for n,o in opt.items():
    subprocess.run(['{}/rns'.format(rns)] + o.replace("eos", "{}/eos/eos".format(rns)).split(' ')) 
 
