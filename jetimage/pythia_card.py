#!/usera/ss2165/anaconda2/bin/python
"""Usage:
    pythia_card.py <process> <file_name> <n_events> [-s=SEED] [--boson-mass=BOSON_MASS] [--pt-hat-min=PTHMIN] [--pt-hat-max=PTHMAX]

Arguments:
    <process>    Process to generate card for
    <file_name>  Card to generate
    <n_events>   Number of events

Options:
    -s=SEED                   Random seed
    --boson-mass=BOSON_MASS   Specify boson mass GeV [default: 800]
    --pt-hat-min=PTHMIN       pthatmin GeV [default: 100]
    --pt-hat-max=PTHMAX       pthatmax GeV [default: 500]
"""
import time
from os import getpid, path, remove
from docopt import docopt
import subprocess

arguments = docopt(__doc__, help=True)
print(arguments)
process = arguments['<process>']
fname = arguments['<file_name>']
seed = arguments['-s']

if seed is None:
    tseed = time.clock()
    seed = abs(((tseed * 181) * ((getpid() - 83) * 359)) % 104729)

home = path.expanduser('~')
template_file = path.abspath(path.join(home, 'pt3proj', 'jetimage', 'template.cmnd'))
with open(template_file, "r") as f:
    lines = f.readlines()

r_index = lines.index('Random:setSeed = on\n')
lines[r_index+1] = 'Random:seed = {}\n'.format(seed)
lines[0] = 'Main:numberOfEvents = {}\n'.format(int(arguments['<n_events>']))
z_index = lines.index('!1MASS\n')
wzl_index = lines.index('!2MASS\n')
wzh_index = lines.index('!3MASS\n')
qcd_index = lines.index('HardQCD:all = on\n')

if process=="wprime":
    lines[wzl_index] = '34:m0={}\n'.format(arguments['--boson-mass'])
    lines = lines[:z_index] + lines[wzl_index:wzh_index]
elif process=="qcd":
    ind = lines.index('ptHatMin.str()\n')
    lines[ind] = 'PhaseSpace:pTHatMin  ={}\n'.format(arguments['--pt-hat-min'])
    lines[ind+1] = 'PhaseSpace:pTHatMin  ={}\n'.format(arguments['--pt-hat-max'])
    lines = lines[:z_index] + lines[qcd_index:]
else:
    raise ValueError("process = wprime/qcd")

pythia_card = path.abspath('/usera/ss2165/tmp/tmppythcard.cmnd')
with open(pythia_card, "w") as f:
    f.writelines(lines)

delphes_card = path.abspath(path.join(home,'pt3proj/configs/ATLAS_me.tcl'))
# remove(fname)
shell_command = ['DelphesPythia8', delphes_card, pythia_card, fname]
# proc = subprocess.Popen(shell_command, stdout=subprocess.PIPE)
t0 = time.time()
print(t0)
subprocess.call(shell_command)
print("Runtime for {}".format(arguments['<n_events>']))
print(time.time()-t0)
remove(pythia_card)