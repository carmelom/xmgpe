#!/usr/bin/python2
import os, shutil
from math import pi, sqrt, log
from lxml import etree
import argparse

# Physics
hbar = 1.05457e-34
kB = 1.38065e-23
# Sodium constants
mass = 3.81924e-26
sclength = 2.75172e-09

N = 100
nu_z = 100
AR = 10

dx = 0.33
B = 2

parser = argparse.ArgumentParser()
parser.add_argument("script", type=str,
                    help="real time simulation runtime")

parser.add_argument("-u", "--u", type=float, default=1.0,
                    help="pin strength (in mu)")
parser.add_argument("-v", "--v", type=float, default=1.0,
                    help="barrier speed (in unit of sound speed)")
parser.add_argument("-w", "--w", type=float, default=0.05,
                    help="barrier starting point")
parser.add_argument("-x", "--x1", type=float, default=0.5,
                    help="barrier starting point")
aux_arguments = ['u', 'v', 'w', 'x1']

######### static arguments

parser.add_argument("-r", "--runtime", type=float, default=10,
                    help="real time simulation runtime (in 1/omega_z)")
parser.add_argument("-s", "--sample-timestep", nargs='+', type=float, default=[0.1, 0.1, 0.05],
                    help="timestep for wavefunction sample output")

parser.add_argument("-D", "--dir", type=str,
                    help="override directory for output storage")
parser.add_argument("-n", "--no-compile", action="store_true",
                    help="Only rewrite .xmds source file")
args = parser.parse_args()




### Set globals
omega_z = 2*pi * nu_z
a_z = sqrt(hbar/(mass*omega_z))

g1 = (4*pi*hbar**2 * sclength/mass)/(hbar*omega_z*a_z**3)
mu1 = 0.5 * (15*AR**2 * N*sclength/a_z)**(2./5)
R_z = sqrt(2*mu1)

const = ['Nparticles = %g'%N,
         'omega_z = 2*M_PI * %g'%nu_z,
         'AR = %g'%AR,
        ]

calc = ['Uint = %g'%g1,
        'mu1 = %g'%mu1,
        'R_z = %g'%R_z, # extra params added here
        ]

sep = '\n\t'
cdata_str = sep.join(['const real %s;'%s for s in const] + [''] +\
                      ['real %s;'%s for s in calc] + [sep+'//arguments'] +\
                      ['real %s = %g;'%(s, getattr(args, s)) for s in aux_arguments]
                     )

### Set geometry
names = ['x', 'y', 'z']

ntilde_z = log(4*mu1*B/dx, 2)
ntilde_r = log(4*mu1*B/dx/AR, 2)

n_z = int(round(ntilde_z))
n_r = int(round(ntilde_r))

stepx = dx * (1/R_z) #dx * xi1
L_z = 2**n_z * stepx
L_r = 2**n_r * stepx

LL = [L_r, L_r, L_z]
nn = [2**n_r, 2**n_r, 2**n_z]
B_actual = L_z/(2*sqrt(2*mu1))

### Estimate simulation time
from operator import mul
volume = reduce(mul, LL, 1)
grid_size = reduce(mul, nn, 1)

time_per_point = 3e-7
runtime = args.runtime + 1
sample_ts = args.sample_timestep[0]

ts = 0.10 * stepx**2 #Heuristic based on Von Neumann criterion

total = time_per_point *runtime * grid_size / ts
# print('ts: %.3e'%ts)
m, s = divmod(total, 60)
h, m = divmod(m, 60)
# print("total: %d:%02d:%02d (%.2f s)" % (h, m, s, total))


### physical_report
gamma = mu1/AR
a_r = 1e6*a_z/AR
Rz = sqrt(2*mu1)*a_z * 1e6 #um
phys_mu = mu1*hbar*omega_z #J
phys_xi = 1e9*a_z/sqrt(2*mu1) #nm


physical_report = """BEC simulation for:
N particles = %g
Trap freq: 2*pi* (%g x %g x %g) Hz (AR = %g)

---- 3d inequality ----
mu3d > AR > 1
(%g, %g, %g)

gamma3d = %.2f
(stable vortices at gamma > 2.65)
(stable rings at gamma > 4 e qualcosa)

---- phisical parameters ----
R_tf: (%.2f x %.2f x %.2f) um
mu1 = %.2f nK
xi = %.2f nm

---- geometry ----
Box size: %.2f R_tf
%d x %d x %d grid points
resolution: %.2f nm (%.2f xi)

runtime: %.1f (1/omega_z)

---- expected simulation time ----
avg timestep: %e
total: %d:%02d:%02d (%.2f s)
------------------------\n
"""%  (N,
       AR*nu_z, AR*nu_z, nu_z, AR,
       mu1, AR, 1,
       gamma,
       Rz/AR, Rz/AR, Rz,
       1e9*phys_mu/kB,
       phys_xi,
       B_actual,
       nn[0], nn[1], nn[2],
       dx*phys_xi, dx,
       args.runtime,
       ts,
       h, m, s, total)

### Prepare folder for run and storage
if not args.no_compile:
    ### Write xml
    parser = etree.XMLParser(strip_cdata=False, remove_blank_text=True)
    with open(args.script, "r") as source:
        et = etree.parse(source, parser=parser)

    def indent(elem, level=0):
        i = "\n" + level*"  "
        i2 = "\n\n" if level < 2 else "\n"
        i2 += level*"  "
        if len(elem):
            if not elem.text or not elem.text.strip():
                elem.text = i + "  "
            if not elem.tail or not elem.tail.strip():
                elem.tail = i2
            for elem in elem:
                indent(elem, level+1)
            if not elem.tail or not elem.tail.strip():
                elem.tail = i2
        else:
            if level and (not elem.tail or not elem.tail.strip()):
                elem.tail = i2

    scriptname = os.path.splitext(args.script)[0]
    name = et.find('name')
    name.text = scriptname

    glob = et.find('features').find('globals')
    glob.text = etree.CDATA(sep+cdata_str+sep)
    # print(etree.tostring(glob, encoding='unicode', pretty_print=True))

    geom = et.find('geometry')
    td = etree.Element('transverse_dimensions')


    for j, name in enumerate(names):
        dim = etree.Element('dimension')
        dim.attrib['name'] = name
        dim.attrib['lattice'] = str(nn[j])
        dim.attrib['domain'] = '(-{0:.3f}, {0:.3f})'.format(LL[j]/2)
        td.append(dim)
    geom.replace(geom.find('transverse_dimensions'), td)
    # print(etree.tostring(geom, encoding='unicode', pretty_print=True))

    sequence = et.find('sequence')
    int_real = list(sequence.iterfind('integrate'))[-1]
    samples = int_real.find('samples')

    nsamples = [str(int(args.runtime / ts)) for ts in args.sample_timestep]
    int_real.attrib['interval'] = str(args.runtime)
    samples.text = ' '.join(nsamples)
    print samples.text

    # dest = 'test_configure.xmds'
    dest = args.script
    print "Rewrite on %s...\n"%dest
    indent(et.getroot(), level=0)

    et.write(dest, )

    print 'recompiling simulation'
    os.system('xmds2 {:s}'.format(dest))

    dirname = '.'.join(['{:s}_{:.2g}'.format(name, getattr(args, name)) for name in aux_arguments ])
    directory = dirname if args.dir is None else args.dir
    print 'configuring %s directory'%directory

    if not os.path.exists(directory):
       os.makedirs(directory)
    shutil.copy2(scriptname, directory)
    with open(os.path.join(directory, 'phys_report.txt'), 'w') as f:
        f.write(physical_report)

# else:
print physical_report
