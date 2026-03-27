from .fileReading import *


def itpwriter(func_name):
    """
    writes gromacs itp file
    :param func_name: takes the name of the built system
    :return: gromacs itp file
    """
    sections = make_dic(f'{func_name}.top')

    atomtypeDict = {}
    for i, l in enumerate(sections['[ atoms ]']):
        a = l.split()[:2]
        atomtypeDict[a[0]] = a[1]

    with open(f'{func_name}.itp', 'w+') as g:
        g.writelines('\n')
        g.writelines('[ moleculetype ]\n')
        g.writelines('; Name       nrexcl\n')
        for i, l in enumerate(sections['[ moleculetype ]']):
            s = l.split()
            g.writelines('{:>6s}{:>10s}\n'.format(s[0], s[1]))
        g.writelines('\n')
        g.writelines('\n')
        g.writelines('[ atoms ]\n')
        g.writelines(
            ';   nr       type  resnr residue  atom   cgnr    charge       mass  typeB    chargeB      massB\n')
        g.writelines('; residue    1 SUR rtp SUR q -0.0\n')
        for i, l in enumerate(sections['[ atoms ]']):
            s = l.split()
            g.writelines(
                '{:>6s}{:>14s}{:>7s}{:>7s}{:>7s}{:>7s}{:>14s}{:>14s}{:>6s}{:>5s}{:>12s}\n'.format(s[0], s[1], s[2],
                                                                                                  s[3],
                                                                                                  s[4], s[5], s[6],
                                                                                                  s[7],
                                                                                                  s[8], s[9], s[10]))
        g.writelines('\n')

        g.writelines('[ bonds ]\n')
        g.writelines(';    ai     aj funct         c0         c1         c2         c3\n')
        for i, l in enumerate(sections['[ bonds ]']):
            s = l.split()
            g.writelines('{:>6s}{:>8s}{:>7s}\n'.format(s[0], s[1], s[2]))  # {:>11s}{:>16s}   , s[3], s[4]
        g.writelines('\n')
        g.writelines('\n')
        g.writelines('[ pairs ]\n')
        g.writelines(';    ai     aj funct         c0         c1         c2         c3\n')
        for i, l in enumerate(sections['[ pairs ]']):
            s = l.split()
            g.writelines('{:>6s}{:>8s}{:>7s}\n'.format(s[0], s[1], s[2]))  # {:>7s}{:>11s}{:>16s} ,,  ,, s[2] s[3], s[4]
        g.writelines('\n')

        g.writelines('[ angles ]\n')
        g.writelines(';    ai     aj     ak funct         c0         c1         c2         c3\n')
        for i, l in enumerate(sections['[ angles ]']):
            s = l.split()
            g.writelines('{:>6s}{:>8s}{:>8s}{:>7s}\n'.format(s[0], s[1], s[2], '5'))  # {:>15s}{:>14s} ,  , s[4], s[5]
        g.writelines('\n')

        g.writelines('[ dihedrals ]\n')
        g.writelines(
            ';    ai     aj     ak     al funct         c0         c1         c2         c3         c4         c5\n')
        for i, l in enumerate(sections['[ dihedrals ]']):
            s = l.split()
            g.writelines('{:>6s}{:>8s}{:>8s}{:>8s}{:>8s}\n'.format(s[0], s[1], s[2], s[3],
                                                                   '9'))  # {:>12s}{:>12s}{:>12s}{:>12s
            # }{:>12s}{:>12s}  , s[5], s[6], s[7], s[8], s[9], s[10]
        g.writelines('\n')


def topologywriter(func_name):
    """
    writes topology file for gromacs
    :param func_name: takes the name of built system.
    :return: gromacs topology file
    """
    sections = make_dic(f'{func_name}.top')
    with open(f'topol.top', 'w+') as p:
        p.writelines(';\n')
        p.writelines(';\n')
        p.writelines(';--------------------------------------------------------\n')
        p.writelines('; Topoloy file \n')
        p.writelines(';--------------------------------------------------------\n')
        p.writelines(';\n')
        p.writelines(';\n')
        p.writelines('#include "charmm36-jul2022.ff/forcefield.itp"\n')
        p.writelines('#include "interface_silica_emami.ff/ff_interface_emami.itp"\n')
        p.writelines('#include "charmm_sam_connection_point.ff/ff_charmm_connection_point.itp"\n')
        p.writelines('\n')
        p.writelines(';--------------------------------------------------------\n')
        p.writelines('\n')
        p.writelines(f'#include "{func_name}.itp"\n')
        p.writelines('\n')
        p.writelines('\n')
        p.writelines(';--------------------------------------------------------\n')
        p.writelines('\n')
        p.writelines('#ifdef POSRES\n')
        p.writelines('#include "POSRES_SUR.itp"\n')
        p.writelines('#endif\n')
        p.writelines('\n')
        p.writelines(';--------------------------------------------------------\n')
        p.writelines('\n')
        p.writelines('[ system ]\n')
        p.writelines('; Name\n')
        for i, l in enumerate(sections['[ system ]']):
            s = l.split()
            p.writelines('{:>6s}{:>6s}\n'.format(s[0], s[1]))
        p.writelines('\n')
        p.writelines(';--------------------------------------------------------\n')
        p.writelines('\n')
        p.writelines('[ molecules ]\n')
        p.writelines('; Compound    #mols\n')
        for i, l in enumerate(sections['[ molecules ]']):
            s = l.split()
            p.writelines('{:>6s}{:>10s}\n'.format(s[0], s[1]))


def indexwriter(SUR, OTS, OTSS, mid_layer):
    """
    writes index file for the built system for gromacs.
    :param SUR: atom IDs of surface, array
    :param OTS: atom IDs of OTS molecules, array
    :param OTSS: atom IDs of OTS molecules individually, 2D array
    :return: gromacs index file.
    """
    with open('index.ndx', 'w+') as f:
        f.writelines("[ Surface ]\n")
        for i in range(1, len(SUR) + 1):
            if i % 15 == 0:
                f.writelines("{:>6d}\n".format(SUR[i - 1]))
            else:
                f.writelines("{:>6d}".format(SUR[i - 1]))
        f.writelines("\n\n[ Ots ]\n")
        for j in range(1, len(OTS) + 1):
            if j % 15 == 0:
                f.writelines("{:>6d}\n".format(OTS[j - 1]))
            else:
                f.writelines("{:>6d}".format(OTS[j - 1]))
        for k in range(len(OTSS)):
            f.writelines(f"\n\n[ Ots{k + 1} ]\n")
            for m in range(1, len(OTSS[k]) + 1):
                if m % 15 == 0:
                    f.writelines("{:>6d}\n".format(OTSS[k][m - 1]))
                else:
                    f.writelines("{:>6d}".format(OTSS[k][m - 1]))
        f.writelines("\n\n[ midSurface ]\n")
        for l in range(1, len(mid_layer) + 1):
            if l % 15 == 0:
                f.writelines("{:>6d}\n".format(mid_layer[l - 1]))
            else:
                f.writelines("{:>6d}".format(mid_layer[l - 1]))
        f.writelines("\n\n")


def posreswriter(SUR, OTS, mid_layer):
    """
    writes position restrains files
    :param SUR: atom IDs of surface
    :param OTS: atom IDs of OTS molecules
    :return: position restrains files for gromacs
    """
    with open('POSRES_SUR.itp', 'w+') as g:
        g.writelines("[ position_restraints ]\n")
        for p in SUR:
            g.writelines("{:>6d}{:>5d}{:>7d}{:>6d}{:>6d}\n".format(p, 1, 1000, 1000, 1000))

    with open('POSRES_OTS.itp', 'w+') as h:
        h.writelines("[ position_restraints ]\n")
        for q in OTS:
            h.writelines("{:>6d}{:>5d}{:>7d}{:>6d}{:>6d}\n".format(q, 1, 1000, 1000, 1000))
            
    with open('POSRES_MIDSUR.itp', 'w+') as p:
        p.writelines("[ position_restraints ]\n")
        for r in mid_layer:
            p.writelines("{:>6d}{:>5d}{:>7d}{:>6d}{:>6d}\n".format(r, 1, 1000, 1000, 1000))            


def gromacssh(func_name):
    """
    gromacs simple run bash file
    :param func_name: system name
    :return: bash run file
    """
    with open('gromacs_func.sh', 'w+') as j:
        j.writelines('#!/bin/bash')
        j.writelines('\n')
        j.writelines(f'gmx grompp -f minimization.mdp  -c {func_name}.gro -p topol.top -o em.tpr ')
        j.writelines('\n')
        j.writelines('gmx mdrun -v -deffnm em ')
        j.writelines('\n')
        j.writelines(f'gmx grompp -f equilibration.mdp  -c em.gro -p topol.top -o nvt.tpr ')
        j.writelines('\n')
        j.writelines('gmx mdrun -v -deffnm nvt')
        j.writelines('\n')
        j.writelines(f'gmx grompp -f production.mdp  -c nvt.gro -p topol.top -o npt.tpr ')
        j.writelines('\n')
        j.writelines('gmx mdrun -v -deffnm npt')
        j.writelines('\n')
