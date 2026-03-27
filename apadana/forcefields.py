from .fileReading import *
from .parameters import *
import parmed as pmd
from parmed import gromacs, amber, unit as u
from pathlib import Path
apadana_package_directory = Path(__file__).parent
charmm36_directory = apadana_package_directory / 'charmm36-jul2022.ff'
ffnonbonded_charmm36 = charmm36_directory.glob('ffnonbonded.itp')
ffbonded_charmm36 = charmm36_directory.glob('ffbonded.itp')


def interface_emami():
    with open('ff_interface_emami.itp', 'w+') as e:
        txt = interfaceemami()
        e.write(txt)
    emami_sections = make_dic('ff_interface_emami.itp')
    return emami_sections


def charmmcc():
    with open('ff_CHARMM_CCCC.itp', 'w+') as charmmcc:
        ch = charmm_cccc()
        charmmcc.write(ch)
    charmmcc_sections = make_dic('ff_CHARMM_CCCC.itp')
    return charmmcc_sections


def charmm_connecting_point():
    charmm_extra()
    charmm_connecting_sections = make_dic('ff_charmm_connection_point.itp')
    return charmm_connecting_sections


def charmm36():
    charmm36nonbonded_sections = make_dic(ffnonbonded_charmm36)
    charmm36bonded_sections = make_dic(ffbonded_charmm36)
    return charmm36nonbonded_sections, charmm36bonded_sections

