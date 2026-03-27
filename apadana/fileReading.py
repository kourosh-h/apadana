HEADERS = list([
    '[ molecules ]',
    '[ system ]',
    '[ angles ]',
    '[ dihedrals ]',
    '[ pairs ]',
    '[ bonds ]',
    '[ atoms ]',
    '[ moleculetype ]',
    '[ atomtypes ]',
    '[ pairtypes ]',
    '[ angletypes ]',
    '[ dihedraltypes ]',
    '[ bondtypes ]',
    '[ defaults ]'])

defaultsTXT = {'[ defaults ]': '1	2	yes	1.000000	1.000000'}


def read_topology(file):
    """
    :param file: text files; gromacs based
    :return: removes comments and spaces
    """
    with open(file) as f:
        for line in f:
            line1 = line
            line = line.partition(';')[0].strip()
            if line:
                yield line1.strip()


def make_dic(file):
    """
    Takes the output of read_topology
    :param file: takes the input text file this ; indicating comments
    :return: a dictionary of sections
    """
    f = list(read_topology(file))
    starts = [i for i, line in enumerate(f) if line in HEADERS]
    starts += [None]
    f_sections = {f[l]: f[l + 1:starts[i + 1]] for i, l in enumerate(starts[:-1])}
    return f_sections
