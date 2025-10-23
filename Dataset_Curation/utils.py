import time
import numpy as np

AA3TO1 = {
    'ALA':'A', 'VAL':'V', 'PHE':'F', 'PRO':'P', 'MET':'M',
    'ILE':'I', 'LEU':'L', 'ASP':'D', 'GLU':'E', 'LYS':'K',
    'ARG':'R', 'SER':'S', 'THR':'T', 'TYR':'Y', 'HIS':'H',
    'CYS':'C', 'ASN':'N', 'GLN':'Q', 'TRP':'W', 'GLY':'G',}

restype_name_to_atom14_names = {
        'ALA': ['N', 'CA', 'C', 'O', 'CB', '',    '',    '',    '',    '',    '',    '',    '',    ''],
        'ARG': ['N', 'CA', 'C', 'O', 'CB', 'CG',  'CD',  'NE',  'CZ',  'NH1', 'NH2', '',    '',    ''],
        'ASN': ['N', 'CA', 'C', 'O', 'CB', 'CG',  'OD1', 'ND2', '',    '',    '',    '',    '',    ''],
        'ASP': ['N', 'CA', 'C', 'O', 'CB', 'CG',  'OD1', 'OD2', '',    '',    '',    '',    '',    ''],
        'CYS': ['N', 'CA', 'C', 'O', 'CB', 'SG',  '',    '',    '',    '',    '',    '',    '',    ''],
        'GLN': ['N', 'CA', 'C', 'O', 'CB', 'CG',  'CD',  'OE1', 'NE2', '',    '',    '',    '',    ''],
        'GLU': ['N', 'CA', 'C', 'O', 'CB', 'CG',  'CD',  'OE1', 'OE2', '',    '',    '',    '',    ''],
        'GLY': ['N', 'CA', 'C', 'O', '',   '',    '',    '',    '',    '',    '',    '',    '',    ''],
        'HIS': ['N', 'CA', 'C', 'O', 'CB', 'CG',  'ND1', 'CD2', 'CE1', 'NE2', '',    '',    '',    ''],
        'ILE': ['N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2', 'CD1', '',    '',    '',    '',    '',    ''],
        'LEU': ['N', 'CA', 'C', 'O', 'CB', 'CG',  'CD1', 'CD2', '',    '',    '',    '',    '',    ''],
        'LYS': ['N', 'CA', 'C', 'O', 'CB', 'CG',  'CD',  'CE',  'NZ',  '',    '',    '',    '',    ''],
        'MET': ['N', 'CA', 'C', 'O', 'CB', 'CG',  'SD',  'CE',  '',    '',    '',    '',    '',    ''],
        'PHE': ['N', 'CA', 'C', 'O', 'CB', 'CG',  'CD1', 'CD2', 'CE1', 'CE2', 'CZ',  '',    '',    ''],
        'PRO': ['N', 'CA', 'C', 'O', 'CB', 'CG',  'CD',  '',    '',    '',    '',    '',    '',    ''],
        'SER': ['N', 'CA', 'C', 'O', 'CB', 'OG',  '',    '',    '',    '',    '',    '',    '',    ''],
        'THR': ['N', 'CA', 'C', 'O', 'CB', 'OG1', 'CG2', '',    '',    '',    '',    '',    '',    ''],
        'TRP': ['N', 'CA', 'C', 'O', 'CB', 'CG',  'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
        'TYR': ['N', 'CA', 'C', 'O', 'CB', 'CG',  'CD1', 'CD2', 'CE1', 'CE2', 'CZ',  'OH',  '',    ''],
        'VAL': ['N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2', '',    '',    '',    '',    '',    '',    ''],
        'UNK': ['',  '',   '',  '',  '',   '',    '',    '',    '',    '',    '',    '',    '',    ''],
    }
    
hydrophobic_residues=['V','I','L','M','F','W','C']
charged_residues=['H','R','K','D','E']

def add_cb(input_array):
    #from protein mpnn
    #The virtual Cβ coordinates were calculated using ideal angle and bond length definitions: b = Cα - N, c = C - Cα, a = cross(b, c), Cβ = -0.58273431*a + 0.56802827*b - 0.54067466*c + Cα.
    N,CA,C,O = input_array
    b = CA - N
    c = C - CA
    a = np.cross(b,c)
    CB = np.around(-0.58273431*a + 0.56802827*b - 0.54067466*c + CA,3)
    return CB #np.array([N,CA,C,CB,O])

aaSMILES = {'G':  'NCC(=O)O',
            'A':  'N[C@@]([H])(C)C(=O)O',
            'R':  'N[C@@]([H])(CCCNC(=N)N)C(=O)O',
            'N':  'N[C@@]([H])(CC(=O)N)C(=O)O',
            'D':  'N[C@@]([H])(CC(=O)O)C(=O)O',
            'C':  'N[C@@]([H])(CS)C(=O)O',
            'E':  'N[C@@]([H])(CCC(=O)O)C(=O)O',
            'Q':  'N[C@@]([H])(CCC(=O)N)C(=O)O',
            'H':  'N[C@@]([H])(CC1=CN=C-N1)C(=O)O',
            'I':  'N[C@@]([H])(C(CC)C)C(=O)O',
            'L':  'N[C@@]([H])(CC(C)C)C(=O)O',
            'K':  'N[C@@]([H])(CCCCN)C(=O)O',
            'M':  'N[C@@]([H])(CCSC)C(=O)O',
            'F':  'N[C@@]([H])(Cc1ccccc1)C(=O)O',
            'P':  'N1[C@@]([H])(CCC1)C(=O)O',
            'S':  'N[C@@]([H])(CO)C(=O)O',
            'T':  'N[C@@]([H])(C(O)C)C(=O)O',
            'W':  'N[C@@]([H])(CC(=CN2)C1=C2C=CC=C1)C(=O)O',
            'Y':  'N[C@@]([H])(Cc1ccc(O)cc1)C(=O)O',
            'V':  'N[C@@]([H])(C(C)C)C(=O)O'}

def timer(func):
    def wrapper(*args, **kw):
        time_start=time.time()  
        result = func(*args, **kw)
        time_end=time.time()
        print('time cost',time_end-time_start,'s')
        return result
    return wrapper

def split_task(task_list,n_groups):
    return [task_list[i::n_groups] for i in range(n_groups)]
