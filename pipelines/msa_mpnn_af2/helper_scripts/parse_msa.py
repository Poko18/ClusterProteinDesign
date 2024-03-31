import sys, os, glob
import pandas as pd
import numpy as np
import string

num2aa=[
    'ALA','ARG','ASN','ASP','CYS',
    'GLN','GLU','GLY','HIS','ILE',
    'LEU','LYS','MET','PHE','PRO',
    'SER','THR','TRP','TYR','VAL',
    ]   
        
aa2num= {x:i for i,x in enumerate(num2aa)}

alpha_1 = list("ARNDCQEGHILKMFPSTWYV-")
aa_N_1 = {n:a for n,a in enumerate(alpha_1)}
aa_1_N = {a:n for n,a in enumerate(alpha_1)}

aa123 = {aa1: aa3 for aa1, aa3 in zip(alpha_1, num2aa)}
aa321 = {aa3: aa1 for aa1, aa3 in zip(alpha_1, num2aa)}


# minimal sc atom representation (Nx8)
aa2short=[
    (" N  "," CA "," C  "," CB ",  None,  None,  None,  None), # ala
    (" N  "," CA "," C  "," CB "," CG "," CD "," NE "," CZ "), # arg
    (" N  "," CA "," C  "," CB "," CG "," OD1",  None,  None), # asn
    (" N  "," CA "," C  "," CB "," CG "," OD1",  None,  None), # asp
    (" N  "," CA "," C  "," CB "," SG ",  None,  None,  None), # cys
    (" N  "," CA "," C  "," CB "," CG "," CD "," OE1",  None), # gln
    (" N  "," CA "," C  "," CB "," CG "," CD "," OE1",  None), # glu
    (" N  "," CA "," C  ",  None,  None,  None,  None,  None), # gly
    (" N  "," CA "," C  "," CB "," CG "," ND1",  None,  None), # his
    (" N  "," CA "," C  "," CB "," CG1"," CD1",  None,  None), # ile
    (" N  "," CA "," C  "," CB "," CG "," CD1",  None,  None), # leu
    (" N  "," CA "," C  "," CB "," CG "," CD "," CE "," NZ "), # lys
    (" N  "," CA "," C  "," CB "," CG "," SD "," CE ",  None), # met
    (" N  "," CA "," C  "," CB "," CG "," CD1",  None,  None), # phe
    (" N  "," CA "," C  "," CB "," CG "," CD ",  None,  None), # pro
    (" N  "," CA "," C  "," CB "," OG ",  None,  None,  None), # ser
    (" N  "," CA "," C  "," CB "," OG1",  None,  None,  None), # thr
    (" N  "," CA "," C  "," CB "," CG "," CD1",  None,  None), # trp
    (" N  "," CA "," C  "," CB "," CG "," CD1",  None,  None), # tyr
    (" N  "," CA "," C  "," CB "," CG1",  None,  None,  None), # val
]

# full sc atom representation (Nx14)
aa2long=[
    (" N  "," CA "," C  "," O  "," CB ",  None,  None,  None,  None,  None,  None,  None,  None,  None), # ala
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD "," NE "," CZ "," NH1"," NH2",  None,  None,  None), # arg
    (" N  "," CA "," C  "," O  "," CB "," CG "," OD1"," ND2",  None,  None,  None,  None,  None,  None), # asn
    (" N  "," CA "," C  "," O  "," CB "," CG "," OD1"," OD2",  None,  None,  None,  None,  None,  None), # asp
    (" N  "," CA "," C  "," O  "," CB "," SG ",  None,  None,  None,  None,  None,  None,  None,  None), # cys
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD "," OE1"," NE2",  None,  None,  None,  None,  None), # gln
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD "," OE1"," OE2",  None,  None,  None,  None,  None), # glu
    (" N  "," CA "," C  "," O  ",  None,  None,  None,  None,  None,  None,  None,  None,  None,  None), # gly
    (" N  "," CA "," C  "," O  "," CB "," CG "," ND1"," CD2"," CE1"," NE2",  None,  None,  None,  None), # his
    (" N  "," CA "," C  "," O  "," CB "," CG1"," CG2"," CD1",  None,  None,  None,  None,  None,  None), # ile
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2",  None,  None,  None,  None,  None,  None), # leu
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD "," CE "," NZ ",  None,  None,  None,  None,  None), # lys
    (" N  "," CA "," C  "," O  "," CB "," CG "," SD "," CE ",  None,  None,  None,  None,  None,  None), # met
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2"," CE1"," CE2"," CZ ",  None,  None,  None), # phe
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD ",  None,  None,  None,  None,  None,  None,  None), # pro
    (" N  "," CA "," C  "," O  "," CB "," OG ",  None,  None,  None,  None,  None,  None,  None,  None), # ser
    (" N  "," CA "," C  "," O  "," CB "," OG1"," CG2",  None,  None,  None,  None,  None,  None,  None), # thr
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2"," CE2"," CE3"," NE1"," CZ2"," CZ3"," CH2"), # trp
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2"," CE1"," CE2"," CZ "," OH ",  None,  None), # tyr
    (" N  "," CA "," C  "," O  "," CB "," CG1"," CG2",  None,  None,  None,  None,  None,  None,  None), # val
]

def parse_a3m(filename):
    '''read A3M and convert letters into integers in the 0..20 range,
    also keep track of insertions
    '''

    # read A3M file line by line
    lab,seq = [],[] # labels and sequences
    for line in open(filename, "r"):
        if line[0] == '>':
            lab.append(line.split()[0][1:])
            seq.append("")
        else:
            seq[-1] += line.rstrip()

    # parse sequences
    msa,ins = [],[]
    table = str.maketrans(dict.fromkeys(string.ascii_lowercase))
    nrow,ncol = len(seq),len(seq[0])

    for seqi in seq:

        # remove lowercase letters and append to MSA
        msa.append(seqi.translate(table))

        # 0 - match or gap; 1 - insertion
        a = np.array([0 if c.isupper() or c=='-' else 1 for c in seqi])
        i = np.zeros((ncol))

        if np.sum(a) > 0:
            # positions of insertions
            pos = np.where(a==1)[0]

            # shift by occurrence
            a = pos - np.arange(pos.shape[0])

            # position of insertions in the cleaned sequence
            # and their length
            pos,num = np.unique(a, return_counts=True)
            i[pos[pos<ncol]] = num[pos<ncol]

        # append to the matrix of insetions
        ins.append(i)

    # convert letters into numbers
    alphabet = np.array(list("ARNDCQEGHILKMFPSTWYV-"), dtype='|S1').view(np.uint8)
    msa = np.array([list(s) for s in msa], dtype='|U1')
    msa_num = msa.copy().astype('|S1').view(np.uint8)

    for i in range(alphabet.shape[0]):
        msa_num[msa_num == alphabet[i]] = i

    # treat all unknown characters as gaps
    msa[msa_num > 20] = '-'
    msa_num[msa_num > 20] = 20

    ins = np.array(ins, dtype=np.uint8)

    return {"msa":msa, "msa_num":msa_num, "labels":lab, "insertions":ins}

def get_conserved_positions(a3m, frac_conserved, min_count=10):
    aln = parse_a3m(a3m)  # You need to make sure `tools` is properly imported.
    msa = aln['msa']
    msa_num = aln['msa_num']
    L = msa.shape[1]

    counts = np.stack([np.bincount(column, minlength=21) for column in msa_num.T]).T
    max_count = np.max(counts, axis=0)

    freq = counts / msa_num.shape[0]
    freq_norm = freq[:20] / freq[:20].sum(axis=0)
    max_freq_norm = np.max(freq_norm, axis=0)  # frequency of the most frequent AA, excluding gaps
    max_freq_norm[max_count < min_count] = 0  # don't choose positions that have too low counts

    conserved = np.argsort(max_freq_norm)[::-1][:int(L * frac_conserved)] + 1  # make 1-indexed to input to MPNN
    
    return conserved

def process_a3m_file_rf2(sequence, sample_a3m_path, processed_file_name):
    if os.path.isfile(sample_a3m_path):
        with open(sample_a3m_path, 'r') as file:
            lines = file.readlines()

        # Change line 2
        lines[1] = sequence + "\n"
        lines[3] = sequence + "\n"

        # Write the processed lines to the output file
        with open(processed_file_name, 'w') as processed_file:
            processed_file.writelines(lines)
    else:
        raise ValueError("a3m file doesnt exist!")

def process_a3m_file_rf2_to_af2(sequences, sample_a3m_path, processed_file_name):
    if os.path.isfile(sample_a3m_path):
        with open(sample_a3m_path, 'r') as file:
            lines = file.readlines()
        
        # Calculate lengths and nums for the sequences
        lengths = [len(seq) for seq in sequences.split(":")]
        nums = [1] * len(lengths)

        # Create the header line
        header_line = f"#{','.join(map(str, lengths))}\t{','.join(map(str, nums))}\n"
        lines.insert(0, header_line)

        # Swap lines 2 and 4 for the input sequence
        sequence_line = ''.join(map(str, sequences.split(':')))
        lines[2] = sequence_line + '\n'
        lines[4] = sequence_line + '\n'

        # Remove all remaining "/" characters in all lines
        lines = [line.replace('/', '') for line in lines]

        # Write the processed lines to the output file
        with open(processed_file_name, 'w') as processed_file:
            processed_file.writelines(lines)
    else:
        raise ValueError("a3m file doesn't exist!")
    
def process_a3m_file_af2(sequences, sample_a3m_path, processed_file_name):
    if os.path.isfile(sample_a3m_path):
        with open(sample_a3m_path, 'r') as file:
            lines = file.readlines()

        # Swap lines 2 and 4 for the input sequence
        sequence_line = ''.join(map(str, sequences.split(':')))
        lines[2] = sequence_line + '\n'
        #lines[4] = sequence_line + '\n'

        # Remove all remaining "/" characters in all lines
        lines = [line.replace('/', '') for line in lines]

        # Write the processed lines to the output file
        with open(processed_file_name, 'w') as processed_file:
            processed_file.writelines(lines)
    else:
        raise ValueError("a3m file doesn't exist!")