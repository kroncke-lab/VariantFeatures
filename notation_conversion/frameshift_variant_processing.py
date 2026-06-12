import pandas as pd

# Inputs:

# 1. Coding sequence for the gene
# 2. Amino acid position at which frameshift variant is observed
# 3. Max steps (left and right) from amino acid position to lookup for near-stop codons

# Output:

# Dictionary split into sub-dictionaries 'left' and/or 'right' based on which side of reference amino acid the near-stop codon is found.
# Each sub-dictionary contains ref_codon, ref_aa, ref_aa_position, n_steps_wrt_ref_aa, and new lof_variant

# Method
def frameshift_variant_processing(
    coding_sequence: str,
    aa_position: int,
    n_steps: int
):

    try:
        
        # Initialization
        list_codons_to_stop = [
            'AAA', 'AAG', 'AGA', 'CAA', 'CAG',
            'CGA', 'GAA', 'GAG', 'GGA', 'TAC',
            'TAT', 'TCA', 'TCG', 'TGC', 'TGG',
            'TGT', 'TTA', 'TTG', 'TAA', 'TAG', 'TGA'
        ]

        dict_codon_aa = {
            'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S',
            'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y', 'TGT': 'C', 'TGC': 'C',
            'TGG': 'W', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P',
            'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q',
            'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I',
            'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T',
            'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S',
            'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V',
            'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D',
            'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
            'GGG': 'G', 'TAA': 'X', 'TAG': 'X', 'TGA': 'X'
        }
        
        # Validation of inputs
        if (aa_position < 1) | ((aa_position * 3) > len(coding_sequence)):
            raise Exception
        else:
            dict_return = dict()

        # Iterate up to n+1 times
        for step in range(0, n_steps + 1):

            # Get l, r AA position (int)
            aa_position_l = aa_position - step
            aa_position_r = aa_position + step
        
            # Compute l, r codon start position (int; position)
            codon_base_position_l = max(1, (3 * aa_position_l) - 2)
            codon_base_position_r = min((3 * aa_position_r) - 2, len(coding_sequence) - 2)
        
            # Get l, r codons (string; codon)
            # Three letters for codon; subtract 1 from index range here (since Python index begins with 0)
            codon_l = coding_sequence[codon_base_position_l - 1: codon_base_position_l + 2]
            codon_r = coding_sequence[codon_base_position_r - 1: codon_base_position_r + 2]
        
            # Check if l, r codon is in list of codons one SNV away from stop codon (boolean; flag)
            codon_l_flag = codon_l in list_codons_to_stop
            codon_r_flag = codon_r in list_codons_to_stop
        
            # Check if codon one SNV away from stop codon found on either side
            if codon_l_flag | codon_r_flag:
                if codon_l_flag:
                    
                    dict_left = dict()
                    
                    aa_position_l = (codon_base_position_l + 2) / 3
                    aa_position_l_int = int(aa_position_l) if aa_position_l.is_integer() else aa_position_l
                    dict_left['ref_aa_position'] = aa_position_l_int
                    
                    aa_l = dict_codon_aa.get(codon_l, pd.NA)
                    dict_left['ref_aa'] = aa_l

                    dict_left['n_steps_wrt_ref_aa'] = step
                    dict_left['near_stop_codon'] = codon_l
                    dict_left['lof_variant'] = f'{aa_l}{aa_position_l_int}X'

                    dict_return['left'] = dict_left
                    
                if codon_r_flag:
                    
                    dict_right = dict()
                    
                    aa_position_r = (codon_base_position_r + 2) / 3
                    aa_position_r_int = int(aa_position_r) if aa_position_r.is_integer() else aa_position_r
                    dict_right['ref_aa_position'] = aa_position_r_int
                    
                    aa_r = dict_codon_aa.get(codon_r, pd.NA)
                    dict_right['ref_aa'] = aa_r

                    dict_right['n_steps_wrt_ref_aa'] = step
                    dict_right['near_stop_codon'] = codon_r
                    dict_right['lof_variant'] = f'{aa_r}{aa_position_r_int}X'

                    dict_return['right'] = dict_right
                break

        return dict_return

    except Exception as e:

        # Return NAs if no near-stop codons found
        return {
            'left': {
                'ref_aa_position': aa_position,
                'ref_aa': pd.NA,
                'n_steps_wrt_ref_aa': pd.NA,
                'near_stop_codon': pd.NA,
                'lof_variant': pd.NA
            },
            'right': {
                'ref_aa_position': aa_position,
                'ref_aa': pd.NA,
                'n_steps_wrt_ref_aa': pd.NA,
                'near_stop_codon': pd.NA,
                'lof_variant': pd.NA
            }
        }



# Sample usage

# frameshift_variant_processing(
#     coding_sequence = coding_sequence_kcnh2, # Coding sequence as string
#     aa_position = 1155,
#     n_steps = 20
# )



# Sample output

# {'left': {'codon': 'TCG',
#   'aa_position': 1155,
#   'aa': 'S',
#   'n_steps': 0,
#   'variant': 'S1155X'},
#  'right': {'codon': 'TCG',
#   'aa_position': 1155,
#   'aa': 'S',
#   'n_steps': 0,
#   'variant': 'S1155X'}}


