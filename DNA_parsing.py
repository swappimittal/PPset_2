def reverse_complement(sequence):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', "[" : "[", "/" :"/", "]":"]"}
    reverse_comp = ''.join([complement_dict.get(base, base) for base in reversed(sequence)])
    return reverse_comp

def complement(sequence):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', "[" : "[", "/" :"/", "]":"]"}
    comp_sequence = ''.join([complement_dict[base] for base in sequence])
    return comp_sequence
  
  def clean_up_input(gblock):
    gblock = gblock.replace(" ", "")
    gblock = gblock.upper()
    acceptable_characters = ['a', 'c', 't', 'g', 'A', 'T', 'G', 'C', '[', ']', '/',]
    # Create a translation table to remove characters not in acceptable_characters
    translation_table = str.maketrans("", "", "".join(c for c in gblock if c not in acceptable_characters))
    # Apply the translation table to gblock
    gblock = gblock.translate(translation_table)
    return gblock
    
def get_variant_regions(gblock):
    seq_before_snp = gblock[gblock.index('/') - 12:gblock.index('/') - 2]
    seq_after_snp = gblock[gblock.index('/') + 3:gblock.index('/') + 13]
    variant_1 = gblock[gblock.index('/') - 1]
    variant_2 = gblock[gblock.index('/') + 1]
    seq_1 = seq_before_snp + variant_1 + seq_after_snp
    seq_2 = seq_before_snp + variant_2 + seq_after_snp
    return {seq_1:variant_1, seq_2:variant_2 }  # Return a dict containing the two sequences

def calculate_tm(sequence):

    # Create a sequence object from the input sequence
    dna_seq = Seq(sequence)

    # Calculate the Tm with custom salt and DNA concentrations
    tm = MeltingTemp.Tm_Wallace(dna_seq)

    return tm
