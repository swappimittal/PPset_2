
import streamlit as st
import pandas as pd
import itertools
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp

def get_variant_regions(gblock):
    seq_before_snp = gblock[gblock.index('/') - 12:gblock.index('/') - 2]
    seq_after_snp = gblock[gblock.index('/') + 3:gblock.index('/') + 13]
    variant_1 = gblock[gblock.index('/') - 1]
    variant_2 = gblock[gblock.index('/') + 1]
    seq_1 = seq_before_snp + variant_1 + seq_after_snp
    seq_2 = seq_before_snp + variant_2 + seq_after_snp
    return [seq_1, seq_2]  # Return a list containing the two sequences

def calculate_tm(sequence):

    # Create a sequence object from the input sequence
    dna_seq = Seq(sequence)

    # Calculate the Tm with custom salt and DNA concentrations
    tm = MeltingTemp.Tm_Wallace(dna_seq)

    return tm

def get_valid_permutations():
    valid_permutations = []
    length_range = range(8, 13)

    def has_5_or_6_consecutive_1s(sequence):
        consecutive_count = 0
        for digit in sequence:
            if digit == '1':
                consecutive_count += 1
                if consecutive_count >= 5:
                    return True
            else:
                consecutive_count = 0
        return False

    for length in length_range:
        for perm in itertools.product([0, 1], repeat=length):
            perm_str = ''.join(map(str, perm))
            if sum(perm) <= 6 and '111' in perm_str and not has_5_or_6_consecutive_1s(perm_str):
                valid_permutations.append('0' + perm_str + '0')
    return valid_permutations

def generate_sub_sequences(sequence):
    base_list = [base for base in sequence]
    base_list[10] = "*" + base_list[10]
    base_list[11] = "*" + base_list[11]
    base_list[9] = "*" + base_list[9]
    sub_sequences = []

    for i in range(0, 6):
        sub_sequence = base_list[4+i:14+i]
        sub_sequences.append(sub_sequence)
    for i in range(0, 7):
        sub_sequence = base_list[3+i:14+i]
        sub_sequences.append(sub_sequence)
    for i in range(0, 8):
        sub_sequence = base_list[2+i:14+i]
        sub_sequences.append(sub_sequence)
    for i in range(0, 9):
        sub_sequence = base_list[1+i:14+i]
        sub_sequences.append(sub_sequence)
    for i in range(0, 10):
        sub_sequence = base_list[i:14+i]
        sub_sequences.append(sub_sequence)
    return sub_sequences

def generate_master_probe_list(sub_sequences, valid_permutations):
    master_probe_list = []

    for sub_sequence in sub_sequences:
        sub_sequence_LNA = []
        for perm in valid_permutations:
            if len(sub_sequence) == len(perm):
                modified_sequence = []
                for i in range(len(sub_sequence)):
                    if perm[i] == '0':
                        modified_sequence.append(sub_sequence[i])
                    if perm[i] == '1':
                        modified_sequence.append("+" + sub_sequence[i])
                master_probe_list.append(modified_sequence)

    filtered_master_probe_list = []
    for sub_sequence in master_probe_list:
        exclude_sub_sequence = False
        for base in sub_sequence:
            if base[0] == "*":
                exclude_sub_sequence = True
                break
        if not exclude_sub_sequence:
            filtered_master_probe_list.append(sub_sequence)

    return filtered_master_probe_list


def remove_3G_3C(probe_list):
    filtered_probe_list = []
    for sub_sequence in probe_list:
        exclude_sub_sequence = False
        for i in range(len(sub_sequence) - 2):
            if (
                sub_sequence[i:i + 3] == ['G', 'G', 'G'] or
                sub_sequence[i:i + 3] == ['C', 'C', 'C']
            ):
                exclude_sub_sequence = True
                break
        if not exclude_sub_sequence:
            filtered_probe_list.append(sub_sequence)
    return filtered_probe_list

def remove_3primeG_5primeG(probe_list):
    filtered_probe_list = []
    for sub_sequence in probe_list:
        exclude_sub_sequence = False
        if sub_sequence[0][-1] == 'G' or sub_sequence[-1][-1] == 'G':
            exclude_sub_sequence = True
        if not exclude_sub_sequence:
            filtered_probe_list.append(sub_sequence)
    return filtered_probe_list

def calculate_Tm_values(probe_list):
    tm_dict = {}
    for sub_sequence in probe_list:
        for i, base in enumerate(sub_sequence):
            if len(base) == 3:
                sub_sequence[i] = "+" + "*" + base[-1].lower()
               
          
        base_sequence = "".join(s for s in sub_sequence)
       
        base_sequence = "".join([char for char in base_sequence if char != "*"])
        G_count = base_sequence.count('G')
        C_count = base_sequence.count('C')
        LNA_seq = ""
        for i in range(len(base_sequence)):
            if base_sequence[i] == "+":
                LNA_seq += base_sequence[i + 1]
        base_sequence = "".join([char for char in base_sequence if char != "+"])
        
        tm = calculate_tm(base_sequence) + calculate_tm(LNA_seq) +14
        tm_dict[''.join(sub_sequence)] = tm  # Use ''.join to create a string key
    return tm_dict

def filter_Tm_probes_2(tm_dict, tm_range=(40, 50)):
    probes_to_remove = [probe for probe, tm in tm_dict.items() if not (tm_range[0] <= tm <= tm_range[1])]

    for probe in probes_to_remove:
        del tm_dict[probe]
    return(tm_dict)

def create_probe_parameter_dict(tm_dict):
    probe_para_dict = {}  # Create a new dictionary to store the modified values
    for probe in tm_dict:
        parameter_dict = {"Tm": tm_dict[probe]}  # Create a dictionary with the temperature parameter
        probe_para_dict[probe] = parameter_dict  # Add the dictionary to the new dictionary
    return probe_para_dict  # Return the new dictionary with values as dictionaries

def add_length_parameter(probe_para_dict):
    for probe in probe_para_dict:
        PROBE = probe.upper()
        probe_para_dict[probe]['probe length'] = PROBE.count("A") + PROBE.count("T") + PROBE.count("C") + PROBE.count("G")
    return probe_para_dict
def add_GC_ratio_parameter(probe_para_dict):
    for probe in probe_para_dict:
        PROBE = probe.upper()
        probe_para_dict[probe]['% GC content'] = int(((PROBE.count("G") + PROBE.count("C"))/(PROBE.count("A") + PROBE.count("T") + PROBE.count("C") + PROBE.count("G")))*100)
    return probe_para_dict
def add_LNA_count_parameter(probe_para_dict):
    for probe in probe_para_dict:
        PROBE = probe.upper()
        probe_para_dict[probe]['LNA count'] = PROBE.count("+")
    return probe_para_dict

def add_snp_distance_parameter(probe_para_dict):
    for probe in probe_para_dict:
        PROBE = ''.join([char for char in probe if char != "*" and char != "+"])
        distances = [PROBE.find(nucleotide) for nucleotide in 'atgc']
        max_distance = max(distances)
        snp_dist = max(0, max_distance)
        probe_para_dict[probe]['snp position'] = snp_dist
    return probe_para_dict
  
def filter_Tm_probes(master_probe_list, tm_dict, tm_range=(40, 50)):
    filtered_master_probe_list = []
    probe_counter = 0
    for sublist in master_probe_list:
        probe_sequence = ''.join(s for s in sublist)
        probe_sequence = ''.join([char for char in probe_sequence if char != "*"])
        tm_value = tm_dict.get(''.join(sublist), None)  # Use ''.join for the key
        if tm_value is not None and tm_range[0] <= tm_value <= tm_range[1]:
            filtered_master_probe_list.append(sublist)
            probe_counter += 1
            print(f"Probe {probe_counter} Count: {len(sublist)}, Sequence: {probe_sequence}, Tm: {tm_value:.2f}°C")
    total_probes_printed = len(filtered_master_probe_list)
    print("Total probes printed:", total_probes_printed)


def display_probe_data(probe_dict):
    probe_data = []
    for probe, parameters in probe_dict.items():
        probe = probe.upper()
        probe = ''.join([char for char in probe if char != "*"])
        probe_info = {"Probe": probe}
        probe_info.update(parameters)
        probe_data.append(probe_info)

    # Display as a table
    st.table(probe_data)
def export_probe_data_to_excel(probe_dict, probe_name):
    data = []
    for probe, parameters in probe_dict.items():
        probe = probe.upper()
        probe = ''.join([char for char in probe if char != "*"])
        row = {"Probe": probe}
        row.update(parameters)  # Add the parameters to the row
        data.append(row)
    df = pd.DataFrame(data)

    # Export to an Excel file
    excel_file = f"probe_data_{probe_name}.xlsx"
    df.to_excel(excel_file, index=False)
    st.write(f"Probe data exported to {excel_file}")

# Your functions here (formatted_probe_dict, excel_probe_dict, main, etc.)

def main():
    st.title("Probe generator")

    # Sidebar for user input
    st.sidebar.header("User Input")
    tm_range_input = st.sidebar.text_input("Enter the desired Tm range (e.g., '62 67'):")
    gblock = st.sidebar.text_input("Enter the gblock seq from ELN")

    if not tm_range_input or not gblock:
        st.warning("Please enter both the Tm range and gblock sequence.")
        return

    tm_range = tm_range_input.split()

    if len(tm_range) != 2:
        st.warning("Invalid Tm range format. Please enter two values separated by a space.")
        return

    valid_permutations = get_valid_permutations()
    input_seq = get_variant_regions(gblock)
    seq_1 = input_seq[0]
    seq_2 = input_seq[1]

    # Process seq_1
    sub_sequences_seq1 = generate_sub_sequences(seq_1)
    master_probe_list_seq1 = generate_master_probe_list(sub_sequences_seq1, valid_permutations)
    master_probe_list_seq1 = remove_3primeG_5primeG(remove_3G_3C(master_probe_list_seq1))
    tm_dict_seq1 = calculate_Tm_values(master_probe_list_seq1)
    filtered_probes_seq1 = filter_Tm_probes_2(tm_dict_seq1, (int(tm_range[0]), int(tm_range[1])))
    probe_dict_seq1 = add_LNA_count_parameter(add_snp_distance_parameter(add_GC_ratio_parameter(add_length_parameter(create_probe_parameter_dict(tm_dict_seq1)))))
    # Display probe data and offer Excel export
    st.header("Probes for seq_1")
    display_probe_data(probe_dict_seq1)
    if st.button("Export Excel for Probe 1"):
        export_probe_data_to_excel(probe_dict_seq1, "Probe1")

    # Process seq_2
    sub_sequences_seq2 = generate_sub_sequences(seq_2)
    master_probe_list_seq2 = generate_master_probe_list(sub_sequences_seq2, valid_permutations)
    master_probe_list_seq2 = remove_3primeG_5primeG(remove_3G_3C(master_probe_list_seq2))
    tm_dict_seq2 = calculate_Tm_values(master_probe_list_seq2)
    filtered_probes_seq2 = filter_Tm_probes_2(tm_dict_seq2, (int(tm_range[0]), int(tm_range[1])))
    probe_dict_seq2 = add_LNA_count_parameter(add_snp_distance_parameter(add_GC_ratio_parameter(add_length_parameter(create_probe_parameter_dict(tm_dict_seq2)))))

    # Display probe data and offer Excel export
    st.header("Probes for seq_2")
    display_probe_data(probe_dict_seq2)
    if st.button("Export Excel for Probe 2"):
        export_probe_data_to_excel(probe_dict_seq2, "Probe2")

if __name__ == "__main__":
    main()
