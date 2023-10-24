import PPSET_functions.py
def main():
    st.title("Probe generator!")

    # Sidebar for user input
    st.sidebar.header("Input Fields here:")
    input_gblock = st.sidebar.text_input("Enter the gblock seq from ELN")
    gblock = clean_up_input(input_gblock)
    

    tm_range = st.sidebar.slider("Tm", 60, 67, (63, 65), 1)
    GC_range = st.sidebar.slider("GC content (%)", 0, 100, (40, 60), 1)
    pos_range = st.sidebar.slider("SNP position from 5' end", 1, 14, (4, 9), 1)    
    len_range = st.sidebar.slider("probe length", 9, 14, (10, 10), 1)
    LNA_range = st.sidebar.slider("Number of LNA", 3, 7, (3, 6), 1)
    aprox_tm_range = (58, 68)
    
    if not input_gblock:
        st.warning("Please enter the gblock sequence.")
        return
    rev_comp = st.sidebar.checkbox('reverse complement', value=False)
    if rev_comp:
        gblock = reverse_complement(gblock) 
    else:
        gblock = gblock

    valid_permutations = get_valid_permutations()
    
    input_seq = get_variant_regions(gblock)
    seq_1 = list(input_seq.keys())[0]
    seq_2 = list(input_seq.keys())[1]
    
    # Process seq_1
    sub_sequences_seq1 = generate_sub_sequences(seq_1)
    master_probe_list_seq1 = generate_master_probe_list(sub_sequences_seq1, valid_permutations)
    master_probe_list_seq1 = remove_3primeG_5primeG(remove_3G_3C(master_probe_list_seq1))
    tm_dict_seq1 = calculate_Tm_values(master_probe_list_seq1) 
    probe_dict_seq1 = add_LNA_count_parameter(add_snp_distance_parameter(add_GC_ratio_parameter(add_length_parameter(create_probe_parameter_dict(tm_dict_seq1)))))
    filtered_probes_seq1 = filter_aprox_Tm_probes(probe_dict_seq1, (int(aprox_tm_range[0]), int(aprox_tm_range[1])))
    filtered_probes_seq1 = filter_GC_probes(probe_dict_seq1, (int(GC_range[0]), int(GC_range[1])))
    filtered_probes_seq1 = filter_snp_pos(probe_dict_seq1, (int(pos_range[0]), int(pos_range[1])))
    filtered_probes_seq1 = filter_length_probe(probe_dict_seq1, (int(len_range[0]), int(len_range[1])))
    filtered_probes_seq1 = filter_LNA_count_probe(probe_dict_seq1, (int(LNA_range[0]), int(LNA_range[1])))
    filtered_probes_seq1 = refine_Tm_values(probe_dict_seq1, token)
    filtered_probes_seq1 = filter_Tm_probes(probe_dict_seq1, (int(tm_range[0]), int(tm_range[1])))
    get_hairpin_values(probe_dict_seq1, token)
    get_mismatch_values(probe_dict_seq1,input_seq[seq_2], token)
    #get_selfdimer_values(probe_dict_seq1, token)
    # Display probe data and offer Excel export
    st.header("Probes for " + input_seq[seq_1] + " allele")
    display_probe_data(probe_dict_seq1)
    probe_name = f"{input_seq[seq_1]}_allele"
    if st.button("Export probes 1 to Excel"):
        excel_name = 'probes1'
        excel_file = export_probe_data_to_excel(probe_dict_seq1, excel_name)
        st.success(f"Data exported to Excel file: {excel_file}")

        # Create a download button for the Excel file
        with open(excel_file, 'rb') as my_file:
            st.download_button(
                label="Download Excel File",
                data=my_file,
                file_name=f"{excel_name}.xlsx",
                mime='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'
            )
    # Process seq_2
    sub_sequences_seq2 = generate_sub_sequences(seq_2)
    master_probe_list_seq2 = generate_master_probe_list(sub_sequences_seq2, valid_permutations)
    master_probe_list_seq2 = remove_3primeG_5primeG(remove_3G_3C(master_probe_list_seq2))
    tm_dict_seq2 = calculate_Tm_values(master_probe_list_seq2)
    
    probe_dict_seq2 = add_LNA_count_parameter(add_snp_distance_parameter(add_GC_ratio_parameter(add_length_parameter(create_probe_parameter_dict(tm_dict_seq2)))))
    filtered_probes_seq2 = filter_aprox_Tm_probes(probe_dict_seq2, (int(aprox_tm_range[0]), int(aprox_tm_range[1])))
    filtered_probes_seq2 = filter_GC_probes(probe_dict_seq2, (int(GC_range[0]), int(GC_range[1])))
    filtered_probes_seq2 = filter_snp_pos(probe_dict_seq2, (int(pos_range[0]), int(pos_range[1])))
    filtered_probes_seq2 = filter_length_probe(probe_dict_seq2, (int(len_range[0]), int(len_range[1])))
    filtered_probes_seq2 = filter_LNA_count_probe(probe_dict_seq2, (int(LNA_range[0]), int(LNA_range[1])))
    filtered_probes_seq2 = refine_Tm_values(probe_dict_seq2, token)
    filtered_probes_seq2 = filter_Tm_probes(probe_dict_seq2, (int(tm_range[0]), int(tm_range[1])))
    get_hairpin_values(probe_dict_seq2, token)
    get_mismatch_values(probe_dict_seq2,input_seq[seq_1], token)
    #get_selfdimer_values(probe_dict_seq2, token)
    # Display probe data and offer Excel export
    # Display probe data and offer Excel export
    st.header("Probes for " + input_seq[seq_2] + " allele")
    display_probe_data(probe_dict_seq2)
    if st.button("Export probes 2 to Excel"):
        excel_name = 'probes2'
        excel_file = export_probe_data_to_excel(probe_dict_seq2, excel_name)
        st.success(f"Data exported to Excel file: {excel_file}")

        # Create a download button for the Excel file
        with open(excel_file, 'rb') as my_file:
            st.download_button(
                label="Download Excel File",
                data=my_file,
                file_name=f"{excel_name}.xlsx",
                mime='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'
            )

if __name__ == "__main__":
    client_id = "swapnil.mittal"
    client_secret = "f669c31a-1817-49ea-b8d4-d666dd1fb8bf"
    idt_username = "Swappi.mittal"
    idt_password = "Swappi_IDT"
    token = get_access_token(client_id, client_secret, idt_username, idt_password)
    main()
