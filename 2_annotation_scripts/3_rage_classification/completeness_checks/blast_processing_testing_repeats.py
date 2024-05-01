import sys

# usage: python3 blast_processing_testing_repeats.py ../traa_test/nucl_blastp_1align_1e-30.txt > test.txt
# how do we want to output this information?

def process_blast_output(blast_file):
    # awk '/^>/ {header = "\"" substr($0, 2) "\""; getline; print header ": " length($0)}' all_full_tra_proteins.fna # code used to make the dict
    gene_lengths = {
        "Full-length_GILLIAM_00070_traA": 2868
        "Full-length_GILLIAM_00148_traC": 2490
        "Full-length_GILLIAM_00149_traW": 660
        "Full-length_GILLIAM_00153_traN": 1692
        "Full-length_GILLIAM_00155_traH": 1302
        "Full-length_GILLIAM_00156_traG": 2625
        "Full-length_GILLIAM_00160_traI": 2670
        "Full-length_GILLIAM_00161_traA": 2820
        "Full-length_GILLIAM_00300_traA": 2715
        "Full-length_GILLIAM_00302_traI": 2670
        "Full-length_GILLIAM_00333_traN": 1689
        "Full-length_GILLIAM_00338_traH": 1302
        "Full-length_GILLIAM_00339_traG": 2652
        "Full-length_GILLIAM_00376_traI": 2670
        "Full-length_GILLIAM_00378_traA": 2733
        "Full-length_GILLIAM_00427_traN": 1689
        "Full-length_GILLIAM_00456_traI": 2670
        "Full-length_GILLIAM_00530_traW": 645
        "Full-length_GILLIAM_00534_traN": 1695
        "Full-length_GILLIAM_00537_traA": 2841
        "Full-length_GILLIAM_00659_traI": 2670
        "Full-length_GILLIAM_00814_traC": 2490
        "Full-length_GILLIAM_00815_traW": 660
        "Full-length_GILLIAM_00819_traN": 1692
        "Full-length_GILLIAM_00821_traH": 1302
        "Full-length_GILLIAM_00822_traG": 2652
        "Full-length_GILLIAM_00826_traI": 2670
        "Full-length_GILLIAM_00827_traA": 2877
        "Full-length_GILLIAM_00882_traA": 2856
        "Full-length_GILLIAM_00883_traI": 2670
        "Full-length_GILLIAM_01182_traC": 2490
        "Full-length_GILLIAM_01183_traW": 660
        "Full-length_GILLIAM_01187_traN": 1692
        "Full-length_GILLIAM_01189_traH": 1302
        "Full-length_GILLIAM_01190_traG": 2640
        "Full-length_GILLIAM_01194_traI": 2670
        "Full-length_GILLIAM_01329_traI": 2670
        "Full-length_GILLIAM_01490_traI": 2610
        "Full-length_GILLIAM_01492_traA": 2733
        "Full-length_GILLIAM_01508_traW": 660
        "Full-length_GILLIAM_01580_traG": 2646
        "Full-length_GILLIAM_01868_traA": 2820
        "Full-length_GILLIAM_01869_traI": 2664
        "Full-length_GILLIAM_01873_traG": 2526
        "Full-length_GILLIAM_01874_traH": 1302
        "Full-length_GILLIAM_01878_traN": 1689
        "Full-length_GILLIAM_01882_traW": 660
        "Full-length_GILLIAM_01883_traC": 2490
        "Full-length_GILLIAM_01924_traN": 1689
        "Full-length_GILLIAM_01928_traW": 660
        "Full-length_GILLIAM_02054_traA": 2841
        "Full-length_GILLIAM_02055_traI": 2670
        "Full-length_GILLIAM_02065_traN": 1692
        "Full-length_GILLIAM_02069_traW": 660
        "Full-length_GILLIAM_02131_traC": 2490
        "Full-length_GILLIAM_02132_traW": 660
        "Full-length_GILLIAM_02136_traN": 1692
        "Full-length_GILLIAM_02138_traH": 1302
        "Full-length_GILLIAM_02139_traG": 2649
        "Full-length_GILLIAM_02143_traI": 2670
        "Full-length_GILLIAM_02144_traA": 2871
        "Full-length_GILLIAM_02466_traA": 2868
        "Full-length_GILLIAM_02473_traG": 2670
        "Full-length_GILLIAM_02479_traN": 1689
        "Full-length_GILLIAM_02482_traW": 660
        "Full-length_GILLIAM_02560_traA": 2859
        "Full-length_GILLIAM_02561_traI": 2670
        "Full-length_GILLIAM_02574_traW": 660
        "Full-length_GILLIAM_02648_traA": 2847
        "Full-length_GILLIAM_02649_traI": 2670
        "Full-length_GILLIAM_02656_traG": 2646
        "Full-length_GILLIAM_02657_traH": 1302
        "Full-length_GILLIAM_02659_traN": 1692
        "Full-length_GILLIAM_02663_traW": 660
        "Full-length_GILLIAM_02664_traC": 2490
        "Full-length_Karp_00081_traN": 1695
        "Full-length_Karp_00187_traW": 660
        "Full-length_Karp_00191_traN": 1695
        "Full-length_Karp_00199_traG": 2583
        "Full-length_Karp_00205_traA": 2838
        "Full-length_Karp_00297_traI": 2424
        "Full-length_Karp_00468_traC": 2493
        "Full-length_Karp_00469_traW": 660
        "Full-length_Karp_00473_traN": 1695
        "Full-length_Karp_00475_traH": 1302
        "Full-length_Karp_00483_traA": 2868
        "Full-length_Karp_00533_traA": 2820
        "Full-length_Karp_00539_traG": 2640
        "Full-length_Karp_00544_traN": 1695
        "Full-length_Karp_00548_traW": 660
        "Full-length_Karp_00549_traC": 2493
        "Full-length_Karp_00686_traW": 660
        "Full-length_Karp_00788_traN": 1695
        "Full-length_Karp_00994_traG": 2652
        "Full-length_Karp_01000_traN": 1692
        "Full-length_Karp_01043_traC": 2493
        "Full-length_Karp_01044_traW": 660
        "Full-length_Karp_01051_traH": 1302
        "Full-length_Karp_01052_traG": 2628
        "Full-length_Karp_01056_traI": 2670
        "Full-length_Karp_01057_traA": 2820
        "Full-length_Karp_01136_traW": 660
        "Full-length_Karp_01137_traC": 2493
        "Full-length_Karp_01309_traG": 2664
        "Full-length_Karp_01314_traN": 1695
        "Full-length_Karp_01318_traW": 660
        "Full-length_Karp_01389_traN": 1692
        "Full-length_Karp_01434_traA": 2841
        "Full-length_Karp_01451_traW": 660
        "Full-length_Karp_01452_traC": 2493
        "Full-length_Karp_01553_traG": 2637
        "Full-length_Karp_01559_traN": 1692
        "Full-length_Karp_01563_traW": 660
        "Full-length_Karp_01564_traC": 2493
        "Full-length_Karp_01666_traC": 2493
        "Full-length_Karp_01667_traW": 660
        "Full-length_Karp_01671_traN": 1692
        "Full-length_Karp_01673_traH": 1302
        "Full-length_Karp_01674_traG": 2637
        "Full-length_Karp_01680_traA": 2820
        "Full-length_Karp_01714_traG": 2658
        "Full-length_Karp_01715_traH": 1302
        "Full-length_Karp_01717_traN": 1695
        "Full-length_Karp_01833_traN": 1692
        "Full-length_Karp_01933_traW": 660
        "Full-length_Karp_02084_traH": 1302
        "Full-length_Karp_02092_traW": 660
        "Full-length_Karp_02286_traW": 660
        "Full-length_Karp_02290_traN": 1689
        "Full-length_Karp_02292_traN": 1695
        "Full-length_Karp_02328_traA": 2823
        "Full-length_Karp_02462_traA": 2802
        "Full-length_Karp_02463_traI": 2670
        "Full-length_Karp_02467_traG": 2634
        "Full-length_Karp_02468_traH": 1302
        "Full-length_Karp_02470_traN": 1695
        "Full-length_Karp_02474_traW": 660
        "Full-length_Karp_02475_traC": 2490
        "Full-length_KATO_00037_traA": 2820
        "Full-length_KATO_00038_traI": 2664
        "Full-length_KATO_00042_traG": 2652
        "Full-length_KATO_00043_traH": 1302
        "Full-length_KATO_00045_traN": 1692
        "Full-length_KATO_00049_traW": 660
        "Full-length_KATO_00050_traC": 2493
        "Full-length_KATO_00116_traC": 2493
        "Full-length_KATO_00117_traW": 660
        "Full-length_KATO_00124_traH": 1302
        "Full-length_KATO_00125_traG": 2640
        "Full-length_KATO_00129_traI": 2664
        "Full-length_KATO_00273_traW": 660
        "Full-length_KATO_00277_traN": 1692
        "Full-length_KATO_00346_traI": 2478
        "Full-length_KATO_00350_traG": 2640
        "Full-length_KATO_00354_traN": 1695
        "Full-length_KATO_00358_traW": 660
        "Full-length_KATO_00359_traC": 2490
        "Full-length_KATO_00432_traA": 2868
        "Full-length_KATO_00516_traI": 2670
        "Full-length_KATO_00560_traI": 2670
        "Full-length_KATO_00564_traG": 2553
        "Full-length_KATO_00568_traN": 1695
        "Full-length_KATO_00763_traN": 1695
        "Full-length_KATO_00767_traG": 2640
        "Full-length_KATO_00811_traI": 2670
        "Full-length_KATO_00816_traG": 2634
        "Full-length_KATO_00994_traI": 2670
        "Full-length_KATO_01090_traW": 615
        "Full-length_KATO_01139_traG": 2652
        "Full-length_KATO_01144_traN": 1695
        "Full-length_KATO_01250_traI": 2670
        "Full-length_KATO_01254_traG": 2658
        "Full-length_KATO_01258_traN": 1692
        "Full-length_KATO_01262_traW": 660
        "Full-length_KATO_01263_traC": 2490
        "Full-length_KATO_01331_traN": 1692
        "Full-length_KATO_01392_traC": 2490
        "Full-length_KATO_01393_traW": 660
        "Full-length_KATO_01444_traI": 2664
        "Full-length_KATO_01617_traC": 2490
        "Full-length_KATO_01618_traW": 660
        "Full-length_KATO_01622_traN": 1692
        "Full-length_KATO_01628_traG": 2553
        "Full-length_KATO_01632_traI": 2670
        "Full-length_KATO_01679_traC": 2490
        "Full-length_KATO_01680_traW": 660
        "Full-length_KATO_01684_traN": 1692
        "Full-length_KATO_01690_traG": 2640
        "Full-length_KATO_01811_traA": 2868
        "Full-length_KATO_01812_traI": 2664
        "Full-length_KATO_01817_traG": 2628
        "Full-length_KATO_01818_traH": 1302
        "Full-length_KATO_01820_traN": 1692
        "Full-length_KATO_01824_traW": 660
        "Full-length_KATO_01825_traC": 2490
        "Full-length_KATO_01841_traI": 2670
        "Full-length_KATO_02049_traA": 2832
        "Full-length_KATO_02165_traW": 660
        "Full-length_KATO_02169_traN": 1692
        "Full-length_KATO_02171_traH": 1302
        "Full-length_KATO_02172_traG": 2628
        "Full-length_KATO_02176_traI": 2673
        "Full-length_KATO_02177_traA": 2829
        "Full-length_KATO_02301_traA": 2841
        "Full-length_KATO_02311_traG": 2619
        "Full-length_KATO_02320_traN": 1689
        "Full-length_KATO_02360_traI": 4101
        "Full-length_KATO_02406_traN": 1695
        "Full-length_KATO_02410_traW": 660
        "Full-length_TA686_00147_traN": 1695
        "Full-length_TA686_00152_traG": 2637
        "Full-length_TA686_00195_traA": 2820
        "Full-length_TA686_00414_traG": 2637
        "Full-length_TA686_00476_traA": 2847
        "Full-length_TA686_00483_traG": 2637
        "Full-length_TA686_00509_traA": 2820
        "Full-length_TA686_00510_traI": 2679
        "Full-length_TA686_00515_traG": 2640
        "Full-length_TA686_00519_traN": 1695
        "Full-length_TA686_00523_traW": 660
        "Full-length_TA686_00524_traC": 2490
        "Full-length_TA686_00700_traI": 2673
        "Full-length_TA686_01254_traI": 2670
        "Full-length_TA686_01293_traC": 2490
        "Full-length_TA686_01294_traW": 660
        "Full-length_TA686_01298_traN": 1689
        "Full-length_TA686_01302_traH": 1302
        "Full-length_TA686_01303_traG": 2643
        "Full-length_TA686_01438_traN": 1689
        "Full-length_TA686_01556_traW": 660
        "Full-length_TA686_01561_traN": 1692
        "Full-length_TA686_01567_traG": 2652
        "Full-length_TA686_01573_traI": 2673
        "Full-length_TA686_01594_traC": 2493
        "Full-length_TA686_01595_traW": 660
        "Full-length_TA686_01599_traN": 1689
        "Full-length_TA686_01605_traG": 2649
        "Full-length_TA686_01609_traI": 2679
        "Full-length_TA686_01610_traA": 2847
        "Full-length_TA686_01659_traW": 660
        "Full-length_TA686_01663_traN": 1689
        "Full-length_TA686_02223_traG": 2784
        "Full-length_TA686_02281_traN": 1695
        "Full-length_TA686_02435_traC": 2493
        "Full-length_TA686_02436_traW": 660
        "Full-length_TA686_02440_traN": 1695
        "Full-length_TA686_02442_traH": 1302
        "Full-length_TA686_02443_traG": 2526
        "Full-length_TA686_02448_traI": 2619
        "Full-length_TA686_02449_traA": 2820
        "Full-length_TA686_02562_traN": 1545
        "Full-length_UT176_00062_traN": 1692
        "Full-length_UT176_00067_traG": 2637
        "Full-length_UT176_00077_traA": 2829
        "Full-length_UT176_00167_traI": 2472
        "Full-length_UT176_00220_traC": 2493
        "Full-length_UT176_00221_traW": 660
        "Full-length_UT176_00225_traN": 1692
        "Full-length_UT176_00232_traG": 2535
        "Full-length_UT176_00238_traI": 2664
        "Full-length_UT176_00239_traA": 2829
        "Full-length_UT176_00339_traA": 2820
        "Full-length_UT176_00507_traI": 2418
        "Full-length_UT176_00870_traC": 2484
        "Full-length_UT176_00961_traC": 2493
        "Full-length_UT176_00962_traW": 660
        "Full-length_UT176_00975_traH": 1302
        "Full-length_UT176_01010_traN": 1692
        "Full-length_UT176_01012_traH": 1302
        "Full-length_UT176_01235_traI": 2664
        "Full-length_UT176_01428_traG": 2544
        "Full-length_UT176_01438_traC": 2484
        "Full-length_UT176_01587_traC": 2484
        "Full-length_UT176_01599_traG": 2664
        "Full-length_UT176_01606_traI": 2664
        "Full-length_UT176_01777_traA": 2820
        "Full-length_UT176_01783_traG": 2634
        "Full-length_UT176_01787_traN": 1692
        "Full-length_UT176_01791_traW": 660
        "Full-length_UT176_01792_traC": 2493
        "Full-length_UT176_02029_traA": 2820
        "Full-length_Boryong_00258_traN": 1701
        "Full-length_Boryong_00723_traN": 1701
        "Full-length_Boryong_00815_traW": 660
        "Full-length_Boryong_01021_traN": 1701
        "Full-length_Boryong_01263_traN": 1701
        "Full-length_Boryong_01491_traN": 1701
        "Full-length_Boryong_01496_traW": 660
        "Full-length_Boryong_01707_traN": 1701
        "Full-length_Boryong_01712_traW": 660
        "Full-length_Boryong_01982_traA": 2487
        "Full-length_Boryong_02264_traA": 2487
        "Full-length_Boryong_02350_traN": 1701
        "Full-length_Boryong_02355_traW": 660
        "Full-length_Ikeda_01647_traA": 2823
        "Full-length_Ikeda_01451_traN": 1692
        "Full-length_Ikeda_01447_traW": 660
        "Full-length_Ikeda_00796_traW": 660
        "Full-length_Ikeda_00561_traH": 1302
        "Full-length_Ikeda_00457_traH": 1302
        "Full-length_Ikeda_00239_traN": 1692
        "Full-length_Ikeda_00117_traW": 660
        "Full-length_Ikeda_01823_traN": 1692
        "Full-length_Ikeda_01936_traW": 660
        "Full-length_Ikeda_01940_traN": 1692
        "Full-length_Ikeda_02064_traC": 2490
        "Full-length_Ikeda_02065_traW": 660
        "Full-length_Ikeda_02069_traN": 1695
        "Full-length_Ikeda_02073_traG": 2628
        "Full-length_Ikeda_02164_traN": 1695
        "Full-length_UT76HP_00205_traA": 2856
        "Full-length_UT76HP_00212_traG": 2640
        "Full-length_UT76HP_00248_traA": 2820
        "Full-length_UT76HP_00249_traI": 2670
        "Full-length_UT76HP_00255_traG": 2652
        "Full-length_UT76HP_00256_traH": 1302
        "Full-length_UT76HP_00258_traN": 1692
        "Full-length_UT76HP_00262_traW": 660
        "Full-length_UT76HP_00263_traC": 2493
        "Full-length_UT76HP_00273_traC": 2493
        "Full-length_UT76HP_00274_traW": 660
        "Full-length_UT76HP_00278_traN": 1689
        "Full-length_UT76HP_00284_traG": 2640
        "Full-length_UT76HP_00288_traI": 2670
        "Full-length_UT76HP_00289_traA": 2820
        "Full-length_UT76HP_00508_traI": 2478
        "Full-length_UT76HP_00594_traW": 660
        "Full-length_UT76HP_00623_traW": 660
        "Full-length_UT76HP_00627_traN": 1692
        "Full-length_UT76HP_00633_traG": 2637
        "Full-length_UT76HP_00639_traA": 2820
        "Full-length_UT76HP_00656_traH": 1302
        "Full-length_UT76HP_00657_traG": 2640
        "Full-length_UT76HP_00661_traI": 2664
        "Full-length_UT76HP_00662_traA": 2817
        "Full-length_UT76HP_00721_traH": 1302
        "Full-length_UT76HP_01034_traA": 2817
        "Full-length_UT76HP_01049_traN": 1692
        "Full-length_UT76HP_01154_traW": 660
        "Full-length_UT76HP_01197_traG": 2442
        "Full-length_UT76HP_01266_traG": 2505
        "Full-length_UT76HP_01384_traA": 2820
        "Full-length_UT76HP_01446_traA": 2820
        "Full-length_UT76HP_01504_traC": 2493
        "Full-length_UT76HP_01505_traW": 660
        "Full-length_UT76HP_01514_traG": 2652
        "Full-length_UT76HP_01827_traN": 1614
        "Full-length_UT76HP_01916_traG": 2658
        "Full-length_UT76HP_01917_traH": 1302
        "Full-length_UT76HP_02045_traC": 2490
        "Full-length_UT76HP_02046_traW": 660
        "Full-length_UT76HP_02055_traG": 2658
        "Full-length_UT76HP_02060_traI": 2664
        "Full-length_UT76HP_02061_traA": 2820
        "Full-length_UT76HP_02207_traC": 2622
        "Full-length_UT76HP_02208_traN": 1695
        "Full-length_UT76HP_02244_traN": 1698
    }

    gene_matches = {}

    with open(blast_file, 'r') as f:
        for line in f:
            columns = line.strip().split('\t')
            gene_name_match = columns[1]
            gene_name_output = columns[0]
            start = int(columns[8])
            end = int(columns[9])

            # Create a list of zeros of the gene length
            gene_length = gene_lengths.get(gene_name_match)
            if gene_length is None:
                continue
            gene_list = [0] * gene_length

            # Mark the matching positions as 1s
            for i in range(start - 1, end):
                gene_list[i] = 1

            # Calculate per_match and completeness
            matches = gene_list.count(1)
            per_match = matches / gene_length
            completeness = 'complete' if per_match >= 0.95 else 'truncated'

            # Update gene_matches dictionary
            if gene_name_output in gene_matches:
                # Modify existing entry
                existing_entry = gene_matches[gene_name_output]
                existing_list = existing_entry['matches']
                for i in range(start - 1, end):
                    if i < len(existing_list):
                        existing_list[i] = 1
                existing_entry['matches'] = existing_list
                existing_entry['per_match'] = existing_list.count(1) / gene_length
                existing_entry['completeness'] = 'complete' if existing_entry['per_match'] >= 0.95 else 'truncated'
            else:
                # Create new entry
                gene_matches[gene_name_output] = {
                    'matches': gene_list,
                    'per_match': per_match,
                    'completeness': completeness
                }

    return gene_matches

def main():
    if len(sys.argv) != 2:
        print("Usage: python script.py blast_output.txt")
        sys.exit(1)

    blast_file = sys.argv[1]
    gene_matches = process_blast_output(blast_file)
    for gene, data in gene_matches.items():
        #print(f"Gene: {gene}")
        ##print(f"Matches: {data['matches']}")
        #print(f"Per Match: {data['per_match']}")
        #print(f"Completeness: {data['completeness']}")
        print(f"Gene: {data['gene']}", end=" | ")
        #print(f"Matches: {data['matches']}", end=" | ")
        print(f"Per Match: {data['per_match']}", end=" | ")
        print(f"Completeness: {data['completeness']}")
        print()

if __name__ == "__main__":
    main()
