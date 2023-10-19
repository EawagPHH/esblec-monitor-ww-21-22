
params.reads = "$baseDir/data/raw/*_3_{1,2}.fastq.gz"
params.outdir = 'reports'

reads = Channel.fromFilePairs(params.reads)
reads.into { fastqc_reads; reads_for_merging }


process merge_reads {
    input:
    tuple val(sample_id), file(reads_file) from reads_for_merging

    output:
    tuple val(sample_id), file("${sample_id}_merged.fastq.gz") into merged

    script:
    """
    nsearch merge \
        --forward ${sample_id}_3_1.fastq.gz \
        --reverse ${sample_id}_3_2.fastq.gz \
        --out ${sample_id}_merged.fastq.gz
    """

}


process filter_reads {
    input:
    tuple val(sample_id), file(merged_file) from merged

    output:
    file "${sample_id}_filtered.fasta.gz" into filtered

    script:
    """
    nsearch filter \
        --in $merged_file \
        --out ${sample_id}_filtered.fasta.gz
    """
}


process classify_reads {
    input:
    file filtered_file from filtered

    output:
    file 'out_file.csv' into read_tbl

    script:
    """
    #! /usr/bin/env python3

    from itertools import groupby
    from collections import defaultdict
    from braceexpand import braceexpand


    def read_fasta(fasta):
        with open(fasta) as file_handle:
            grouped = groupby(file_handle, lambda x: x[0] == ">")
            for cond, entry in grouped:
                if cond:
                    fasta_id = next(entry)
                    _, seq_iter = next(grouped)
                    seq = ''.join([line.strip() for line in seq_iter]).upper()
                    yield ([fasta_id, seq])


    def reverse_complement(sequence):
        intab = "ATGCYRSWKMBDHVN-"
        outtab = "TACGRYSWMKVHDBN-"
        trans_table = str.maketrans(intab, outtab)
        complement = sequence.translate(trans_table)
        rev_complement = complement[::-1]
        return rev_complement


    f_file = read_fasta("$filtered_file")


    def expand_primers(primer):
        deg_dic = {
            'W': '{A,T}',
            'S': '{C,G}',
            'M': '{A,C}',
            'K': '{G,T}',
            'R': '{A,G}',
            'Y': '{C,T}',
            'B': '{C,G,T}',
            'D': '{A,G,T}',
            'H': '{A,C,G}',
            'V': '{A,C,T}',
            'N': '{A,C,G,T}'}
        expand_template = ""
        for letter in primer:
            if letter in {'A', 'T', 'C', 'G', '-'}:
                expand_template += letter
            else:
                expand_template += deg_dic[letter]
        expanded_primers = braceexpand(expand_template)
        return expanded_primers


    C0aT_templates = [reverse_complement(primer) for primer in
                    expand_primers("AAAACAGGAABGCAAAAWGCCGCAARAAAGGGAAWAAGGG")]
    C0bT_templates = [reverse_complement(primer) for primer in
                    expand_primers("ACYCMCCGTCRTGTAVATAAMTACGATACGGGAGGGCTTA")]
    C1aT_templates = [reverse_complement(primer) for primer in
                    expand_primers("GACATRAGYTCAAGYCCAATACGACGAGSTAAATCYTGAT")]
    C1bT_templates = [reverse_complement(primer) for primer in
                    expand_primers("AGGCYKAACAAYCSATYCAGTTAACCARCCYACTTGTBGG")]
    C2aT_templates = [reverse_complement(primer) for primer in
                    expand_primers("GCARATTGGCGRCGCTGTTATCGCTCAYGGTAATGRCGGC")]
    C2bT_templates = [reverse_complement(primer) for primer in
                    expand_primers("YAGCTCCGGYCYTATCGGCGATAAACCAGCCCGCCGRCAG")]
    C3aT_templates = [reverse_complement(primer) for primer in
                    expand_primers("ACTRACCGANCCKARYTCAAACAGMGTWTGYTGMGTGACT")]
    C3bT_templates = [reverse_complement(primer) for primer in
                    expand_primers("AAMAGACCAATGCTGGAGTTAGMRTAAARWCKYTTAGYGC")]
    C4aT_templates = [reverse_complement(primer) for primer in
                    expand_primers("CCATAGCTGAAGTAATGCGGTTYTCCTTTCAGGCTGATGG")]
    C4bT_templates = [reverse_complement(primer) for primer in
                    expand_primers("CAAGCCGCCGGCGGTATAGGTCGCGAGGTCGAGCAGGYTG")]
    C5aT_templates = [reverse_complement(primer) for primer in
                    expand_primers("AAYAACAGCGTGACGGYTGCCGTCGCCATSAGCGTGAACT")]
    C5bT_templates = [reverse_complement(primer) for primer in
                    expand_primers("GCACGATCTTTTGGCCARATCACCGCGATATCRTTGGTGG")]
    C6aT_templates = [reverse_complement(primer) for primer in
                    expand_primers("TATTYTGAATRACVCCCACRGCCATACCNGGYACATYATA")]
    C6bT_templates = [reverse_complement(primer) for primer in
                    expand_primers("TRGCRGTCACTTTATTAGGYTTCATCACRATYTGYTCTGA")]
    C7aT_templates = [reverse_complement(primer) for primer in
                    expand_primers("TTSGGVACGTTAATCCADGTATGGTYCAGNYTRAGDGGYT")]
    C7bT_templates = [reverse_complement(primer) for primer in
                    expand_primers("GCYTTACCVTCRCGRTABCCCCWGGCGWAATGYKVYKCTT")]
    C8aT_templates = [reverse_complement(primer) for primer in
                    expand_primers("TCTCATCGCCGATCKMGCGGGCAAAAGCCGTCACGCCTYC")]
    C8bT_templates = [reverse_complement(primer) for primer in
                    expand_primers("GAATGGCGGTATTCWGCGTABVTTCAGTGCGATCCAGACG")]
    C9aT_templates = [reverse_complement(primer) for primer in
                    expand_primers("RTTAAGCCAYYCWATTCCNSCYGHRCTRTCRCYATGRAAA")]
    C9bT_templates = [reverse_complement(primer) for primer in
                    expand_primers("CATAYGTKGRRATBGAYYGAGARTTAAGCCAYYCWATTCC")]
    C10aT_templates = [reverse_complement(primer) for primer in
                    expand_primers("TGCTTYTCAATCTCCGCGAGAAGKGCYRCYGTGTTYTTCG")]
    C10bT_templates = [reverse_complement(primer) for primer in
                    expand_primers("AGGRCATCAACKCCRCCGACKCGRTCGTCATGRAAGTGCG")]
    C12aT_templates = [reverse_complement(primer) for primer in
                    expand_primers("ATRTCYTTTTYCCARTCRGGAAAYAARCGYTTYTGDCCRT")]
    C12bT_templates = [reverse_complement(primer) for primer in
                    expand_primers("GCDAYAATTTYTCCYTGHKRTTGAAYNAYCCARCCBGTWA")]
    C13aT_templates = [reverse_complement(primer) for primer in
                    expand_primers("TCCTTGGGGATCGACGATCGCAACACCGATCTSAGCTGCT")]
    C13bT_templates = [reverse_complement(primer) for primer in
                    expand_primers("TGCTCAGGATGAGTTGTGTAATAACTTGACCGACAGAGGC")]
    C14aT_templates = [reverse_complement(primer) for primer in
                    expand_primers("TGCTSACSGABCCWATMTCGAACAGSGTCTGCTCRYTGAC")]
    C14bT_templates = [reverse_complement(primer) for primer in
                    expand_primers("GWTRGGRKAGTTGCGRTTGGCCAGCATGACGAYGSCBATC")]
    C15aT_templates = [reverse_complement(primer) for primer in
                    expand_primers("CAGTTGAGYAAACCATCAGCCCARCCAACAGGAATAACTT")]
    C15bT_templates = [reverse_complement(primer) for primer in
                    expand_primers("TTTTGTTGAGTTTGCCGCTVTCCATRTCGCTTAGCATGGT")]
    C18aT_templates = [reverse_complement(primer) for primer in
                    expand_primers("CTAAATAGCAGGGGTAGCGTCGCCATCACCGTTAACATTG")]
    C18bT_templates = [reverse_complement(primer) for primer in
                    expand_primers("AAGCCGACCTCCCGAACTTTTCTCCAGGGCTTYCAGCTGC")]
    C20aT_templates = [reverse_complement(primer) for primer in
                    expand_primers("WTGYCACCAHHAAYACAGYCGATAAGGCTAAYRCACGCAT")]
    C20bT_templates = [reverse_complement(primer) for primer in
                    expand_primers("GCCCGYTTAAGATYRTTGGTAAAWCCYTGCTGYKTRTTCT")]
    C27aT_templates = [reverse_complement(primer) for primer in
                    expand_primers("ATGGGCAATCAGCTTRTTCATGGCMGTATTRTCGCTGTAC")]
    C27bT_templates = [reverse_complement(primer) for primer in
                    expand_primers("GCACGCTGAGTKTCACCYARGGCAYTGCCCAAYGTSAGAT")]

    with open("out_file.csv", 'w') as of:
        for _, f_line in f_file:
            target = "NA"
            if "AAAGA" in f_line:
                f_line = reverse_complement(f_line)
            try:
                _, f_line = f_line.split("TCTTTTCGCAGGCTGGAGCCCAGGTCTTCCTAT")
                f_line, _ = f_line.split("GAATGAGTGTGCGTGCACTC")
                bc, target = f_line.split("TGGGCCCAATTTTCCGTGACAATTAATT")
                for C0aT in C0aT_templates:
                    if C0aT in target:
                        of.write(bc + "," + "C0aT" + "\\n")
                        continue
                for C0bT in C0bT_templates:
                    if C0bT in target:
                        of.write(bc + "," + "C0bT" + "\\n")
                        continue
                for C1aT in C1aT_templates:
                    if C1aT in target:
                        of.write(bc + "," + "C1aT" + "\\n")
                        continue
                for C1bT in C1bT_templates:
                    if C1bT in target:
                        of.write(bc + "," + "C1bT" + "\\n")
                        continue
                for C2aT in C2aT_templates:
                    if C2aT in target:
                        of.write(bc + "," + "C2aT" + "\\n")
                        continue
                for C2bT in C2bT_templates:
                    if C2bT in target:
                        of.write(bc + "," + "C2bT" + "\\n")
                        continue
                for C3aT in C3aT_templates:
                    if C3aT in target:
                        of.write(bc + "," + "C3aT" + "\\n")
                        continue
                for C3bT in C3bT_templates:
                    if C3bT in target:
                        of.write(bc + "," + "C3bT" + "\\n")
                        continue
                for C4aT in C4aT_templates:
                    if C4aT in target:
                        of.write(bc + "," + "C4aT" + "\\n")
                        continue
                for C4bT in C4bT_templates:
                    if C4bT in target:
                        of.write(bc + "," + "C4bT" + "\\n")
                        continue
                for C5aT in C5aT_templates:
                    if C5aT in target:
                        of.write(bc + "," + "C5aT" + "\\n")
                        continue
                for C5bT in C5bT_templates:
                    if C5bT in target:
                        of.write(bc + "," + "C5bT" + "\\n")
                        continue
                for C6aT in C6aT_templates:
                    if C6aT in target:
                        of.write(bc + "," + "C6aT" + "\\n")
                        continue
                for C6bT in C6bT_templates:
                    if C6bT in target:
                        of.write(bc + "," + "C6bT" + "\\n")
                        continue
                for C7aT in C7aT_templates:
                    if C7aT in target:
                        of.write(bc + "," + "C7aT" + "\\n")
                        continue
                for C7bT in C7bT_templates:
                    if C7bT in target:
                        of.write(bc + "," + "C7bT" + "\\n")
                        continue
                for C8aT in C8aT_templates:
                    if C8aT in target:
                        of.write(bc + "," + "C8aT" + "\\n")
                        continue
                for C8bT in C8bT_templates:
                    if C8bT in target:
                        of.write(bc + "," + "C8bT" + "\\n")
                        continue
                for C9aT in C9aT_templates:
                    if C9aT in target:
                        of.write(bc + "," + "C9aT" + "\\n")
                        continue
                for C9bT in C9bT_templates:
                    if C9bT in target:
                        of.write(bc + "," + "C9bT" + "\\n")
                        continue
                for C10aT in C10aT_templates:
                    if C10aT in target:
                        of.write(bc + "," + "C10aT" + "\\n")
                        continue
                for C10bT in C10bT_templates:
                    if C10bT in target:
                        of.write(bc + "," + "C10bT" + "\\n")
                        continue
                for C12aT in C12aT_templates:
                    if C12aT in target:
                        of.write(bc + "," + "C12aT" + "\\n")
                        continue
                for C12bT in C12bT_templates:
                    if C12bT in target:
                        of.write(bc + "," + "C12bT" + "\\n")
                        continue
                for C13aT in C13aT_templates:
                    if C13aT in target:
                        of.write(bc + "," + "C13aT" + "\\n")
                        continue
                for C13bT in C13bT_templates:
                    if C13bT in target:
                        of.write(bc + "," + "C13bT" + "\\n")
                        continue
                for C14aT in C14aT_templates:
                    if C14aT in target:
                        of.write(bc + "," + "C14aT" + "\\n")
                        continue
                for C14bT in C14bT_templates:
                    if C14bT in target:
                        of.write(bc + "," + "C14bT" + "\\n")
                        continue
                for C15aT in C15aT_templates:
                    if C15aT in target:
                        of.write(bc + "," + "C15aT" + "\\n")
                        continue
                for C15bT in C15bT_templates:
                    if C15bT in target:
                        of.write(bc + "," + "C15bT" + "\\n")
                        continue
                for C18aT in C18aT_templates:
                    if C18aT in target:
                        of.write(bc + "," + "C18aT" + "\\n")
                        continue
                for C18bT in C18bT_templates:
                    if C18bT in target:
                        of.write(bc + "," + "C18bT" + "\\n")
                        continue
                for C20aT in C20aT_templates:
                    if C20aT in target:
                        of.write(bc + "," + "C20aT" + "\\n")
                        continue
                for C20bT in C20bT_templates:
                    if C20bT in target:
                        of.write(bc + "," + "C20bT" + "\\n")
                        continue
                for C27aT in C27aT_templates:
                    if C27aT in target:
                        of.write(bc + "," + "C27aT" + "\\n")
                        continue
                for C27bT in C27bT_templates:
                    if C27bT in target:
                        of.write(bc + "," + "C27bT" + "\\n")
                        continue
                of.write(bc + "," + target + "\\n")
            except ValueError:
                continue
        """
    }


process count_targets {
    publishDir "$params.outdir", mode: 'copy'

    input:
    file 'out_file.csv' from read_tbl

    output:
    file 'sequencing_targets.csv' into targets

    script:
    """
    #! /usr/bin/env Rscript

    library(tidyverse)

    tbl <-
      read_csv("out_file.csv",
      col_names=FALSE)

    tbl_counts <-
      tbl %>%
      count(X1, X2) %>%
      arrange(desc(n))

    bcs <-
      c("GAACACA",
        "ACATTCA",
        "TCGCTAG",
        "ATCATTA",
        "TGTATGT",
        "TAAGATA",
        "CGTTTCA",
        "ACGTTGC",
        "TAGATGA",
        "CCACTAG",
        "GTCTTCT",
        "ACAAAGC",
        "AAAAGGC",
        "ACTGTGT",
        "CTGGTAC",
        "GTTGCTA",
        "CCATTCA",
        "GATATCG",
        "ACAGGAT",
        "CGCATAC",
        "GGACCTA",
        "GAACTGA",
        "AGTCGTG",
        "TTCCAGG",
        "TAGAGCG",
        "AAACCCT",
        "ATCCAGT",
        "TCTTCGT",
        "GTGTCCT",
        "CACAGAT",
        "ATAGAGC",
        "GGTCATG",
        "CCATCTC",
        "TATGCAG",
        "AACCGCT",
        "TGACTCT",
        "AAGCACA",
        "TGTGAAC",
        "AACAGTG",
        "GCCTATA",
        "GTCCTCT",
        "TAGAACG",
        "TGGGTGA",
        "ACACGTG",
        "TGCCCAA",
        "TCGTACA",
        "AACCAAG",
        "GTGCTAT",
        "TCCTCAT",
        "GTCGTAT",
        "TCTGGAA",
        "CAGTAGG",
        "ACGTCAT",
        "ACGGAAC",
        "CGGATGA",
        "CCGCATA",
        "GTGAAGA",
        "AACGTGT",
        "AGAAGAC",
        "TGTAGGG",
        "AATGCCT",
        "GCTTTCT",
        "CATGAAG",
        "CGATCAC",
        "TCGCGTT",
        "GATTGGC",
        "CATCGGT",
        "CTTTCCA",
        "ACCAGAT",
        "GTGATAC",
        "TGCATCC",
        "TGCAAAC",
        "GCAACGA",
        "GCCATAC",
        "TACCTTC",
        "GAATCGA",
        "CTACGTT",
        "GTTTCGG",
        "GCAAATG",
        "CTGACAC",
        "GCACTTT",
        "AGCAAAC",
        "ATCCTCG",
        "CTGCTGT",
        "CAATCAC",
        "GTTCTGG",
        "TTGCTAC",
        "CATACGT",
        "GATCGTT",
        "CACGTTT",
        "AACACCT",
        "AGCTGTA",
        "GGTCTTG",
        "CCGACAT",
        "TATCGTC",
        "ATGGTTG")

    tbl_counts %>%
      filter(X2 %in% c("C0aT", "C0bT", "C1aT", "C1bT", "C2aT", "C2bT", "C3aT", "C3bT", "C4aT", "C4bT", "C5aT", "C5bT", "C6aT", "C6bT", "C7aT", "C7bT", "C8aT", "C8bT", "C9aT", "C9bT", "C10aT", "C10bT", "C12aT", "C12bT", "C13aT", "C13bT", "C14aT", "C14bT", "C15aT", "C15bT", "C18aT", "C18bT", "C20aT", "C20bT", "C27aT", "C27bT"),
             X1 %in% bcs) %>%
      write_csv("sequencing_targets.csv")
    """
    }

