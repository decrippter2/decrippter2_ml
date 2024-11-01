def parse_prot_seq(file):
    """
    parses sequence from FASTA file
    :param file: str, filename of FASTA sequence
    :return: str, containing solely aminoacid composition
    """
    with open(file) as fasta:
        prot_seq = ""
        fasta_lines = fasta.readlines()
        for line in fasta_lines:
            if not line.startswith(">"):
                prot_seq += line[:-1]
    return prot_seq


def classify_aa_seq(fasta, core_list):
    """
    Extracts leader, core and follower sequences of RiPP
    :param fasta: str containing aa sequence
    :param core: list containing core sequences
    :return: 3 strings, leader core and follower sequences
    """
    if core_list != []:
        for core in core_list:
            if core in fasta:
                start = fasta.index(core)
                leader_seq = fasta[:start]
                core_seq = fasta[start : start + len(core)]
                follower_seq = fasta[start + len(core) :]
                return leader_seq, core_seq, follower_seq
            else:
                continue
        leader_seq = "TBA"
        core_seq = "TBA"
        follower_seq = "TBA"
        return leader_seq, core_seq, follower_seq
    else:
        leader_seq = "TBA"
        core_seq = "TBA"
        follower_seq = "TBA"
    return leader_seq, core_seq, follower_seq
