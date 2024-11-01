"""
Script containing all necessary functions to parse .gbk files
"""

import re


def read_gbk(gbk):
    """
    reads .gbk file
    :param gbk: str, contains .gbk filename
    :return: list of strings, containing all read lines from .gbk
    """
    with open(gbk) as gbk_lines:
        gbk_lines = gbk_lines.readlines()
    return gbk_lines


def parse_numbers(s):
    """
    Extracts start and end positions of CDS
    :param s: string containing CDS positions
    :return: start, end values
    """
    match = re.search(r"(\d+)\.\.[^\d]?(\d+)", s)

    num1 = int(match.group(1))
    num2 = int(match.group(2))
    return num1, num2


def find_max_elem(list):
    """
    Finds maximum element index in a list
    :param list: list of values
    :return: maximum value's index in list
    """
    max_val = max(list)
    idx_max = list.index(max_val)
    return idx_max


def parse_lines(gbk_file, ncomp):
    """
    Parses .gbk file extracting shortest proteins
    :param gbk_file: str, contains .gbk filename
    :param ncomp: number of proteins to extract
    :return: prot_lengths, list of integers containing the length of the
            extracted proteins
            prot_ids, list of str containing protein accessions
    """
    prot_ids = [None] * ncomp
    prot_lengths = [None] * ncomp
    gbk_lines = read_gbk(gbk_file)
    flag = False
    for line in gbk_lines:
        line = line.strip()
        if (
            line.startswith("CDS ")
            and "<" not in line
            and ">" not in line
            and "::" not in line
        ):
            line_list = line.split()
            flag = True
            cds = [None, None]
            try:
                num1, num2 = parse_numbers(line_list[1])
                cds[0] = num1
                cds[1] = num2
                aa_length = (int(cds[1]) - int(cds[0]) + 1) / 3 - 1
                aa_length = int(aa_length)
            except AttributeError:
                continue

        elif line.startswith("/protein_id") and flag:
            prot_entry = line.split(sep="=")
            prot_id = prot_entry[1][1:-1]

            for i in range(ncomp):
                if prot_id not in prot_ids:
                    if prot_lengths[i] == None:
                        prot_lengths[i] = aa_length
                        prot_ids[i] = prot_id
                    elif aa_length < prot_lengths[i] and None not in prot_lengths:
                        max_pos = find_max_elem(prot_lengths)
                        prot_lengths[max_pos] = aa_length
                        prot_ids[max_pos] = prot_id
                else:
                    continue
            flag = False
    return prot_lengths, prot_ids
