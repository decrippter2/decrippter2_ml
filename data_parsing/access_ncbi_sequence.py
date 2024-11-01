from subprocess import run


def download_gbk(accession, start, end, folder):
    if start != 0 and end != 0:
        try:
            run(
                f"ncbi-acc-download {accession} --range {start}:{end} -v -o {folder}/{accession}.gbk"
            )
        except Exception as e:
            print(accession, "not found")
    else:
        try:
            run(f"ncbi-acc-download {accession} -v -o ncbi_gbk/{accession}.gbk")
        except Exception as e:
            print(accession, "not found")


def download_prot(accession, folder):
    run(f"ncbi-acc-download {accession} -m protein -F fasta -o {folder}/{accession}.fa")
