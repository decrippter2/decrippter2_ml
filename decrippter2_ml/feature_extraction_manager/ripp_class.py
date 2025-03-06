
import contextlib

from numpy import log2


class RiPP:
    def __init__(self, sequence):
        self.sequence = sequence
        self.cys30 = ""
        self.cys20 = ""
        self.cys_ser30 = ""
        self.cys_ser20 = ""
        self.aafreq = {
            "A": 0,
            "R": 0,
            "N": 0,
            "D": 0,
            "C": 0,
            "E": 0,
            "Q": 0,
            "G": 0,
            "H": 0,
            "I": 0,
            "L": 0,
            "K": 0,
            "M": 0,
            "F": 0,
            "P": 0,
            "S": 0,
            "T": 0,
            "W": 0,
            "Y": 0,
            "V": 0,
        }
        self.clfreq = {
            "RHK": 0.0,
            "DE": 0.0,
            "STNQ": 0.0,
            "CGP": 0.0,
            "AVIL": 0.0,
            "MFYW": 0.0,
        }
        self.charge = ""
        self.avgcharge = ""
        self.avghydrop = ""
        self.length = len(self.sequence)
        self.entropy = ""
        self.entropyratio = ""
        self.boman = ""
        self.instability = ""
        self.aliphatic = ""

    def calculate_features(self):
        """
        Wimley-White whole residue hydrophobicity interface scale
        """
        hydrophobicity_dict = {
            "A": 0.17,
            "R": 0.81,
            "N": 0.42,
            "D": 1.23,
            "C": -0.24,
            "Q": 0.58,
            "E": 2.02,
            "G": 0.01,
            "H": 0.96,
            "I": -0.31,
            "L": -0.56,
            "K": 0.99,
            "M": -0.23,
            "F": -1.13,
            "P": 0.45,
            "S": 0.13,
            "T": 0.14,
            "W": -1.85,
            "Y": -0.94,
            "V": 0.07,
            "X": 0.0,
        }

        charge_dict = {
            "A": 0.0,
            "R": 1.0,
            "N": 0.0,
            "D": -1.0,
            "C": 0.0,
            "Q": 0.0,
            "E": -1.0,
            "G": 0.0,
            "H": 0.5,
            "I": 0.0,
            "L": 0.0,
            "K": 1.0,
            "M": 0.0,
            "F": 0.0,
            "P": 0.0,
            "S": 0.0,
            "T": 0.0,
            "W": 0.0,
            "Y": 0.0,
            "V": 0.0,
            "X": 0.0,
        }

        boman_dict = {
            "A": 1.81,
            "R": -14.92,
            "N": -6.64,
            "D": -8.72,
            "C": 1.28,
            "Q": -5.54,
            "E": -6.81,
            "G": 0.94,
            "H": -4.66,
            "I": 4.92,
            "L": 4.92,
            "K": -5.55,
            "M": 2.35,
            "F": 2.98,
            "P": 0.0,
            "S": -3.4,
            "T": -2.57,
            "W": 2.33,
            "Y": -0.14,
            "V": 4.04,
        }

        instability_dict = {
            "W": {
                "W": 1.0,
                "C": 1.0,
                "M": 24.68,
                "H": 24.68,
                "Y": 1.0,
                "F": 1.0,
                "Q": 1.0,
                "N": 13.34,
                "I": 1.0,
                "R": 1.0,
                "D": 1.0,
                "P": 1.0,
                "T": -14.03,
                "K": 1.0,
                "E": 1.0,
                "V": -7.49,
                "S": 1.0,
                "G": -9.37,
                "A": -14.03,
                "L": 13.34,
            },
            "C": {
                "W": 24.68,
                "C": 1.0,
                "M": 33.6,
                "H": 33.6,
                "Y": 1.0,
                "F": 1.0,
                "Q": -6.54,
                "N": 1.0,
                "I": 1.0,
                "R": 1.0,
                "D": 20.26,
                "P": 20.26,
                "T": 33.6,
                "K": 1.0,
                "E": 1.0,
                "V": -6.54,
                "S": 1.0,
                "G": 1.0,
                "A": 1.0,
                "L": 20.26,
            },
            "M": {
                "W": 1.0,
                "C": 1.0,
                "M": -1.88,
                "H": 58.28,
                "Y": 24.68,
                "F": 1.0,
                "Q": -6.54,
                "N": 1.0,
                "I": 1.0,
                "R": -6.54,
                "D": 1.0,
                "P": 44.94,
                "T": -1.88,
                "K": 1.0,
                "E": 1.0,
                "V": 1.0,
                "S": 44.94,
                "G": 1.0,
                "A": 13.34,
                "L": 1.0,
            },
            "H": {
                "W": -1.88,
                "C": 1.0,
                "M": 1.0,
                "H": 1.0,
                "Y": 44.94,
                "F": -9.37,
                "Q": 1.0,
                "N": 24.68,
                "I": 44.94,
                "R": 1.0,
                "D": 1.0,
                "P": -1.88,
                "T": -6.54,
                "K": 24.68,
                "E": 1.0,
                "V": 1.0,
                "S": 1.0,
                "G": -9.37,
                "A": 1.0,
                "L": 1.0,
            },
            "Y": {
                "W": -9.37,
                "C": 1.0,
                "M": 44.94,
                "H": 13.34,
                "Y": 13.34,
                "F": 1.0,
                "Q": 1.0,
                "N": 1.0,
                "I": 1.0,
                "R": -15.91,
                "D": 24.68,
                "P": 13.34,
                "T": -7.49,
                "K": 1.0,
                "E": -6.54,
                "V": 1.0,
                "S": 1.0,
                "G": -7.49,
                "A": 24.68,
                "L": 1.0,
            },
            "F": {
                "W": 1.0,
                "C": 1.0,
                "M": 1.0,
                "H": 1.0,
                "Y": 33.6,
                "F": 1.0,
                "Q": 1.0,
                "N": 1.0,
                "I": 1.0,
                "R": 1.0,
                "D": 13.34,
                "P": 20.26,
                "T": 1.0,
                "K": -14.03,
                "E": 1.0,
                "V": 1.0,
                "S": 1.0,
                "G": 1.0,
                "A": 1.0,
                "L": 1.0,
            },
            "Q": {
                "W": 1.0,
                "C": -6.54,
                "M": 1.0,
                "H": 1.0,
                "Y": -6.54,
                "F": -6.54,
                "Q": 20.26,
                "N": 1.0,
                "I": 1.0,
                "R": 1.0,
                "D": 20.20,
                "P": 20.26,
                "T": 1.0,
                "K": 1.0,
                "E": 20.26,
                "V": -6.54,
                "S": 44.94,
                "G": 1.0,
                "A": 1.0,
                "L": 1.0,
            },
            "N": {
                "W": -9.37,
                "C": -1.88,
                "M": 1.0,
                "H": 1.0,
                "Y": 1.0,
                "F": -14.03,
                "Q": -6.54,
                "N": 1.0,
                "I": 44.94,
                "R": 1.0,
                "D": 1.0,
                "P": -1.88,
                "T": -7.49,
                "K": 24.68,
                "E": 1.0,
                "V": 1.0,
                "S": 1.0,
                "G": -14.03,
                "A": 1.0,
                "L": 1.0,
            },
            "I": {
                "W": 1.0,
                "C": 1.0,
                "M": 1.0,
                "H": 13.34,
                "Y": 1.0,
                "F": 1.0,
                "Q": 1.0,
                "N": 1.0,
                "I": 1.0,
                "R": 1.0,
                "D": 1.0,
                "P": -1.88,
                "T": 1.0,
                "K": -7.49,
                "E": 44.94,
                "V": -7.49,
                "S": 1.0,
                "G": 1.0,
                "A": 1.0,
                "L": 20.26,
            },
            "R": {
                "W": 58.28,
                "C": 1.0,
                "M": 1.0,
                "H": 20.26,
                "Y": -6.54,
                "F": 1.0,
                "Q": 20.26,
                "N": 13.34,
                "I": 1.0,
                "R": 58.28,
                "D": 1.0,
                "P": 20.26,
                "T": 1.0,
                "K": 1.0,
                "E": 1.0,
                "V": 1.0,
                "S": 44.94,
                "G": -7.49,
                "A": 1.0,
                "L": 1.0,
            },
            "D": {
                "W": 1.0,
                "C": 1.0,
                "M": 1.0,
                "H": 1.0,
                "Y": 1.0,
                "F": -6.54,
                "Q": 1.0,
                "N": 1.0,
                "I": 1.0,
                "R": -6.54,
                "D": 1.0,
                "P": 1.0,
                "T": -14.03,
                "K": -7.49,
                "E": 1.0,
                "V": 1.0,
                "S": 20.26,
                "G": 1.0,
                "A": 1.0,
                "L": 1.0,
            },
            "P": {
                "W": -1.88,
                "C": -6.54,
                "M": -6.54,
                "H": 1.0,
                "Y": 1.0,
                "F": 20.26,
                "Q": 20.26,
                "N": 1.0,
                "I": 1.0,
                "R": -6.54,
                "D": -6.54,
                "P": 20.26,
                "T": 1.0,
                "K": 1.0,
                "E": 18.38,
                "V": 20.26,
                "S": 20.26,
                "G": 1.0,
                "A": 20.26,
                "L": 1.0,
            },
            "T": {
                "W": -14.03,
                "C": 1.0,
                "M": 1.0,
                "H": 1.0,
                "Y": 1.0,
                "F": 13.34,
                "Q": -6.54,
                "N": -14.03,
                "I": 1.0,
                "R": 1.0,
                "D": 1.0,
                "P": 1.0,
                "T": 1.0,
                "K": 1.0,
                "E": 20.26,
                "V": 1.0,
                "S": 1.0,
                "G": -7.49,
                "A": 1.0,
                "L": 1.0,
            },
            "K": {
                "W": 1.0,
                "C": 1.0,
                "M": 33.6,
                "H": 1.0,
                "Y": 1.0,
                "F": 1.0,
                "Q": 24.68,
                "N": 1.0,
                "I": -7.49,
                "R": 33.6,
                "D": 1.0,
                "P": -6.54,
                "T": 1.0,
                "K": 1.0,
                "E": 1.0,
                "V": -7.49,
                "S": 1.0,
                "G": -7.49,
                "A": 1.0,
                "L": -7.49,
            },
            "E": {
                "W": -14.03,
                "C": 44.94,
                "M": 1.0,
                "H": -6.54,
                "Y": 1.0,
                "F": 1.0,
                "Q": 20.26,
                "N": 1.0,
                "I": 20.26,
                "R": 1.0,
                "D": 20.26,
                "P": 20.26,
                "T": 1.0,
                "K": 1.0,
                "E": 33.6,
                "V": 1.0,
                "S": 20.26,
                "G": 1.0,
                "A": 1.0,
                "L": 1.0,
            },
            "V": {
                "W": 1.0,
                "C": 1.0,
                "M": 1.0,
                "H": 1.0,
                "Y": -6.54,
                "F": 1.0,
                "Q": 1.0,
                "N": 1.0,
                "I": 1.0,
                "R": 1.0,
                "D": -14.03,
                "P": 20.26,
                "T": -7.49,
                "K": -1.88,
                "E": 1.0,
                "V": 1.0,
                "S": 1.0,
                "G": -7.49,
                "A": 1.0,
                "L": 1.0,
            },
            "S": {
                "W": 1.0,
                "C": 33.6,
                "M": 1.0,
                "H": 1.0,
                "Y": 1.0,
                "F": 1.0,
                "Q": 20.26,
                "N": 1.0,
                "I": 1.0,
                "R": 20.26,
                "D": 1.0,
                "P": 44.94,
                "T": 1.0,
                "K": 1.0,
                "E": 20.26,
                "V": 1.0,
                "S": 20.26,
                "G": 1.0,
                "A": 1.0,
                "L": 1.0,
            },
            "G": {
                "W": 13.34,
                "C": 1.0,
                "M": 1.0,
                "H": 1.0,
                "Y": -7.49,
                "F": 1.0,
                "Q": 1.0,
                "N": -7.49,
                "I": -7.49,
                "R": 1.0,
                "D": 1.0,
                "P": 1.0,
                "T": -7.49,
                "K": -7.49,
                "E": -6.54,
                "V": 1.0,
                "S": 1.0,
                "G": 13.34,
                "A": -7.49,
                "L": 1.0,
            },
            "A": {
                "W": 1.0,
                "C": 44.94,
                "M": 1.0,
                "H": -7.49,
                "Y": 1.0,
                "F": 1.0,
                "Q": 1.0,
                "N": 1.0,
                "I": 1.0,
                "R": 1.0,
                "D": -7.49,
                "P": 20.26,
                "T": 1.0,
                "K": 1.0,
                "E": 1.0,
                "V": 1.0,
                "S": 1.0,
                "G": 1.0,
                "A": 1.0,
                "L": 1.0,
            },
            "L": {
                "W": 24.68,
                "C": 1.0,
                "M": 1.0,
                "H": 1.0,
                "Y": 1.0,
                "F": 1.0,
                "Q": 33.6,
                "N": 1.0,
                "I": 1.0,
                "R": 20.26,
                "D": 1.0,
                "P": 20.26,
                "T": 1.0,
                "K": -7.49,
                "E": 1.0,
                "V": 1.0,
                "S": 1.0,
                "G": 1.0,
                "A": 1.0,
                "L": 1.0,
            },
        }

        # --- cys20/30; includes normalization in case a peptide is shorter than 30 aa
        if self.length < 20:
            self.cys30 = self.sequence.count("C") / float(self.length)
            self.cys20 = self.sequence.count("C") / float(self.length)
            self.cys_ser30 = (
                self.sequence.count("C") + self.sequence.count("S")
            ) / float(self.length)
            self.cys_ser20 = (
                self.sequence.count("C") + self.sequence.count("S")
            ) / float(self.length)

        elif 20 <= self.length < 30:
            self.cys30 = self.sequence.count("C") / float(self.length)
            self.cys_ser30 = (
                self.sequence.count("C") + self.sequence.count("S")
            ) / float(self.length)
            self.cys20 = max(
                [
                    self.sequence[rng : rng + 20].count("C")
                    / float(len(self.sequence[rng : rng + 20]))
                    for rng in range(0, self.length, 1)
                    if len(self.sequence[rng : rng + 20]) == 20
                ]
            )
            self.cys_ser20 = max(
                [
                    (
                        self.sequence[rng : rng + 20].count("C")
                        + self.sequence[rng : rng + 20].count("S")
                    )
                    / float(len(self.sequence[rng : rng + 20]))
                    for rng in range(0, self.length, 1)
                    if len(self.sequence[rng : rng + 20]) == 20
                ]
            )

        else:
            self.cys30 = max(
                [
                    self.sequence[rng : rng + 30].count("C")
                    / float(len(self.sequence[rng : rng + 30]))
                    for rng in range(0, self.length, 1)
                    if len(self.sequence[rng : rng + 30]) == 30
                ]
            )
            self.cys20 = max(
                [
                    self.sequence[rng : rng + 20].count("C")
                    / float(len(self.sequence[rng : rng + 20]))
                    for rng in range(0, self.length, 1)
                    if len(self.sequence[rng : rng + 20]) == 20
                ]
            )
            self.cys_ser30 = max(
                [
                    (
                        self.sequence[rng : rng + 30].count("C")
                        + self.sequence[rng : rng + 30].count("S")
                    )
                    / float(len(self.sequence[rng : rng + 30]))
                    for rng in range(0, self.length, 1)
                    if len(self.sequence[rng : rng + 30]) == 30
                ]
            )
            self.cys_ser20 = max(
                [
                    (
                        self.sequence[rng : rng + 20].count("C")
                        + self.sequence[rng : rng + 20].count("S")
                    )
                    / float(len(self.sequence[rng : rng + 20]))
                    for rng in range(0, self.length, 1)
                    if len(self.sequence[rng : rng + 20]) == 20
                ]
            )

        # --- aafreq
        for aa in self.aafreq:
            self.aafreq[aa] = self.sequence.count(aa) / float(self.length)

        # --- clfreq
        for aas in self.clfreq:
            self.clfreq[aas] = sum([self.sequence.count(aa) for aa in aas]) / float(
                self.length
            )

        # --- charge & avgcharge
        self.charge = sum(
            [charge_dict[aa] for aa in self.sequence if aa in self.aafreq]
        )
        self.avgcharge = self.charge / float(self.length)

        # --- avghydrop
        self.avghydrop = sum(
            [hydrophobicity_dict[aa] for aa in self.sequence if aa in self.aafreq]
        ) / float(self.length)

        # --- k-tuplet entropy
        es = []
        for rng in range(self.length):
            s = 0.0
            seqtmp = self.sequence[rng : rng + 10]
            if 1:  # len(seqtmp)==10:
                for i in range(1, len(seqtmp)):
                    with contextlib.suppress(KeyError):
                        s -= (
                            self.aafreq[seqtmp[i]]
                            * self.aafreq[seqtmp[i]]
                            * self.aafreq[seqtmp[i - 1]]
                            * log2(self.aafreq[seqtmp[i]] * self.aafreq[seqtmp[i - 1]])
                        )
            es.append(s)
        self.entropy = max(es)

        # --- entropy ratio
        es = []
        for rng in range(self.length):
            s = 0.0
            seqtmp = self.sequence[rng : rng + 10]
            if len(seqtmp) < 5:
                continue
            aatmp = set(seqtmp)
            stotal = -sum(
                [
                    self.aafreq[i] * log2(self.aafreq[i])
                    for i in aatmp
                    if i in self.aafreq
                ]
            )
            stemp = -sum(
                [
                    (seqtmp.count(i) / 10.0) * log2(seqtmp.count(i) / 10.0)
                    for i in aatmp
                    if i in self.aafreq
                ]
            )
            try:
                es.append(stemp / stotal)
            except ZeroDivisionError:
                es.append(2.0)
        try:
            self.entropyratio = min(es)
        except ValueError:
            print(f"Error with min(Es) on smORF with sequence {self.sequence}")

        # --- boman index
        self.boman = sum(
            [boman_dict[aa] for aa in self.sequence if aa in self.aafreq]
        ) / float(self.length)

        # --- inestability
        s =0
        for i in range(self.length - 1):
            try:
                s+=instability_dict[self.sequence[i]][self.sequence[i + 1]]
            except KeyError:
                s+=0
        self.instability = (10 / self.length) * s

        # --- aliphatic index
        # count aliphatic residues
        ala = self.sequence.count("A") / self.length
        val = self.sequence.count("V") / self.length
        leu = self.sequence.count("L") / self.length
        ile = self.sequence.count("I") / self.length
        # support unknown Leu/Ile residues
        xle = self.sequence.count("J") / self.length
        # return aliphatic index
        self.aliphatic = (ala + 2.9 * val + 3.9 * (leu + ile + xle)) * 100

        return self  # for use in list comprehension

    def make_list(self):
        aalist = [
            "A",
            "R",
            "N",
            "D",
            "C",
            "E",
            "Q",
            "G",
            "H",
            "I",
            "L",
            "K",
            "M",
            "F",
            "P",
            "S",
            "T",
            "W",
            "Y",
            "V",
        ]
        cllist = ["RHK", "DE", "STNQ", "CGP", "AVIL", "MFYW"]

        l = []
        l += [self.aafreq[aa] for aa in aalist]
        l += [self.clfreq[aa] for aa in cllist]
        l += [
            self.cys30,
            self.cys20,
            self.charge,
            self.avgcharge,
            self.avghydrop,
            self.entropy,
            self.entropyratio,
        ]

        return l

    def get_features(self):
        feature_list = [
            self.sequence,
            self.cys30,
            self.cys20,
            self.cys_ser30,
            self.cys_ser20,
            self.charge,
            self.avgcharge,
            self.avghydrop,
            self.length,
            self.entropy,
            self.entropyratio,
        ]
        amino_acids = self.aafreq.keys()
        clusters = self.clfreq.keys()
        feature_list += [self.aafreq[aa] for aa in amino_acids]
        feature_list += [self.clfreq[cl] for cl in clusters]
        return feature_list
