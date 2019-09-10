import os
from Bio.SeqUtils import MeltingTemp as mt

from Bio.SeqUtils import IsoelectricPoint as isp
from Bio.Seq import Seq, MutableSeq
from Bio.Alphabet import generic_rna
from Bio import SeqIO


from Bio import Alphabet
from Bio.Data import IUPACData

def molecular_weight(seq, seq_type=None, double_stranded=False, circular=False,
                     monoisotopic=False):
    """Calculate the molecular mass of DNA, RNA or protein sequences as float.
    Only unambiguous letters are allowed. Nucleotide sequences are assumed to
    have a 5' phosphate.
    Arguments:
     - seq: String or Biopython sequence object.
     - seq_type: The default (None) is to take the alphabet from the seq
       argument, or assume DNA if the seq argument is a string. Override this
       with a string 'DNA', 'RNA', or 'protein'.
     - double_stranded: Calculate the mass for the double stranded molecule?
     - circular: Is the molecule circular (has no ends)?
     - monoisotopic: Use the monoisotopic mass tables?
    Note that for backwards compatibility, if the seq argument is a string,
    or Seq object with a generic alphabet, and no seq_type is specified
    (i.e. left as None), then DNA is assumed.
    >>> print("%0.2f" % molecular_weight("AGC"))
    949.61
    >>> print("%0.2f" % molecular_weight(Seq("AGC")))
    949.61
    However, it is better to be explicit - for example with strings:
    >>> print("%0.2f" % molecular_weight("AGC", "DNA"))
    949.61
    >>> print("%0.2f" % molecular_weight("AGC", "RNA"))
    997.61
    >>> print("%0.2f" % molecular_weight("AGC", "protein"))
    249.29
    Or, with the sequence alphabet:
    >>> from Bio.Seq import Seq
    >>> from Bio.Alphabet import generic_dna, generic_rna, generic_protein
    >>> print("%0.2f" % molecular_weight(Seq("AGC", generic_dna)))
    949.61
    >>> print("%0.2f" % molecular_weight(Seq("AGC", generic_rna)))
    997.61
    >>> print("%0.2f" % molecular_weight(Seq("AGC", generic_protein)))
    249.29
    Also note that contradictory sequence alphabets and seq_type will also
    give an exception:
    >>> from Bio.Seq import Seq
    >>> from Bio.Alphabet import generic_dna
    >>> print("%0.2f" % molecular_weight(Seq("AGC", generic_dna), "RNA"))
    Traceback (most recent call last):
      ...
    ValueError: seq_type='RNA' contradicts DNA from seq alphabet
    """
    # Rewritten by Markus Piotrowski, 2014

    # Find the alphabet type
    tmp_type = ""
    if isinstance(seq, (Seq, MutableSeq)):
        base_alphabet = Alphabet._get_base_alphabet(seq.alphabet)
        if isinstance(base_alphabet, Alphabet.DNAAlphabet):
            tmp_type = "DNA"
        elif isinstance(base_alphabet, Alphabet.RNAAlphabet):
            tmp_type = "RNA"
        elif isinstance(base_alphabet, Alphabet.ProteinAlphabet):
            tmp_type = "protein"
        elif isinstance(base_alphabet, Alphabet.ThreeLetterProtein):
            tmp_type = "protein"
            # Convert to one-letter sequence. Have to use a string for seq1
            seq = Seq(seq1(str(seq)), alphabet=Alphabet.ProteinAlphabet())
        elif not isinstance(base_alphabet, Alphabet.Alphabet):
            raise TypeError("%s is not a valid alphabet for mass calculations"
                            % base_alphabet)
        else:
            tmp_type = "DNA"  # backward compatibity
        if seq_type and tmp_type and tmp_type != seq_type:
            raise ValueError("seq_type=%r contradicts %s from seq alphabet"
                             % (seq_type, tmp_type))
        seq_type = tmp_type
    elif isinstance(seq, str):
        if seq_type is None:
            seq_type = "DNA"  # backward compatibity
    else:
        raise TypeError("Expected a string or Seq object, not seq=%r" % seq)

    seq = "".join(str(seq).split()).upper()  # Do the minimum formatting

    if seq_type == "DNA":
        if monoisotopic:
            weight_table = IUPACData.monoisotopic_unambiguous_dna_weights
        else:
            weight_table = IUPACData.unambiguous_dna_weights
    elif seq_type == "RNA":
        if monoisotopic:
            weight_table = IUPACData.monoisotopic_unambiguous_rna_weights
        else:
            weight_table = IUPACData.unambiguous_rna_weights
    elif seq_type == "protein":
        if monoisotopic:
            weight_table = IUPACData.monoisotopic_protein_weights
        else:
            weight_table = IUPACData.protein_weights
    else:
        raise ValueError("Allowed seq_types are DNA, RNA or protein, not %r"
                         % seq_type)

    if monoisotopic:
        water = 18.010565
    else:
        water = 18.0153

    try:
        weight = sum(weight_table[x] for x in seq) - (len(seq) - 1) * water
        if circular:
            weight -= water
    except KeyError as e:
        #raise ValueError("%s is not a valid unambiguous letter for %s"
        #                % (e, seq_type))
        return -1488
        #break

    if seq_type in ("DNA", "RNA") and double_stranded:
        seq = str(Seq(seq).complement())
        weight += sum(weight_table[x] for x in seq) - (len(seq) - 1) * water
        if circular:
            weight -= water
    elif seq_type == "protein" and double_stranded:
        raise ValueError("double-stranded proteins await their discovery")

    return weight

def nucl_cont(dic):
    
    G = float(dic.count('G'))
    C = float(dic.count('C'))
    A = float(dic.count('A'))
    L = len(dic)
    
    return G,C,A,L


def temperatures(dic):
    Tw = round(mt.Tm_Wallace(dic, strict=False), 2)
    Tgc = round(mt.Tm_GC(dic, strict=False), 2)
    Tnn = round(mt.Tm_NN(dic, strict=False), 2)
    return Tw, Tgc, Tnn

def weights(dic,G,C,A,L):
    mw = molecular_weight(dic)
    if (mw == -1488):
        mw = 0.00
    mass = (A * 329.2) + ((L - A - G - C) * 306.2) + (C * 305.2) + (G * 345.2) + 159.0
    return round(mw, 2), round(mass, 2)
