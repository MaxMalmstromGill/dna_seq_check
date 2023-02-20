# A class for DNA sequence
class DnaSeq:
    def __init__(self, accession, seq):
        # Check if accession and sequence are not empty
        if not accession or not seq:
            raise ValueError("Accession and sequence cannot be empty")
        # Set the attributes
        self.accession = accession
        self.seq = seq

    # A method to return the length of the sequence
    def __len__(self):
        return len(self.seq)

    # A method to return a string representation of the object
    def __str__(self):
        return f"<DnaSeq accession='{self.accession}'>"

# A function to read DNA sequences from a file
def read_dna(file):
    # Open the file and split the sequences
    with open(file, 'r') as f:
        data = f.read().split('>')
        # Remove the empty first element
        data.pop(0)
        dna_seqs = []
        # Loop over the sequences
        for seq in data:
            # Split the sequence into the accession and the DNA string
            seq = seq.strip().split('\n')
            accession = seq[0]
            seq = ''.join(seq[1:])
            # Create a DnaSeq object and add it to the list
            dna_seq = DnaSeq(accession, seq)
            dna_seqs.append(dna_seq)
        # Return the list of DnaSeq objects
        return dna_seqs

# A function to check the exact overlap between two DNA sequences
def check_exact_overlap(seq1, seq2, min_length=10):
    """
    Detect exact overlaps between two DNA sequences.
    """
    # Get the lengths of the sequences
    len1, len2 = len(seq1), len(seq2)
    max_overlap = 0

    # Convert the DnaSeq objects to strings
    seq1_str, seq2_str = seq1.seq, seq2.seq

    # Check for overlaps of length at least min_length
    for i in range(min_length, min(len1, len2) + 1):
        if seq1_str[len1-i:] == seq2_str[:i]:
            max_overlap = i

    return max_overlap

# A function to find overlaps between DNA sequences
def overlaps(data, overlap_fn):
    results = {}
    # Loop over the sequences
    for seq1 in data:
        for seq2 in data:
            # Skip comparing a sequence with itself
            if seq1 == seq2:
                continue
            # Get the overlap length using the provided function
            overlap_len = overlap_fn(seq1, seq2)
            if overlap_len:
                # Add the result to the dictionary
                if seq1.accession not in results:
                    results[seq1.accession] = {}
                results[seq1.accession][seq2.accession] = overlap_len
    # Return the dictionary of overlap results
    return results



#
# Testing code. You should not change any code below here. To run the tests
# uncomment the last line in the file.
#
def test_class_DnaSeq():
    s1 = DnaSeq('s1', 'ACGT')
    s2 = DnaSeq('s2', 'ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTTAATCTTACAACCAGAACTCAAT')
    assert len(s1) == 4, 'Your length method (__len__) is not correct.'
    assert len(s2) == 70, 'Your length method (__len__) is not correct.'

    assert str(s1) == "<DnaSeq accession='s1'>", 'The __str__ method is not following the specification.'
    assert str(s2) == "<DnaSeq accession='s2'>", 'The __str__ method is not following the specification.'

    # The rest of this function is verifying that we are indeed raising an exception.
    status = 0
    try:
        s3 = DnaSeq('', 'ACGT')
    except ValueError:
        status += 1
    try:
        s3 = DnaSeq('s3', None)
    except ValueError:
        status += 1

    try:
        s3 = DnaSeq(None, '')
    except ValueError:
        status += 1
    if status != 3:
        raise Exception('class DnaSeq does not raise a ValueError '
                        'exception with initialised with empty '
                        'accession and sequence.')
    print('DnaSeq passed')


def test_reading():
    dna1 = read_dna('ex1.fa')
    assert len(dna1) == 6, 'The file "ex1.fa" has exactly 6 sequences, but your code does not return that.'
    assert list(map(lambda x: x.accession, dna1)) == [f's{i}' for i in range(6)], 'The accessions are not read correctly'

    dna2 = read_dna('ex2.fa')
    assert len(dna2) == 6, 'The file "ex2.fa" has exactly 6 sequences, but your code does not return that.'

    covid = read_dna('sars_cov_2.fa')
    assert len(covid[0].seq) == 29903, 'The length of the genome in "sars_cov_2.fa" is 29903, but your code does not return that.'

    print('read_dna passed')


def test_overlap():
   s0 = DnaSeq('s0', 'AAACCC')
   s1 = DnaSeq('s1', 'CCCGGG')
   s2 = DnaSeq('s2', 'TTATCC')
   s3 = DnaSeq('s3', 'CCAGGG')
   s4 = DnaSeq('s4', 'GGGGGGGGAAAGGGGG')
   s5 = DnaSeq('s5', 'AAATTTTTTTTTTTTTTTTTAT')

   data1 = [s0, s1, s2, s3]
   assert check_exact_overlap(s0, s1, 2) == 3
   assert check_exact_overlap(s0, s1) == 0
   assert check_exact_overlap(s0, s3, 2) == 2
   assert check_exact_overlap(s1, s2, 2) == 0
   assert check_exact_overlap(s2, s1, 2) == 2
   assert check_exact_overlap(s2, s3, 2) == 2
   assert check_exact_overlap(s4, s5, 1) == 0, 'Do not allow "internal" substrings to overlap. s4 and s5 does not have an overlap.'
   assert check_exact_overlap(s4, s5, 2) == 0
   assert check_exact_overlap(s4, s5, 3) == 0
   assert check_exact_overlap(s5, s2, 1) == 4

   res0 = overlaps(data1, lambda s1, s2: check_exact_overlap(s1, s2, 2))
   assert len(res0) == 2, 'You get the wrong number of overlaps'
   assert res0 == {'s0': {'s1': 3, 's3': 2}, 's2': {'s1': 2, 's3': 2}}

   dna_data = read_dna('ex1.fa')
   res1 = overlaps(dna_data, check_exact_overlap)
   assert len(res1) == 5
   for left, right in [('s0', 's1'), ('s1', 's2'), ('s2', 's3'), ('s3', 's4'),
                       ('s4', 's5')]:
       assert res1[left][right], f'Missing overlap of {left} and {right} (in that order)'
   print('overlap code passed')


def test_all():
    test_class_DnaSeq()
    test_reading()
    test_overlap()
    print('Yay, all good')

# Uncomment this to test everything:
test_all()
