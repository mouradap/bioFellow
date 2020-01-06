aminoacids = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
 'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

def translate(cds):
	codons = {
	
		"TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
		"TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
		"TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
		"TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
		"CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
		"CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
		"CAT": "H", 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
		'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
		'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
		'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
		'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
		'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
		'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
		'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
		'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
		'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
	}

	translation = ""

	length = int(len(cds)/3)

	for pos in range(length):
		codon = cds[pos*3:pos*3+3].upper()

		if codon in codons.keys():
			aminoacid = codons[codon]

			translation += aminoacid
		else:
			return "Invalid Codon"


	return translation

def revTranslation(aminoacidSeq):
	codons = {
	
		"TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
		"TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
		"TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
		"TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
		"CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
		"CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
		"CAT": "H", 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
		'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
		'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
		'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
		'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
		'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
		'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
		'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
		'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
		'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
	}


	revTrans = ""

	for pos in range(len(aminoacidSeq)):
		aa = aminoacidSeq[pos].upper()
		aaList = []
		for key in codons:
			aaList.append(codons[key])

		if aa not in aaList:
			return "Invalid Aminoacid"

		for key in codons:
			if codons[key] == aa:
				codon = key
				revTrans += codon
				break

	return revTrans

def mRNAtoDNA(mRNA):

	DNA = ""
	validChars = ['A', 'U', 'C', 'G']
	for char in mRNA.strip().upper():

		if char not in validChars:
			return "Invalid Codon"
		else:

			if char == "U":
				char = "T"
			DNA += char

	return DNA

def DNAtomRNA(cds):
	mRNA = ""
	validChars = ['A', 'T', 'C', 'G']
	for char in cds.strip().upper():

		if char not in validChars:
			return "Invalid Codon"
		else:

			if char == "T":
				char = "U"
			mRNA += char

	return mRNA


def orfFinder(sequence):
	possibleStarts = []
	terminators = ['TAA', 'TGA', 'TAG']
	possibleOrfs = []
	validChars = ['T', 'A', 'C', 'G']


	sequence = sequence.upper()

	for char in sequence.strip():
		if char not in validChars:
			return "Invalid Nucleotide"

	for i in range(len(sequence)):
		isStartCodon = sequence[i:i+3]
		if isStartCodon == "ATG":
			possibleStarts.append(i)

	for startCodon in possibleStarts:
		startedSeq = sequence[startCodon:]
		for i in range(int(len(startedSeq)/3)):
			isOrf = startedSeq[:i*3]
			if isOrf[-3:] in terminators:
				possibleOrfs.append([startCodon, (startCodon + i*3)])

	# print(possibleOrfs)

	foundOrfs = dict()

	for orf in enumerate(possibleOrfs):
		cds = sequence[orf[1][0]:orf[1][1]]
		# print(cds)
		foundOrfs[orf[0]] = [orf[1], cds]

	return foundOrfs