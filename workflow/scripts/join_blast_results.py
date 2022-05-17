import sys
from dataclasses import dataclass
from pprint import pprint
from typing import Optional


# string constants for reporting, can change for testing
ERRORS = {
    'high_e_val': '0',
    'unset': 'ERROR',
    'mismatch': 'N/A',
    'missing': 'N/A'
}


@dataclass
class BlastRecord():
    query: str
    match: str
    bitscore: float
    e_val: float
    reported_value: Optional[str]

    @classmethod
    def from_blast_line(cls, line: str, threshold: float):
        query, match, bitscore, e_val = line.split()
        reported_value = None
        if float(e_val) > threshold:
            # set reported value to 0
            reported_value = ERRORS['high_e_val']
        return cls(
            query=query,
            match=match,
            bitscore=float(bitscore),
            e_val=float(e_val),
            reported_value=reported_value,
        )

    @property
    def report(self):
        if self.reported_value is None:
            return ERRORS['unset']
        return self.reported_value

    def found_match(self):
        if self.reported_value is None:
            self.reported_value = str(self.bitscore)

    def found_mismatch(self):
        if self.reported_value is None:
            self.reported_value = ERRORS['mismatch']

class BlastRecords():
    def __init__(self, species, threshold=0.001):
        self.species = species
        self.genes = []
        self.db = {}
        self.threshold = threshold

    def read_self_match(self, infile: str):
        """Given self_match blast results, record genes and fill db."""
        for record in self._read_file(infile):
            self.genes.append(record.query)
            record.found_match()  # always self match
            self._update_record(self.species[0], record.query, record)

    def read_forward_match(self, infile: str, species: str):
        """Given query=target and match=species, record all values."""
        for record in self._read_file(infile):
            self._update_record(species, record.query, record)

    def read_reciprocal_match(self, infile: str, species: str):
        """Validate db records based on reciprocal search."""
        for record in self._read_file(infile):
            # check if reciprocal match is present, this is fast
            if (record.match in self.db
                    and species in self.db[record.match]
                    and self.db[record.match][species].match == record.query):
                self.db[record.match][species].found_match()
                continue

            # need to find db record and report mismatch
            for species_db in self.db.values():
                if (species in species_db  # for this gene
                        and species_db[species].match == record.query):
                    species_db[species].found_mismatch()
                    break
            else:  # if not found
                raise ValueError('Unable to match reciprocal search to genes preivously found: '
                                 f'{record}')

    def _read_file(self, infile):
        with open(infile) as records:
            for line in records:
                yield BlastRecord.from_blast_line(line, self.threshold)

    def _update_record(self, species: str, gene: str, record: BlastRecord):
        if gene not in self.db:
            # add first record
            self.db[gene] = {species: record}
            return

        if species in self.db[gene]:
            # record already exists
            if self.db[gene][species].e_val > record.e_val:
                # override if e_val is smaller
                self.db[gene][species] = record
            return

        self.db[gene][species] = record

    def report(self, outfile):
        """Produce tsv file of distances for absense."""
        outfile.write('\t'.join(['gene'] + self.species))
        for gene in self.genes:
            outfile.write(f'\n{gene}\t')
            outfile.write('\t'.join(
                self.db[gene][species].report if species in self.db[gene] else ERRORS['missing']
                for species in self.species))


def main():
    if 'snakemake' in globals():
        self_match = snakemake.input['self_match']
        forwards = snakemake.input['forward']
        reciprocal = snakemake.input['reciprocal']
        species = snakemake.params['species'].split(',')
        e_val_threshold = snakemake.params['e_val_threshold']
        output = snakemake.output[0]
    else:
        base = 'scoring/blast/{query}/{database}.txt'
        species = ['S_cer', 'S_par', 'S_mik', 'S_kud', 'S_bay', 'S_cas',
                   'K_wal', 'A_gos', 'K_lac', 'A_nid', 'S_pom', 'Y_lip']
        self_match = base.format(query=species[0], database=species[0])
        forwards = [base.format(query=species[0], database=db) for db in species[1:]]
        reciprocal = [base.format(query=q, database=species[0]) for q in species[1:]]
        e_val_threshold = 0.001
        output = sys.stdout
        # raise ValueError('Only supported through snakemake')

    records = BlastRecords(species)
    records.read_self_match(self_match)

    for specie, forward, reciprocal in zip(species[1:], forwards, reciprocal):
        records.read_forward_match(forward, specie)
        records.read_reciprocal_match(reciprocal, specie)

    records.report(sys.stdout)

if __name__ == '__main__':
    main()
