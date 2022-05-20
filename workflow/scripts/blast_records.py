"""Library for blast records and database."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Generator

# string constants for reporting, can change for testing
ERRORS = {"high_e_val": "0", "unset": "N/A", "mismatch": "N/A", "missing": "0"}


@dataclass
class BlastRecord:
    """Class for a single blast record."""

    query: str
    match: str
    bitscore: float
    e_val: float
    reported_value: str | None

    @classmethod
    def from_blast_line(cls, line: str, threshold: float) -> BlastRecord:
        """Produce a new blast record from the provided line."""
        query, match, bitscore, e_val = line.split()
        reported_value = None
        if float(e_val) > threshold:
            # set reported value to 0
            reported_value = ERRORS["high_e_val"]
        return cls(
            query=query,
            match=match,
            bitscore=float(bitscore),
            e_val=float(e_val),
            reported_value=reported_value,
        )

    @property
    def report(self) -> str:
        """Report the value of this record."""
        if self.reported_value is None:
            return ERRORS["unset"]
        return self.reported_value

    def found_match(self) -> None:
        """Tell the record a matching reciprocal record has been found."""
        if self.reported_value is None:
            self.reported_value = str(self.bitscore)

    def found_mismatch(self) -> None:
        """Tell the record a reciprocal search did not match."""
        if self.reported_value is None:
            self.reported_value = ERRORS["mismatch"]


class BlastRecords:
    """A database of blast records."""

    def __init__(self, species: list[str], threshold=0.001):
        self.species = species
        self.genes: list[str] = []
        self.database: dict[str, BlastRecord] = {}
        self.threshold = threshold

    def read_self_match(self, infile: str) -> None:
        """Given self_match blast results, record genes and fill database."""
        for record in self._read_file(infile):
            self.genes.append(record.query)
            record.found_match()  # always self match
            self._update_record(self.species[0], record.query, record)

    def read_forward_match(self, infile: str, species: str) -> None:
        """Given query=target and match=species, record all values."""
        for record in self._read_file(infile):
            self._update_record(species, record.query, record)

    def read_reciprocal_match(self, infile: str, species: str) -> None:
        """Validate database records based on reciprocal search."""
        for record in self._read_file(infile):
            # check if reciprocal match is present, this is fast
            if (
                record.match in self.database
                and species in self.database[record.match]
                and self.database[record.match][species].match == record.query
            ):
                self.database[record.match][species].found_match()
                continue

            # need to find database record and report mismatch
            for species_db in self.database.values():
                if (
                    species in species_db  # for this gene
                    and species_db[species].match == record.query
                ):
                    species_db[species].found_mismatch()
                    break
            else:  # if not found
                raise ValueError(
                    "Unable to match reciprocal search to genes preivously found: "
                    f"{record}"
                )

    def _read_file(self, infile) -> Generator[BlastRecord, None, None]:
        """Open provided file and return a generator of blast records."""
        with open(infile, encoding="UTF-8") as records:
            for line in records:
                yield BlastRecord.from_blast_line(line, self.threshold)

    def _update_record(self, species: str, gene: str, record: BlastRecord) -> None:
        if gene not in self.database:
            # add first record
            self.database[gene] = {species: record}
            return

        if species in self.database[gene]:
            # record already exists
            if self.database[gene][species].e_val > record.e_val:
                # override if e_val is smaller
                self.database[gene][species] = record
            return

        self.database[gene][species] = record

    def report(self, outfile) -> None:
        """Produce tsv file of distances for abSENSE."""
        outfile.write("\t".join(["gene"] + self.species))
        for gene in self.genes:
            outfile.write(f"\n{gene}\t")
            outfile.write(
                "\t".join(
                    self.database[gene][species].report
                    if species in self.database[gene]
                    else ERRORS["missing"]
                    for species in self.species
                )
            )
