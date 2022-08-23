"""Main method for joining blast records called from snakemake."""

from blast_records import BlastRecords


def main() -> None:
    """Main method."""
    if "snakemake" in globals():
        self_match = snakemake.input["self_match"][0]  # noqa: F821
        forwards = snakemake.input["forward"]  # noqa: F821
        reciprocal = snakemake.input["reciprocal"]  # noqa: F821
        species = snakemake.params["species"].split(",")  # noqa: F821
        e_val_threshold = snakemake.params["e_val_threshold"]  # noqa: F821
        output = snakemake.output[0]  # noqa: F821
    else:
        raise ValueError("Only supported through snakemake")

    records = BlastRecords(species, threshold=e_val_threshold)
    records.read_self_match(self_match)

    for specie, forward, recip in zip(species[1:], forwards, reciprocal):
        records.read_forward_match(forward, specie)
        records.read_reciprocal_match(recip, specie)

    with open(output, "w", encoding="UTF-8") as outfile:
        records.report(outfile)


if __name__ == "__main__":
    main()
