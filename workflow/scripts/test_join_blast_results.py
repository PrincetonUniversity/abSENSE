from __future__ import annotations

from io import StringIO

import pytest
from blast_records import ERRORS, BlastRecords


@pytest.fixture
def records(mocker):
    # make open just return the object so StringIO objects can be passed in
    mocker.patch("blast_records.open", side_effect=lambda x, **kwargs: x)
    # override error strings as most will be 0/NA and hard to determine cause
    ERRORS["high_e_val"] = "HEV"
    ERRORS["unset"] = "U"
    ERRORS["mismatch"] = "MM"
    ERRORS["missing"] = "MS"
    return BlastRecords("s1 s2".split(), threshold=0.1)


def perform_matching(db, self_match, forward, reciprocal):
    db.read_self_match(self_match)
    db.read_forward_match(forward, "s2")
    db.read_reciprocal_match(reciprocal, "s2")

    output = StringIO()
    db.report(output)
    return output.getvalue()


def test_matching(records):
    """Test handling when genes are found, match and e_val passes."""
    self_match = StringIO("gene1 gene1 100 0.0")
    forward = StringIO("gene1 geneA 101 0.0")
    recip = StringIO("geneA gene1 102 0.0")

    output = perform_matching(records, self_match, forward, recip).split("\n")

    assert output[0] == "gene\ts1\ts2"
    assert output[1] == "gene1\t100.0\t101.0"


def test_matches_high_eval(records):
    """Test handling when genes are found, match and e_val fails."""
    self_match = StringIO("gene1 gene1 100 0.0")
    forward = StringIO("gene1 geneA 101 0.5")
    recip = StringIO("geneA gene1 102 0.0")

    output = perform_matching(records, self_match, forward, recip).split("\n")

    assert output[1] == "\t".join("gene1 100.0 HEV".split())


def test_mismatches_high_eval(records):
    """Test handling when genes are found, mismatch and e_val fails."""
    self_match = StringIO("gene1 gene1 100 0.0")
    forward = StringIO("gene1 geneA 101 0.5")
    recip = StringIO("geneA gene2 102 0.0")

    output = perform_matching(records, self_match, forward, recip).split("\n")

    assert output[1] == "\t".join("gene1 100.0 HEV".split())


def test_no_matches_forward(records):
    """Test handling when genes aren't found in reciprocal search."""
    self_match = StringIO("gene1 gene1 100 0.0")
    forward = StringIO("")
    recip = StringIO("")

    output = perform_matching(records, self_match, forward, recip).split("\n")

    assert output[1] == "gene1\t100.0\tMS"


def test_no_matches_reciprocal(records):
    """Test handling when genes aren't found in reciprocal search."""
    self_match = StringIO("gene1 gene1 100 0.0")
    forward = StringIO("gene1 geneA 101 0.0")
    recip = StringIO("")

    output = perform_matching(records, self_match, forward, recip).split("\n")

    assert output[1] == "gene1\t100.0\tU"


def test_mismatches(records):
    """Test handling when genes mismatch in reciprocal search."""
    self_match = StringIO("gene1 gene1 102 0.0")
    forward = StringIO("gene1 geneA 101 0.0")
    recip = StringIO("geneA gene3 102 0.0")

    output = perform_matching(records, self_match, forward, recip).split("\n")

    assert output[1] == "gene1\t102.0\tMM"


def test_unmatched_reciprocal(records):
    """Test handling when reciprocal gene is not found in forward search."""
    self_match = StringIO("gene1 gene1 102 0.0")
    forward = StringIO("gene1 geneA 101 0.0")
    recip = StringIO("geneA gene3 102 0.0\ngene4 gene3 103 0.0")

    with pytest.raises(ValueError) as error:
        perform_matching(records, self_match, forward, recip).split("\n")
    assert str(error.value) == (
        "Unable to match reciprocal search to genes preivously found: "
        "BlastRecord(query='gene4', match='gene3', bitscore=103.0, e_val=0.0, reported_value=None)"
    )


def test_many_forward_match_recip(records):
    """Test handling when many forward genes map to same reverse."""
    self_match = StringIO("gene1 gene1 102 0.0\ngene2 gene2 103 0.0\n")
    forward = StringIO("gene1 geneA 104 0.0\ngene2 geneA 105 0.0\n")
    recip = StringIO("geneA gene2 102 0.0")

    output = perform_matching(records, self_match, forward, recip).split("\n")

    # gene1 is left unset because gene2 was a match
    assert output[1] == "gene1\t102.0\tU"
    assert output[2] == "gene2\t103.0\t105.0"


def test_many_forward_mismatch_recip(records):
    """Test handling when many forward genes map to same reverse."""
    self_match = StringIO("gene1 gene1 102 0.0\ngene2 gene2 103 0.0\n")
    forward = StringIO("gene1 geneA 104 0.0\ngene2 geneA 105 0.0\n")
    recip = StringIO("geneA gene3 102 0.0")

    output = perform_matching(records, self_match, forward, recip).split("\n")

    # gene1 the first found and so is a mismatch
    assert output[1] == "gene1\t102.0\tMM"
    # after first mismatch, break and leave gene2 unset
    assert output[2] == "gene2\t103.0\tU"


def test_many_recip_match_forard(records):
    """Test handling when many reciprocal genes map to same forward."""
    self_match = StringIO("gene1 gene1 102 0.0\ngene2 gene2 103 0.0\n")
    forward = StringIO("gene1 geneA 104 0.0\ngene2 geneB 105 0.0\n")
    recip = StringIO("geneA gene2 102 0.0\ngeneB gene2 102 0.0\n")

    output = perform_matching(records, self_match, forward, recip).split("\n")

    assert output[1] == "gene1\t102.0\tMM"
    assert output[2] == "gene2\t103.0\t105.0"


def test_missing_forward(records):
    """Test handling when many reciprocal genes map to same forward."""
    self_match = StringIO("gene1 gene1 102 0.0\ngene2 gene2 103 0.0\n")
    forward = StringIO("gene1 geneA 104 0.0\ngene2 geneB 105 0.0\n")
    recip = StringIO("geneA gene2 102 0.0\ngeneB gene2 102 0.0\n")

    output = perform_matching(records, self_match, forward, recip).split("\n")

    assert output[1] == "gene1\t102.0\tMM"
    assert output[2] == "gene2\t103.0\t105.0"
