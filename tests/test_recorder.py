from abSENSE.recorder import File_Recorder
import pytest
from io import StringIO
import itertools


@pytest.fixture()
def faked_recorder(mocker):
    mocker.patch('abSENSE.recorder.os.makedirs')
    recorder = File_Recorder('.', ['s1', 's2', 's3'], predict_all=True)
    # replace files with stringIO for easier testing
    files = {}
    def new_stringIO(*args, **kwargs):
        result = StringIO()
        files[args[0]] = result
        return result
    mocker.patch('abSENSE.recorder.open', side_effect=new_stringIO)
    return recorder, files


def test_open(faked_recorder):
    rec, files = faked_recorder

    assert len(files) == 0

    with rec.open() as _:
        assert len(files) == 6
        for file in files.values():
            assert file.closed is False

    for file in files.values():
        assert file.closed


def test_write_headers(faked_recorder):
    rec, files = faked_recorder
    with rec.open() as recorder:
        recorder.write_headers()

        header = 'Gene\ts1\ts2\ts3\n'
        assert files['./predicted_bitscores.tsv'].getvalue() == (
            "# maximum likelihood bitscore predictions "
            "for each tested gene in each species\n"
            + header
        )
        assert files['./99PI_lower_prediction.tsv'].getvalue() == (
            "# the lower bound of the 99% bitscore prediction "
            "interval for each tested gene in each species\n"
            + header
        )
        assert files['./99PI_high_prediction.tsv'].getvalue() == (
            "# the upper bound of the 99% bitscore prediction "
            "interval for each tested gene in each species\n"
            + header
        )
        assert files['./failure_probabilities.tsv'].getvalue() == (
            "# the probability of a homolog being undetected "
            "at the specified significance threshold (see run info file) in each "
            "tested gene in each species\n"
            + header
        )
        assert files['./parameters.tsv'].getvalue() == (
            "# the best-fit parameters "
            "(performed using only bitscores from species not omitted from the fit; "
            "see run info file) for a and b for each gene\n"
            "Gene\ta\tb\n"  # not species header
        )
        assert files['./run_info.txt'].getvalue() == ''


@pytest.mark.parametrize('genelenfile', [None, 'gene len file'])
@pytest.mark.parametrize('dblenfile', [None, 'db len file'])
@pytest.mark.parametrize('pred_species', [[], ['a', 'b', 'c'], ['a']])
def test_write_info(faked_recorder, mocker, genelenfile, dblenfile, pred_species):
    rec, files = faked_recorder

    args = mocker.MagicMock()
    args.scorefile = 'score file'
    args.distfile = 'distance file'
    args.genelenfile = genelenfile
    args.dblenfile = dblenfile
    args.Eval = 'e value'

    with rec.open() as recorder:
        recorder.write_info('start time',
                            args,
                            'default gene length',
                            'default db size',
                            pred_species=pred_species,
                            )

        for key, file in files.items():
            if key != './run_info.txt':
                assert file.getvalue() == ''
                continue

            genelen = "Gene length (for all genes): default gene length (default)\n"
            if genelenfile is not None:
                genelen = f"Gene length file: {genelenfile}\n"

            dblen = "Database length (for all species): default db size (default)\n"
            if dblenfile is not None:
                dblen = f"Database length file: {dblenfile}\n"

            pred = "Species used in fit: All (default)\n"
            if pred_species != []:
                pred = f"Species used in fit: {' '.join(pred_species)}\n"

            assert file.getvalue() == (
                "abSENSE analysis run on start time\n"
                "Input bitscore file: score file\n"
                "Input distance file: distance file\n"
                f"{genelen}"
                f"{dblen}"
                f"{pred}"
                "E-value threshold: e value\n"
            )


def test_write_gene(faked_recorder):
    rec, files = faked_recorder
    with rec.open() as recorder:
        recorder.write_gene('my gene name')
        for key, file in files.items():
            if key == './run_info.txt':
                assert file.getvalue() == ''
                continue

            assert file.getvalue() == 'my gene name'


def test_analysis_error(faked_recorder):
    rec, files = faked_recorder
    with rec.open() as recorder:
        recorder.analysis_error()
        string = 'analysis_error'
        for key, file in files.items():
            if key == './run_info.txt':
                assert file.getvalue() == ''
                continue

            if key == './parameters.tsv':
                assert file.getvalue() == f'\t{string}\t{string}\n'
                continue

            assert file.getvalue() == f'\t{string}\t{string}\t{string}\n'


def test_not_enough_data(faked_recorder):
    rec, files = faked_recorder
    with rec.open() as recorder:
        recorder.not_enough_data()
        string = 'not_enough_data'
        for key, file in files.items():
            if key == './run_info.txt':
                assert file.getvalue() == ''
                continue

            if key == './parameters.tsv':
                assert file.getvalue() == f'\t{string}\t{string}\n'
                continue

            assert file.getvalue() == f'\t{string}\t{string}\t{string}\n'


def old_write_result(mloutputfile, highboundoutputfile, lowboundoutputfile,
                     pvaloutputfile, outputfileparams, predall, is_considered,
                     is_ambiguous, realscores, a, b, predictions,
                     highs, lows, pvals):
    for j in range(0, len(predictions)):
        prediction = predictions[j]
        lowprediction, highprediction, pval = round(lows[j], 2), round(highs[j], 2), round(pvals[j], 2)
        if (not is_considered[j] and not is_ambiguous[j]):
            mloutputfile.write(f"\t{prediction}")
            highboundoutputfile.write(f"\t{highprediction}")
            lowboundoutputfile.write(f"\t{lowprediction}")
            pvaloutputfile.write(f"\t{pval}")
        elif is_considered[j]:
            if predall:
                realscore = realscores[j]
                mloutputfile.write(
                    f"\t{prediction}(Ortholog_detected:{realscore})"
                )
                highboundoutputfile.write(
                    f"\t{highprediction}(Ortholog_detected)"
                )
                lowboundoutputfile.write(
                    f"\t{lowprediction}(Ortholog_detected)"
                )
                pvaloutputfile.write(
                    f"\t{pval}(Ortholog_detected)"
                )
            else:
                mloutputfile.write("\tOrtholog_detected")
                highboundoutputfile.write("\tOrtholog_detected")
                lowboundoutputfile.write("\tOrtholog_detected")
                pvaloutputfile.write("\tOrtholog_detected")
        elif is_ambiguous[j]:
            if predall:
                mloutputfile.write(
                    f"\t{prediction}"
                    "(Homolog_detected(orthology_ambiguous))"
                )
                highboundoutputfile.write(
                    f"\t{highprediction}"
                    "(Homolog_detected(orthology_ambiguous))"
                )
                lowboundoutputfile.write(
                    f"\t{lowprediction}"
                    "(Homolog_detected(orthology_ambiguous))"
                )
                pvaloutputfile.write(
                    f"\t{pval}"
                    "(Homolog_detected(orthology_ambiguous))"
                )
            else:
                mloutputfile.write(
                    "\tHomolog_detected(orthology_ambiguous)"
                )
                highboundoutputfile.write(
                    "\tHomolog_detected(orthology_ambiguous)"
                )
                lowboundoutputfile.write(
                    "\tHomolog_detected(orthology_ambiguous)"
                )
                pvaloutputfile.write(
                    "\tHomolog_detected(orthology_ambiguous)"
                )

    mloutputfile.write("\n")
    highboundoutputfile.write("\n")
    lowboundoutputfile.write("\n")
    pvaloutputfile.write("\n")
    outputfileparams.write("\t" + str(a) + "\t" + str(b))
    outputfileparams.write("\n")


all_true_false = list(itertools.combinations_with_replacement([True, False], 3))
numbers = [[0.1, 0.2, 0.3], [-0.12345, .12345, .2468]]


@pytest.mark.parametrize('predall', [False, True])
@pytest.mark.parametrize('predictions', numbers)
@pytest.mark.parametrize('highs', numbers)
@pytest.mark.parametrize('lows', numbers)
@pytest.mark.parametrize('pvals', numbers)
@pytest.mark.parametrize('is_ambiguous', all_true_false)
@pytest.mark.parametrize('is_considered', all_true_false)
def test_write_result(faked_recorder, predall, highs, predictions,
                      lows, pvals, is_ambiguous, is_considered):
    # get old results
    a = 'a value'
    b = 'b value'
    realscores = [0, 1, 2]
    test_files = [StringIO() for _ in range(5)]
    old_write_result(test_files[0], test_files[1], test_files[2],
                     test_files[3], test_files[4], predall,
                     is_considered, is_ambiguous, realscores, a, b,
                     predictions, highs, lows, pvals)

    rec, files = faked_recorder
    rec.predict_all = predall
    with rec.open() as recorder:
        for prediction, high, low, pval, realscore, considered, ambiguous in zip(
                predictions, highs, lows, pvals, realscores, is_considered, is_ambiguous):
            recorder.write_result(prediction, high, low, pval, realscore, considered, ambiguous)
        recorder.write_params(a, b)
        recorder.finalize_row()

        assert files['./predicted_bitscores.tsv'].getvalue() == \
            test_files[0].getvalue()
        assert files['./99PI_high_prediction.tsv'].getvalue() == test_files[1].getvalue()
        assert files['./99PI_lower_prediction.tsv'].getvalue() == test_files[2].getvalue()
        assert files['./failure_probabilities.tsv'].getvalue() == test_files[3].getvalue()
        assert files['./parameters.tsv'].getvalue() == test_files[4].getvalue()

    for file in test_files:
        file.close()
