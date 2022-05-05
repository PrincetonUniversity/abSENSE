from abSENSE.recorder import FileRecorder
from abSENSE.parameters import AbsenseParameters
import pytest
import pandas as pd
from io import StringIO
import itertools


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
def test_write_info(faked_recorder, genelenfile, dblenfile, pred_species):
    rec, files = faked_recorder

    include_only = None if pred_species == [] else ','.join(pred_species)
    dist = StringIO()
    dist.name = 'distance file'
    score = StringIO()
    score.name = 'score file'
    db_lens = None
    if dblenfile is not None:
        db_lens = StringIO()
        db_lens.name = dblenfile
    gen_lens = None
    if genelenfile is not None:
        gen_lens = StringIO()
        gen_lens.name = genelenfile
    score.name = 'score file'
    params = AbsenseParameters(
        distances=dist,
        bitscores=score,
        e_value='e value',
        predict_all=False,
        plot_all=False,
        include_only=include_only,
        gene_lengths=gen_lens,
        db_lengths=db_lens,
        out_dir=None,
        start_time='start time',
    )

    with rec.open() as recorder:
        recorder.write_info(params)

        for key, file in files.items():
            if key != './run_info.txt':
                assert file.getvalue() == ''
                continue

            genelen = "Gene length (for all genes): 400 (default)\n"
            if genelenfile is not None:
                genelen = f"Gene length file: {genelenfile}\n"

            dblen = "Database length (for all species): 8000000 (default)\n"
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

def test_plot(default_params, mocker):
    mocker.patch('abSENSE.recorder.os.makedirs')
    default_params.plot_all = True
    recorder = FileRecorder(default_params, [])
    result = StringIO(
        ',distance,score,p_values,low_interval,high_interval,ambiguous,in_fit,prediction\n'
        'S_cer,0.0,2284.0,0.0,2260.2,2467.68,False,True,2359.91\n'
        'S_par,0.05,2254.0,0.0,2187.99,2378.63,False,True,2279.92\n'
        'S_mik,0.09,2234.0,0.0,2136.21,2316.95,False,True,2223.59\n'
        'S_kud,0.104,2197.0,0.0,2113.86,2291.22,False,True,2199.67\n'
        'S_bay,0.108,2200.0,0.0,2108.34,2284.71,False,True,2193.73\n'
        'S_castellii,0.36,1895.0,0.0,1765.16,1928.92,False,True,1846.31\n'
        'K_waltii,0.49,1767.0,0.0,1600.4,1778.84,False,True,1689.81\n'
        'A_gossyppi,0.52,1736.0,0.0,1571.01,1753.6,False,True,1662.61\n'
        'K_lactis,0.56,1691.0,0.0,1525.83,1711.85,False,True,1619.34\n'
        'A_nidulans,0.9,1099.0,0.0,1164.61,1395.01,False,True,1281.56\n'
        'S_pombe,0.92,1194.0,0.0,1147.89,1379.32,False,True,1265.21\n'
        'Y_lipolytica,0.95,1297.0,0.0,1119.85,1353.19,False,True,1238.13\n'
    )
    df = pd.read_csv(result, index_col=0)

    result = StringIO()
    recorder.plot('GENE', df, 2359.9123, 0.6761,
                  correlation=0.9, bit_threshold=40, outfile=result)
    assert 'DOCTYPE svg' in result.getvalue()
