from abSENSE.recorder import FileRecorder


class FitResult():
    def __init__(self, gene):
        self.gene = gene

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def record_to(self, recorder: FileRecorder):
        recorder.write_gene(self.gene)

    def status(self):
        return ''

class ErrorResult(FitResult):
    def __init__(self, gene, predictions=None):
        super().__init__(gene)
        self.predictions = predictions

    def record_to(self, recorder: FileRecorder):
        super().record_to(recorder)
        recorder.analysis_error(predictions=self.predictions)

    def status(self):
        return 'Analysis Error'


class NotEnoughDataResult(FitResult):
    def __init__(self, gene):
        super().__init__(gene)

    def record_to(self, recorder: FileRecorder):
        super().record_to(recorder)
        recorder.not_enough_data()

    def status(self):
        return 'Not Enough Data'


class SampledResult(FitResult):
    def __init__(self, gene, result, a_fit, b_fit, bit_threshold, correlation):
        super().__init__(gene)
        self.result = result
        self.a_fit = a_fit
        self.b_fit = b_fit
        self.bit_threshold = bit_threshold
        self.correlation = correlation

    def record_to(self, recorder: FileRecorder):
        super().record_to(recorder)
        recorder.plot(
            gene=self.gene,
            result=self.result,
            a_fit=self.a_fit,
            b_fit=self.b_fit,
            correlation=self.correlation,
            bit_threshold=self.bit_threshold,
        )


        recorder.write_params(self.a_fit, self.b_fit)

        for _, row in self.result.round(2).iterrows():
            recorder.write_result(
                prediction=row.prediction,
                high=row.high_interval,
                low=row.low_interval,
                pval=row.p_values,
                realscore=row.score,
                is_considered=row.in_fit,
                is_ambiguous=row.ambiguous,
            )

        recorder.finalize_row()

    def status(self):
        return ''
