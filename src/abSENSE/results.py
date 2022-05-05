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
    def __init__(self, gene, result, a_fit, b_fit):
        super().__init__(gene)
        self.result = result
        self.a_fit = a_fit
        self.b_fit = b_fit

    def record_to(self, recorder: FileRecorder):
        super().record_to(recorder)
        recorder.write_params(self.a_fit, self.b_fit)

        for _, row in self.result.iterrows():
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

        recorder.plot(self.gene, self.result, self.a_fit, self.b_fit)

    def status(self):
        return ''
