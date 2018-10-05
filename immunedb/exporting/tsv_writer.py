import csv
import io


class StreamingTSV(object):
    def __init__(self, fieldnames):
        self.out = io.StringIO()
        self.fieldnames = fieldnames
        self.writer = csv.DictWriter(self.out, fieldnames=self.fieldnames,
                                     delimiter='\t', extrasaction='ignore')

    def writeheader(self):
        return self._stream(self.writer.writeheader)

    def writerow(self, row):
        return self._stream(self.writer.writerow, row)

    def _stream(self, action, *args, **kwargs):
        self.out.truncate(0)
        self.out.seek(0)
        action(*args, **kwargs)
        return self.out.getvalue().strip() + '\n'
