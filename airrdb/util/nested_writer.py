import cStringIO as StringIO
import csv


class NestedCSVWriter(object):
    def __init__(self, headers, mapping=None, streaming=False):
        self._headers = headers
        self._mapping = mapping
        self._streaming = streaming
        self._out = StringIO.StringIO()
        self._csv = csv.DictWriter(self._out, fieldnames=headers)
        self._csv.writeheader()

    def add_row(self, row, write_default=False, default=None,
                write_if_stream=True):
        if write_default:
            write_dict = {k: default for k in self._headers}
        else:
            write_dict = {}
        for header in self._headers:
            if self._mapping is not None and header in self._mapping:
                write_dict[header] = self._mapping[header](row)
            elif header in row:
                write_dict[header] = row[header]
        return self.add_raw_row(write_dict, write_if_stream=write_if_stream)

    def add_raw_row(self, raw_dict, write_if_stream=True):
        self._csv.writerow(raw_dict)
        if self._streaming and write_if_stream:
            return self.get_value()

    def get_value(self):
        self._out.seek(0)
        v = self._out.read()
        if self._streaming:
            self._out.seek(0)
            self._out.truncate()
        return v
