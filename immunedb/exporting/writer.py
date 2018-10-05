import io
import types
import zipfile


class ExportWriter:
    def __init__(self, zipped=False):
        self.zipped = zipped
        self.filename = None

    def get_handle(self, fn, mode='s'):
        class WriterFH:
            def __init__(self, writer):
                self.writer = writer

            def __enter__(self):
                self.buffer = io.BytesIO() if mode == 'b' else io.StringIO()
                return self.buffer

            def __exit__(self, exc_type, exc_value, traceback):
                self.writer.set_filename(fn)
                self.writer.write(self.buffer.getvalue())

        return WriterFH(self)

    def set_filename(self, fn, mode='w+'):
        if self.zipped:
            self.filename = fn
        else:
            self.fh = open(fn, mode)

    def write(self, data):
        if self.zipped:
            self.zip.writestr(
                self.filename,
                ''.join(data) if type(data) == str else data
            )
        else:
            if type(data) == str:
                for chunk in data:
                    self.fh.write(chunk)
            else:
                self.fh.write(data)

    def get_zip_value(self, write_fn=None):
        if self.zipped:
            self.zip.close()
            val = self.out.getvalue()
            if write_fn:
                with open(write_fn, 'wb+') as fh:
                    fh.write(val)
            else:
                return val

    def __enter__(self):
        if self.zipped:
            self.out = io.BytesIO()
            self.zip = zipfile.ZipFile(self.out, 'w')
        else:
            self.filename = None
            self.fh = None
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if self.zipped:
            self.zip.close()
        else:
            if self.fh:
                self.fh.close()
