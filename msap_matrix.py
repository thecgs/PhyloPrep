"""Disk-backed alignment access for large concatenated FASTA matrices."""

import tempfile


CHUNK_SIZE = 65535  # A whole number of codons.


class AlignmentSpool:
    """Validate and spool FASTA in bounded chunks, including unwrapped records.

    Callers place the spool beside their outputs, so browser working files
    stay on the mounted disk filesystem rather than Pyodide's /tmp MEMFS.
    Only IDs, text seek cookies and sequence lengths are retained in memory.
    """

    def __init__(self, infile, codon=False, transform=None, directory=None):
        self.records = []
        self.length = None
        self.file = tempfile.TemporaryFile(mode='w+t', encoding='utf-8',
                                           dir=directory)
        identifier, start, length = None, 0, 0
        ids = set()

        def finish():
            if identifier is None:
                return
            if not length:
                raise ValueError(f'{identifier}: sequence is empty')
            if self.length is None:
                self.length = length
            elif length != self.length:
                raise ValueError(f'Alignment length mismatch: {identifier} has length {length}, expected {self.length}')
            self.records.append((identifier, start, length))

        try:
            with open(infile, encoding='utf-8') as source:
                line_start = True
                while True:
                    part = source.readline(CHUNK_SIZE)
                    if not part:
                        break
                    if line_start and part.startswith('>'):
                        finish()
                        # Headers are normally short; sequence lines are always bounded.
                        if not part.endswith('\n'):
                            part += source.readline()
                        header = part[1:].strip()
                        if not header:
                            raise ValueError('FASTA record has an empty identifier')
                        identifier = header.split()[0]
                        if identifier in ids:
                            raise ValueError(f'Duplicate sequence ID: {identifier}')
                        ids.add(identifier)
                        start, length = self.file.tell(), 0
                    else:
                        chunk = ''.join(part.split())
                        if chunk:
                            if identifier is None:
                                raise ValueError('FASTA sequence data appears before its first header')
                            original_length = len(chunk)
                            if transform:
                                chunk = transform(identifier, chunk, length)
                            self.file.write(chunk)
                            length += original_length
                    line_start = part.endswith('\n')
            finish()
            if not self.records:
                raise ValueError('No sequences found in input FASTA.')
            if codon and self.length % 3:
                raise ValueError(f'Codon alignment length ({self.length}) is not divisible by 3.')
            self.file.flush()
        except BaseException:
            self.close()
            raise

    def chunks(self, index, size=CHUNK_SIZE):
        _, position, remaining = self.records[index]
        while remaining:
            self.file.seek(position)
            chunk = self.file.read(min(size, remaining))
            if not chunk:
                raise OSError('Incomplete matrix spool')
            position = self.file.tell()
            remaining -= len(chunk)
            yield chunk

    def close(self):
        self.file.close()

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()
