# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Transcript module"""

__all__ = [ 'Transcript' ]

from typing import Union
from . import Analysis, Biotype, CoordinateSystem, Exon, Location, Slice, Feature, Strand, Translation

class Transcript(Feature):
    """Representation of a transcript.

    TO DO!!!

    The location of a transcript is always provided 5' -> 3' on the forward strand.
    Thus the start should always be lower than the end.

    Attributes:
      exons: List of exons constituting the transcript.
      start: Start position of the transcript.
      end: End position of the transcript.
      strand: Strand of the transcript, either + or -.
      location_name: Name of the region the transcript is on.
      fasta_file: Path to a FASTA file containing the DNA of the region.
      sequence: DNA sequence of the transcript.
      external_id: Name of the transcript.
      introns: List of introns constituting the transcript, should be 0 if exons size is one.
      cds_genomic_start: Genomic start position of the CDS, always be lower than cds_genomic_end.
      cds_genomic_end: Genomic end position of the CDS, always be greater than cds_genomic_start.
      cds_sequence: Translateable cDNA sequence.
      translation_sequence: Protein sequence translated from cds_sequence.
    """

    __type = 'transcript'

    def __init__(self,
                 stable_id: str,
                 version: int = 1,
                 internal_id: Union[str, int] = None,
                 gen_slice: Slice = None,
                 location: Location = None,
                 strand: Strand = None,
                 exons: list[Exon] = None,
                 analysis: Analysis = None,
                 biotype: Biotype = None,
                 source: str = 'ensembl',
                 is_current: bool = True,
                 is_canonical: bool = False
                ) -> None:

        self._stable_id = stable_id
        self._version = version
        self._exons = [] if not exons else exons
        self._internal_id = internal_id
        self._biotype = biotype
        self._source = source
        self._is_current = is_current
        self._is_canonical = is_canonical
        self._translation = None
        self._seq: str = None
        self._coding_region_start: int = None
        self._coding_region_end: int = None

        super().__init__(location, strand, gen_slice, analysis, internal_id)

    def __repr__(self) -> str:
        if self._stable_id:
            if self._attributes.get('is_canonical'):
                return f'{self.__class__.__name__}({self.stable_id} - canonical)'
            return f'{self.__class__.__name__}({self.stable_id})'
        return f'{self.__class__.__name__}(internal_id{self._internal_id})'

    @property
    def stable_id_version(self) -> str:
        return f"{self._stable_id}.{self._version}"

    @property
    def stable_id(self) -> str:
        return self._stable_id

    @stable_id.setter
    def stable_id(self, value: str) -> None:
        if value.find('.') > 0:
            (self._stable_id, self._version) = value.split('.')
        else:
            self._stable_id = value

    @property
    def type(self) -> str:
        return self.__type

    @property
    def version(self):
        return self._version

    @version.setter
    def version(self, value) -> None:
        self._version = value

    @property
    def internal_id(self) -> int:
        return self._internal_id

    @property
    def biotype(self) -> Biotype:
        return self._biotype

    @biotype.setter
    def biotype(self, value: Biotype) -> None:
        self._biotype = value

    @property
    def source(self) -> str:
        return self._source

    @source.setter
    def source(self, value: str) -> None:
        self._source = value

    @property
    def is_current(self) -> bool:
        return self._is_current

    @is_current.setter
    def is_current(self, value: bool) -> None:
        self._is_current = value

    @property
    def translation(self) -> Translation:
        return self._translation

    @translation.setter
    def translation(self, value: Translation) -> None:
        if not isinstance(value, Translation):
            raise AttributeError(f'Provided value must be a Translation, instead of {type(value)}')
        self._translation = value

    def set_translation(self, tl: Translation) -> None:
        self.translation = tl

    def get_exons(self, constitutive: bool = False) -> list[Exon]:
        if constitutive:
            return [ e for e in self._exons if e.is_constitutive ]
        return self._exons

    def set_exons(self, exons: list[Exon]) -> None:
        self._exons = exons

    def set_attribs(self, attribs: dict[str, str]) -> None:
        self._attributes = attribs

    def get_attribs(self, att_code: str = None) -> Union[dict, tuple]:
        if att_code is None:
            return self._attributes
        return tuple(att_code, self._attributes.get(att_code))

    def add_attrib(self, code: str, value: str) -> None:
        self._attributes[code] = value

    def is_canonical(self) -> bool:
        return self._is_canonical

    def is_mane_select(self) -> bool:
        return 'MANE_Select' in self._attributes

    def is_mane(self) -> bool:
        return 'mane' in self._attributes.keys().lower()

    @property
    def seq(self) -> str:
        raise NotImplementedError()

    @seq.setter
    def seq(self, value) -> None:
        raise NotImplementedError()

    def coding_region_start(self) -> int:
        if self._coding_region_start:
            return self._coding_region_start
        if not self._translation:
            return None
        if self._translation.coord_system != CoordinateSystem.GENOMIC:
            raise NotImplementedError

        start = 0
        strand = self._translation.start_exon.strand
        if strand == Strand.FORWARD:
            start = self._translation.start_exon.start
            start += self._translation.start - 1
        else:
            start = self._translation.end_exon.end
            start -= self._translation.end - 1
        self._coding_region_start = start
        return start

    def coding_region_end(self) -> int:
        if self._coding_region_end:
            return self._coding_region_end
        if not self._translation:
            return None
        if self._translation.coord_system != CoordinateSystem.GENOMIC:
            raise NotImplementedError

        end = 0
        strand = self._translation.start_exon.strand
        if strand == Strand.FORWARD:
            end = self._translation.end_exon.start
            end += self._translation.end - 1
        else:
            end = self._translation.start_exon.end
            end -= self._translation.start - 1
        self._coding_region_end = end
        return end

    def get_all_translateable_exons(self) -> list[Exon]:
        if not self._translation:
            return []
        trl = self._translation
        ex_list = []
        for e in self._exons:
            adjust = 0
            if e == trl.start_exon:
                adjust = trl.start - 1
                if adjust != 0:
                    ex_list.append(e.adjust_start_end(adjust, 0))
                    continue
            if e == trl.end_exon:
                adjust = trl.end - e.length
                if adjust != 0:
                    ex_list.append(e.adjust_start_end(0, adjust))
                    break
            ex_list.append(e)
        return ex_list
