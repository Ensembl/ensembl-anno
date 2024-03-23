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
"""Exon module"""

__all__ = ["Exon"]

from typing import Union, Self
import warnings
from . import Feature, Location, Slice, Strand

class Exon(Feature):
    """Representation of an exon.
    The location of an exon is always provided 5' -> 3' on the forward strand.
    Thus the start should always be lower than the end.
    We use phase and end phase to provide information about the coding potential
    and where the next codon start at the start of the exon.
    Attributes:
      start: Start position of the exon.
      end: End position of the exon.
      strand: Strand of the exon, either + or -.
      location_name: Name of the region the exon is on.
      fasta_file: Path to a FASTA file containing the DNA of the region.
      sequence: DNA sequence of the exon.
      public_identifier: Name of the exon.
      exon_start_phase: Phase of the exon, either -1, 0, 1 or 2.
      exon_end_phase: End phase of the exon, either -1, 0, 1 or 2.
    Raises:
      Exception: when start is greater than end
    """

    __type = 'exon'

    def __init__(self,
                 location: Location,
                 strand: Strand,
                 phase: int,
                 end_phase: int,
                 gen_slice: Slice = None,
                 internal_id: Union[str, int] = None,
                 stable_id: str = None,
                 version: int = None,
                 analysis: str = None,
                 is_constitutive: bool = False,
                 is_current: bool = True
                ) -> None:
        if not isinstance(phase, int):
            raise ValueError("phase argument must be int")
        if phase not in (-1, 0, 1, 2):
            raise ValueError(f"Bad value {phase} for exon phase: it must be any of (-1, 0, 1, 2)")
        if end_phase is None:
            raise ValueError("No end phase set in Exon. You must set it explicitly.")
        self._stable_id = stable_id
        self._version = version
        self._slice = gen_slice
        self._phase = phase
        self._end_phase = end_phase
        self._is_constitutive = is_constitutive
        self._is_current = is_current
        self._internal_id = internal_id
        super().__init__(location, strand, gen_slice, analysis, internal_id)

    @classmethod
    def fastinit(cls, start: int, end: int, strand: int, gen_slice: Slice,
                 phase: int = -1, end_phase: int = -1) -> Self:
        loc = Location(start, end)
        return cls(loc, Strand(strand), phase, end_phase, gen_slice)

    def __repr__(self) -> str:
        eid = self.stable_id if self._stable_id else self._slice.region.name
        return f'{self.__class__.__name__}({eid}:{self.start}:{self.end}:{self.strand.value}\
            -{self.phase}:{self._end_phase})'

    @property
    def location(self) -> Location:
        return self._location

    @property
    def start(self) -> int:
        return self._location.start

    @property
    def end(self) -> int:
        return self._location.end

    @property
    def stable_id_version(self) -> str:
        if self._stable_id:
            return f"{self._stable_id}.{self._version}"
        return None

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
    def version(self):
        return self._version

    @version.setter
    def version(self, value) -> None:
        self._version = value

    @property
    def phase(self) -> int:
        return self._phase

    @phase.setter
    def phase(self, value: int) -> None:
        if value not in (-1, 0, 1, 2):
            raise ValueError(f"Bad value {value} for exon phase: it must be any of (-1, 0, 1, 2)")
        self._phase = value

    @property
    def end_phase(self) -> int:
        if self._end_phase:
            return self._end_phase
        exon_id = self._stable_id if self._stable_id else self._internal_id
        msg = f"No end phase set in Exon {exon_id}. You must set it explicitly."
        warnings.warn(msg, UserWarning)
        return None

    @property
    def frame(self) -> int:
        if self._phase == -1:
            return '.' # gff convention for no frame info
        if self._phase == 0:
            return self._start%3
        if self._phase == 1:
            return (self._start + 2)%3
        if self._phase == 2:
            return (self._start + 1)%3
        raise ValueError(f"bad phase in exon {self._phase}")

    @property
    def is_constitutive(self) -> bool:
        return self._is_constitutive

    @is_constitutive.setter
    def is_constitutive(self, value: bool) -> None:
        if not isinstance(value, bool):
            raise ValueError("New value must be boolean")
        self._is_constitutive = value

    @property
    def exon_string(self) -> str:
        return f"{self.start}:{self.end}:{self.strand.value}"

    @property
    def seq(self) -> str:
        return self._sequence.seq

    @seq.setter
    def seq(self, value) -> None:
        raise NotImplementedError()

    def adjust_start_end(self, start_adj: int, end_adj: int) -> Self:
        """
        What's the meaning of a new exon with same IDs and different
        coordinates? Wouldn't be better to reset the IDs or change the
        initial exon?
        To be confirmed!!
        """
        cur_loc = self._location
        if self._strand == Strand.FORWARD:
            new_loc = Location(cur_loc.start + start_adj,
                               cur_loc.end + end_adj)
        else:
            new_loc = Location(cur_loc.start - end_adj,
                               cur_loc.end - start_adj)
        return Exon(
            new_loc,
            self._strand,
            self._phase,
            self._end_phase,
            self._slice,
            self._internal_id,
            self._stable_id,
            self._version,
            self._analysis,
            self._is_constitutive,
            self._is_current
        )
