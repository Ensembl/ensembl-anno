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
"""Generic Feature module"""

__all__ = [ "Feature", "Strand" ]

from enum import Enum
from typing import Union
import warnings
from . import Analysis, Location, RegionType, Sequence, Slice

class Strand(Enum):
    """
    Enum class to represent the strand-ness
    """
    FORWARD = 1
    UNDEFINED = 0
    REVERSE = -1

class Feature():
    """
    Ensembl specific sequence feature.
    This is the Base feature class from which all Ensembl features inherit.
    It provides a bare minimum functionality that all features require.  It
    basically describes a location on a sequence in an arbitrary coordinate
    system.
    """
    def __init__(self,
            location: Location,
            strand: Strand,
            reg_slice: Slice = None,
            analysis: Analysis = None,
            internal_id: int = None,
            sequence: Sequence = None
        ) -> None:

        if reg_slice and not isinstance(reg_slice, Slice):
            raise ValueError(f"Slice argument must be of type {type(Slice)}")
        # if analysis and not isinstance(analysis, Slice):
        #     raise ValueError('Analysis argument must be a Bio::EnsEMBL::Analysis')
        if not location or not isinstance(location, Location):
            raise ValueError(f"Location argument must be defined and of type {type(Location)}")
        if not strand or not isinstance(strand, Strand):
            raise ValueError(f"Location argument must be defined and of type {type(Strand)}")
        if location.start > location.end:
            raise ValueError(f"Invalid Location: start({location.start}) > end({location.end})")
        try:
            if not reg_slice.is_circular and location not in reg_slice.location:
                raise ValueError(f"Feature {internal_id} location not contained in given Slice")
        except TypeError:
            pass

        self._location = location
        self._slice = reg_slice
        self._strand = strand
        self._analysis = analysis
        self._internal_id = internal_id
        self._sequence = sequence
        self._attributes = {}

    @classmethod
    def fastinit(cls, start: int, end: int, length: int, analysis_name: str,
                 strand: int, reg_slice: Slice = None, internal_id: int = None,
                 sequence: str = None):
        if not start:
            raise ValueError("Feature start must be specified")
        if not end and not length:
            raise ValueError("Either feature end or length must be specified")
        if not end:
            end = start + length - 1
        an = Analysis(analysis_name) if analysis_name else None
        loc = Location(start, end)
        st = Strand(strand)
        seq = Sequence(seq_id=None, seq=sequence) if sequence else None
        return cls(loc, st, reg_slice, an, internal_id, seq)

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}(internal_id{self._internal_id}:\
            {self._slice.region.name}:{self.start}:{self.end})'

    @property
    def start(self) -> int:
        return self._location.start

    @start.setter
    def start(self, value: int) -> None:
        if not isinstance(value, int):
            raise ValueError('start argument must be int')
        self._location.start = value

    @property
    def end(self) -> int:
        return self._location.end

    @end.setter
    def end(self, value: int) -> None:
        if not isinstance(value, int):
            raise ValueError('end argument must be int')
        self._location.end = value

    @property
    def strand(self) -> Strand:
        return self._strand

    @strand.setter
    def strand(self, value: Union[Strand, int]) -> None:
        if not value:
            warnings.warn("Provided value is empty or None", UserWarning)
            return
        if not isinstance(value, Strand) or not isinstance(value, int):
            raise ValueError('strand argument must be Strand or int')
        self._strand = value if isinstance(value, Strand) else Strand(value)

    @property
    def analysis(self) -> Analysis:
        return self._analysis

    @analysis.setter
    def analysis(self, value) -> None:
        self._analysis = value

    @property
    def internal_id(self) -> int:
        return self._internal_id

    @internal_id.setter
    def internal_id(self, value: int) -> None:
        self._internal_id = value

    @property
    def length(self) -> int:
        raw_len = self._location.end - self._location.start
        if raw_len >= 0:
            return raw_len + 1
        if self._slice.is_circular():
            # if circular, we can work out the length of an origin-spanning
            # feature using the size of the underlying region.
            return  raw_len%self._slice.length + 1
        raise ValueError('Cannot determine length of non-circular feature where start > end')

    def set_attribs(self, attribs: dict[str, str]) -> None:
        self._attributes = attribs

    def get_attrib(self, att_code: str = None) -> dict:
        if att_code is None:
            return self._attributes
        return {att_code: self._attributes.get(att_code)}

    def add_attrib(self, code: str, value: str) -> None:
        self._attributes[code] = value

    def get_slice(self) -> Slice:
        return self._slice

    def set_slice(self, gen_slice: Slice) -> None:
        self._slice = gen_slice

    def feature_slice(self) -> Slice:
        if not self._slice:
            warnings.warn("Cannot create feature Slice for feature without an attached slice")
            return None
        return Slice(self._slice.region, self._location)

    @property
    def seq_region_name(self) -> str:
        return self._slice.region.name if self._slice else ''


    #  FROM SLICE - TO BE CONFIRMED!!
    @property
    def seq_region_length(self) -> int:
        if self._slice:
            return self._slice.length
        warnings.warn(f"There is no slice for this feature {self.internal_id}")
        return None

    @property
    def seq_region_type(self) -> RegionType:
        if self._slice:
            return self._slice.region_type
        return None

    @property
    def coord_system_name(self) -> str:
        if self._slice:
            return self._slice.location.coordinate_system.name
        return None

    # def get_summary(self) -> dict:
    #     summary = {}
    #     summary['id'] = self.internal_id
    #     summary['version'] = self._version if self._version else ''
    #     summary['start'] = self.seq_region_start()
    #     summary['end'] = self.seq_region_end()
    #     summary['strand'] = self._strand
    #     summary['seq_region_name'] = self.seq_region_name()
    #     return summary
