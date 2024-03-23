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
"""Slice module"""

__all__ = ["Slice", "Location", "LocationModifier", "Region", 
           "RegionType", "Topology", "CoordinateSystem"]

from dataclasses import dataclass, field
from enum import Enum
import warnings
from . import Sequence, Assembly, Alphabet, Strand

class Topology(Enum):
    """
    Enum class to represent a region's topology
    """
    LINEAR = 1
    CIRCULAR = 2

class CoordinateSystem(Enum):
    """
    Enum class to represent the coordinate system
    """
    GENOMIC = 1

class RegionType(Enum):
    """
    Enum class to represent the region type
    """
    CHROMOSOME = 1
    SCAFFOLD = 2
    PLASMID = 3
    CONTIG = 4

@dataclass
class LocationModifier():
    """
    Captures if a Feature is incomplete on either its
    5' or 3' end
    """
    partial_start: bool
    partial_end: bool

class Location():
    """
    Dataclass to store coordinates to position a Feature
    on a Region
    """
    def __init__(self, start: int, end: int,
                 coord_sys_type: CoordinateSystem = CoordinateSystem.GENOMIC,
                 location_modifier: LocationModifier = None) -> None:
        if not start or not end:
            raise ValueError("Both Location start and end must be specified")
        if self.end + 1 < self.start:
            raise ValueError("Location start must be less than end")
        if coord_sys_type != CoordinateSystem.GENOMIC:
            raise NotImplementedError("Can only manage genomic coordinates!")
        self._start = start
        self._end = end
        self._coord_system = coord_sys_type
        self._location_modifier = location_modifier

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({self.start}-{self.end}-{self.coordinate_system.name})'

    def __str__(self) -> str:
        return f'{self.__class__.__name__}({self.start}-{self.end})'

    def __contains__(self, other) -> bool:
        if isinstance(other, Location):
            if other._coord_system != self._coord_system:
                warnings.warn("Cannot compare Location objects on different coordinate systems")
                return False
            return self.start <= other.start and other.end <= self.end
        if isinstance(other, int):
            return self.start <= other <= self.end
        warnings.warn(f"Comparing Location with another object ({type(other)})")
        return False

    def __eq__(self, other) -> bool:
        if not isinstance(other, Location):
            warnings.warn(f"Cannot compare Location to object type {type(other)}")
            return False
        if other._coord_system != self._coord_system:
            warnings.warn("Cannot compare Location objects on different coordinate systems")
            return False
        if self.start == other.start and self.end == other.end:
            return True
        return False

    def __lt__(self, other) -> bool:
        if not isinstance(other, Location):
            warnings.warn(f"Cannot compare Location to object type {type(other)}")
            return False
        return self.end < other.start

    @property
    def coordinate_system(self) -> CoordinateSystem:
        return self._coord_system

    @property
    def start(self) -> int:
        return self._start

    @start.setter
    def start(self, value) -> None:
        try:
            if int(value) > 0:
                self._start = int(value)
        except TypeError:
            warnings.warn("Start value must be positive integer")
        except ValueError:
            warnings.warn("Start value must be positive integer")

    @property
    def end(self) -> int:
        return self._end

    @end.setter
    def end(self, value) -> None:
        try:
            if int(value) > self._start + 1:
                self._end = int(value)
        except TypeError:
            warnings.warn(f"End value ({value}) must be greater than start ({self._start})")
        except ValueError:
            warnings.warn(f"End value ({value}) must be greater than start ({self._start})")

    @property
    def length(self) -> int:
        return self.end - self.start + 1


@dataclass
class Region():
    """
    It describes the coordinate system that contains Features
    """
    name: str
    region_type: RegionType
    length: int
    genetic_code: int = field(repr=False, default=1)
    is_top_level: bool = field(default=True, compare=False)
    rank: int = field(init=False, repr=False, default=None)
    topology: Topology = field(repr=False, default=Topology.LINEAR)
    accession_id: str = field(repr=False, default='', compare=False)
    sequence: Sequence = field(init=False, repr=False, compare=False, default=None)
    assembly: Assembly = field(init=False, repr=False, default=None)
    metadata: dict[str,str] = field(repr=False, compare=False, default={})

class Slice():
    """
    Slice represent a part of a region, inclusive of an orientation (strand),
    in a coordinate system (only GENOMIC for now)
    """
    def __init__(self, region: Region, location: Location, strand: Strand) -> None:
        self._region = region
        self._location = location
        self._strand = strand

    def __repr__(self) -> str:
        pass

    @classmethod
    def fastinit(cls,
                 seq_region_type: str,
                 seq_region_name: str,
                 seq_region_length: int,
                 seq_region_start: int,
                 seq_region_end: int,
                 seq_region_strand: int,
                 seq: str = None,
                 alphabet: str = "DNA"):
        if not seq_region_type:
            raise ValueError("Region type is required")
        if not seq_region_name:
            raise ValueError("Region name argument is required")
        if seq_region_start is None or seq_region_start <= 0:
            raise ValueError("Invalid Region start: must be positive")
        if seq_region_end is None or seq_region_end <= 0:
            raise ValueError("Invalid Region end: must be positive")
        if seq_region_length <= 0:
            raise ValueError("Invalid Region length: must be positive or None")
        seq_len = seq_region_end - seq_region_start + 1
        if not seq_region_length:
            seq_region_length = seq_len
        elif seq_region_length != seq_len:
            raise ValueError(f"Inconsistent provided seq_region_length ({seq_region_length}) \
                             does not match with start and end ({seq_len})")
        if seq and not alphabet:
            raise ValueError("Alphabet is required, when sequence is provided")
        try:
            if len(seq) != seq_len:
                raise ValueError(f"len(seq)={len(seq)} differs from seq_region_length={seq_region_length}")
        except TypeError:
            pass

        reg = Region(seq_region_name, RegionType[seq_region_type], seq_region_length,
                     genetic_code=None)
        loc = Location(seq_region_start, seq_region_end)
        if seq:
            reg.sequence = Sequence(None, Alphabet(alphabet.upper()), seq)
        return cls(reg, loc, Strand(seq_region_strand))

    @property
    def region(self) -> Region:
        """
        Get region
        """
        return self._region

    @property
    def location(self) -> Location:
        """
        Get Location
        """
        return self._location

    @property
    def start(self) -> int:
        return self._location.start

    @property
    def end(self) -> int:
        return self._location.end
    
    @property
    def strand(self) -> Strand:
        return self._strand

    @property
    def length(self) -> int:
        return self._location.length

    @property
    def region_type(self) -> RegionType:
        return self._region.region_type

    @property
    def karyotype_rank(self) -> int:
        return self._region.rank

    @property
    def is_toplevel(self) -> bool:
        return self._region.is_top_level

    @property
    def is_circular(self) -> bool:
        if self._region.topology == Topology.CIRCULAR:
            return True
        return False
    