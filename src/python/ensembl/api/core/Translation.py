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
"""Translation module"""

__all__ = [ "Translation" ]

from . import CoordinateSystem, Exon, Location

class Translation():
    def __init__(self,
                 location: Location,
                 start_exon: Exon,
                 end_exon: Exon,
                 internal_id: int = None) -> None:
        self._location = location
        self._start_exon = start_exon
        self._end_exon = end_exon
        self._internal_id = internal_id

    def __repr__(self) -> str:
        return ""

    @property
    def internal_id(self):
        return self._internal_id

    @property
    def start(self) -> int:
        return self._location.start

    @property
    def end(self) -> int:
        return self._location.end

    @property
    def coord_system(self) -> CoordinateSystem:
        return self._location.coordinate_system

    @property
    def start_exon(self) -> Exon:
        return self._start_exon

    @property
    def end_exon(self) -> Exon:
        return self._end_exon
