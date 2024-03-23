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
"""Biotype module"""

__all__ = [ 'Biotype' ]

from dataclasses import dataclass, field

@dataclass
class Biotype():

    __type = 'biotype'

    internal_id: int = field(repr=False)
    name: str
    object_type: str
    biotype_group: str = field(repr=False)
    so_acc: str
    so_term: str
    description: str = field(repr=False)
    db_type: str = field(repr=False)
    attrib_type_id: int = field(repr=False, default=None)

    @property
    def type(self) -> str:
        self.__type



