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
"""Sequence module"""

__all__ = ["Alphabet", "Sequence"]

from dataclasses import dataclass, field
from enum import Enum
from functools import cache

@cache
def _get_sequence() -> str:
    raise NotImplementedError

class Alphabet(Enum):
    """
    Enum class to represent a Sequence's alphabet
    """
    DNA: 1
    RNA: 2
    PROTEIN: 3

@dataclass
class Sequence():
    """
    The Sequence data type describes a contiguous set of
    nucleic acid or amino acid residues in biological molecules,
    such as DNA, RNA, or proteins
    """
    seq_id: str
    seq: str = field(repr=False)
    alphabet: Alphabet = field(default=Alphabet.DNA)
