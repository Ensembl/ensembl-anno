.. See the NOTICE file distributed with this work for additional information
   regarding copyright ownership.

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

API Setup and installation
===========================

Requirements
--------------

.. _install:

An Ensembl API checkout including:

- ensembl-production `ensembl-production <https://github.com/Ensembl/ensembl-production>`_.
- ensembl-analysis `ensembl-analysis <https://github.com/Ensembl/ensembl-analysis/tree/dev/hive_master>`_. (on dev/hive_master branch)
- ensembl-taxonomy `ensembl-taxonomy <https://github.com/Ensembl/ensembl-taxonomy>`_.
- ensembl-orm `ensembl-orm <https://github.com/Ensembl/ensembl-orm>`_.

Software
^^^^^^^^

#. Python 3.8+
#. Bioperl 1.6.9+

Python Modules
^^^^^^^^^^^^^^
#. argschema



Installation
------------
Directly from GitHub:

.. code-block:: none
   :linenos:

   git clone https://github.com/Ensembl/ensembl-analysis -b experimental/gbiab
   git clone https://github.com/Ensembl/ensembl-production
   git clone https://github.com/Ensembl/ensembl-hive
   git clone https://github.com/Ensembl/ensembl-taxonomy
   git clone https://github.com/Ensembl/ensembl-orm