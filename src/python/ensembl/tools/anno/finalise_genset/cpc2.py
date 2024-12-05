# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
check_file(cpc2_output_path)
cpc2_output_path = validation_dir / "cpc2.tsv"
    
    cpc2_volume = f"{validation_dir}/:/app:rw"
    cpc2_cmd = [
        "singularity",
        "exec",
        "--bind",
        str(cpc2_volume),
        str(cpc2_bin),
        "python3",
        "/CPC2_standalone-1.0.1/bin/CPC2.py",
        "-i",
        str(cdna_file),
        "--ORF",
        "-o",
        str(cpc2_output_path),
    ]
    logger.info(" ".join(cpc2_cmd))
    subprocess.run(cpc2_cmd, check=True)
    cpc2_output_path = f"{str(cpc2_output_path)}.txt"
    
    
def read_cpc2_results(file_path):

    results = []

    file_in = open(file_path)
    line = file_in.readline()
    while line:
        line = line.rstrip()
        match = re.search(r"^#ID", line)
        if match:
            line = file_in.readline()
            continue

        eles = line.split("\t")
        if not len(eles) == 9:
            line = file_in.readline()
            continue

        transcript_id = eles[0]
        transcript_length = eles[1]
        peptide_length = eles[2]
        coding_probability = eles[7]
        coding_potential = eles[8]
        results.append(
            [
                transcript_id,
                coding_probability,
                coding_potential,
                transcript_length,
                peptide_length,
            ]
        )
        line = file_in.readline()
    file_in.close()

    return results    