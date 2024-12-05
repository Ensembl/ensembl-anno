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
check_file(rnasamba_output_path)
rnasamba_output_path = validation_dir / "rnasamba.tsv.txt"
rnasamba_volume = f"{validation_dir}/:/app:rw"
    rnasamba_cmd = [
        "singularity",
        "exec",
        "--bind",
        str(rnasamba_volume),
        str(rnasamba_bin),
        "rnasamba",
        "classify",
        str(rnasamba_output_path),
        str(cdna_file),
        str(rnasamba_weights),
    ]
    logger.info(" ".join(rnasamba_cmd))
    subprocess.run(rnasamba_cmd, check=True)
    
    
    def read_rnasamba_results(file_path):

    results = []

    file_in = open(file_path)
    line = file_in.readline()
    while line:
        line = line.rstrip()
        match = re.search(r"^sequence_name", line)
        if match:
            line = file_in.readline()
            continue

        eles = line.split("\t")
        if not len(eles) == 3:
            line = file_in.readline()
            continue

        transcript_id = eles[0]
        coding_probability = eles[1]
        coding_potential = eles[2]
        results.append([transcript_id, coding_probability, coding_potential])
        line = file_in.readline()
    file_in.close()

    return results