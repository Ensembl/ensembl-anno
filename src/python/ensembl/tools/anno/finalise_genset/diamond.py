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

 logger.info("diamond validation")
    diamond_results = None
    if diamond_validation_db is not None:
        diamond_output_dir = create_dir(validation_dir, "diamond_output")
        diamond_validation(
            diamond_validation_db,
            amino_acid_file,
            diamond_output_dir,
            num_threads,
        )
        diamond_results = read_diamond_results(diamond_output_dir)
        
        
        def diamond_validation(
    diamond_validation_db: Path,
    amino_acid_file: Path, 
    diamond_output_dir: Path, 
    num_threads: int
)-> None:

    batched_protein_files = split_protein_file(amino_acid_file, diamond_output_dir, 100)

    pool = multiprocessing.Pool(int(num_threads)) # pylint: disable=consider-using-with
    for batched_protein_file in batched_protein_files:
        pool.apply_async(
            multiprocess_diamond,
            args=(
                batched_protein_file,
                diamond_output_dir,
                diamond_validation_db,
            ),
        )
    pool.close()
    pool.join()


def multiprocess_diamond(
    batched_protein_file,
    diamond_output_dir,
    diamond_validation_db,
)->None:

    #batch_num = os.path.splitext(batched_protein_file)[0]
    #batch_dir = os.path.dirname(batched_protein_file)
    diamond_output_file = f"{batched_protein_file}.dmdout"
    logger.info("Running diamond on %s :", batched_protein_file)

    diamond_cmd = [
        "diamond",
        "blastp",
        "--query",
        batched_protein_file,
        "--db",
        diamond_validation_db,
        "--out",
        diamond_output_file,
    ]

    logger.info(" ".join(diamond_cmd))
    subprocess.run(diamond_cmd)
    subprocess.run(["mv", diamond_output_file, diamond_output_dir])
    
    
def read_diamond_results(diamond_output_dir):

    results = []
    diamond_files = glob.glob(diamond_output_dir + "/*.dmdout")
    for file_path in diamond_files:
        file_in = open(file_path)
        line = file_in.readline()
        while line:
            line = line.rstrip()

            eles = line.split("\t")
            if not len(eles) == 12:
                line = file_in.readline()
                continue

            transcript_id = eles[0]
            e_value = eles[10]
            results.append([transcript_id, e_value])
            line = file_in.readline()
    file_in.close()

    return results
    