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


def create_dir(main_output_dir, dir_name):

    if dir_name:
        target_dir = os.path.join(main_output_dir, dir_name)
    else:
        target_dir = main_output_dir

    if os.path.exists(target_dir):
        logger.warning("Directory already exists, will not create again")
        return target_dir

    logger.info("Attempting to create target dir: %s" % target_dir)

    try:
        os.mkdir(target_dir)

    except OSError:
        logger.error("Creation of the dir failed, path used: %s" % target_dir)
    else:
        logger.info(
            "Successfully created the dir on the following path: %s" % target_dir
        )

    return target_dir
