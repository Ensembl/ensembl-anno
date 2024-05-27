import argparse
import pymysql
import csv
from typing import List, Dict, Any, Tuple
from pathlib import Path

# beta meta db details
meta_db_name = "ensembl_genome_metadata"
meta_db_host = "mysql-ens-stage-1"
meta_db_port = 4732
meta_db_user = "ensro"
meta_db_pass = ""

# Function to execute the query and format the results
 execute_query(accession, meta_db_name, meta_db_host, meta_db_port,
                                                        meta_db_user, meta_db_pass
def execute_query(accession:str,database:str, host:str, port:int, user:str,password:str )-> Tuple[str, List[Dict[str, Any]]]:
    """Executes query on db and if the accession is present get some statistics values.

    Args:
        accession (str): gca_accession
        database (str): db name
        host (str): db host
        port (int): db port
        user (str): db user
        password (str): db pass

    Returns:
        Tuple[str, List[Dict[str, Any]]]: accession, query_results
    """ 

    try:
        # Connect to the MySQL database
        connection = pymysql.connect(host=host,
                                    user=user,
                                    port=port,
                                    password=password,
                                    database=database,
                                    cursorclass=pymysql.cursors.DictCursor)
        
        # Create a cursor object
        cursor = connection.cursor()

        # Define the SQL query
        sql_query = """
            SELECT 
                assembly.accession, 
                attribute.name, 
                dataset_attribute.value 
            FROM 
                assembly 
            RIGHT JOIN 
                genome 
                ON assembly.assembly_id = genome.assembly_id 
            RIGHT JOIN  
                genome_dataset 
                ON genome.genome_id = genome_dataset.genome_id 
            RIGHT JOIN 
                dataset 
                ON genome_dataset.dataset_id = dataset.dataset_id 
            JOIN 
                dataset_attribute 
                ON dataset.dataset_id = dataset_attribute.dataset_id 
            JOIN 
                attribute 
                ON dataset_attribute.attribute_id = attribute.attribute_id 
            WHERE 
                accession  = %s
                AND attribute.name IN (
                    'assembly.contig_n50', 
                    'genebuild.method', 
                    'genebuild.annotation_source',
                    'genebuild.coding_genes',
                    'genebuild.nc_non_coding_genes',
                    'genebuild.annotation_source'
                );
        """

        # Execute the SQL query
        cursor.execute(sql_query, (accession,))
        
        # Fetch all the results
        results = cursor.fetchall()
        if not results:
            return "",[]
        # Write the results to a CSV file
        #with open('stats_beta.csv', 'a', newline='') as csvfile:
        #    fieldnames = ['accession', 'name', 'value']
        #    writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter='\t')
        #    writer.writerows(results)

        # Return the list of accessions present in the results
        return accession,results

    finally:
        # Close the database connection
        connection.close()

# Function to read accessions from a text file
def read_accessions(filename:Path)-> list:
    """Read file and retrieve gca accessions

    Args:
        filename (Path): accession list 

    Returns:
        list: accessions
    """    
    with open(filename, 'r') as file:
        accessions = [line.strip() for line in file]
    return accessions


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Get Beta stats")
    parser.add_argument(
        "--accession_list", required=True, help="Accessions file path"
    )
    parser.add_argument("--output_stats_file", required=True, help="Output statistics file path")
    parser.add_argument("--output_accession_list", required=True, help="Output accession list of db handed over file path")
    return parser.parse_args()

# Main function
def main():
    args = parse_args()
    # Read accessions from the input text file
    accessions = read_accessions('input_accessions.txt')  # Change to your input file name

    # Write the list of retrieved accessions to a text file
    with open(Path(args.output_accession_list), 'w') as handed_over_accessions, open(Path(args.output_stats_file), 'a', newline='') as stats_file:
        for accession in accessions:
            accession_beta, statistics = execute_query(accession, meta_db_name, meta_db_host, meta_db_port,
                                                        meta_db_user, meta_db_pass)
            if accession_beta:
                handed_over_accessions.write(f"{str(accession_beta)}\n")
                if statistics:  
                    fieldnames = ['accession', 'name', 'value']
                    writer = csv.DictWriter(stats_file, fieldnames=fieldnames, delimiter='\t')
                    writer.writerows(statistics)

if __name__ == "__main__":
    main()