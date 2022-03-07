
import xml.etree.ElementTree as ET
import requests
import xmlschema
import pandas as pd
import argparse
import csv
"""
Script to filter diseases associated genes with given features in unitprot. 
requires a csv file containing disease,gene,unitprotID columns
Optionally can specify maximum protein length in unitprot (default 600)
Optionally can specify which keyword to search unitprot keywords for (default Apoptosis)

example run command: 
healx_unitprot_mapping.py -f "Data Content 2022\healx-test-rare-disease-genes.csv"
"""

def fetch_unitprot_xml(protein_unitprot_id):
    # Fetch the unitprot xml from API, confirm not 404 error and return the content
    url = f'https://www.uniprot.org/uniprot/{protein_unitprot_id}.xml'
    print(f'Fetching {url}')
    response = requests.get(url)
    if response.status_code == 404:
        print(f'{protein_unitprot_id} not found at {url}')
        return 'xml_not_found'
    else:
        string_xml = response.content
        return string_xml


def parse_unitprot_xml(string_xml):
    """
    maps xml content into dictionary format using the xsd. 
    """
    unitprot_schema = xmlschema.XMLSchema('https://www.uniprot.org/docs/uniprot.xsd')    

    # create a dictionary and return content of record
    entry_dict = unitprot_schema.to_dict(string_xml)
    content = entry_dict['entry'][0]
    return content


def get_unitprot_sequence_len(xml_content):
    """Get length of protein sequence element of xml from a given unitprot XML"""
    prot_sequence = xml_content['sequence']
    length_seq = prot_sequence['@length']
    return int(length_seq)


def check_keywords(xml_content, keyword):
    """
    Checks keyword elements of xml from given unitprot ID for given annotation presence
    returns True/False indicator of presence in keyword elements
    """
    annotation = keyword
    keywords = xml_content['keyword']
    keyword_present = False
    for entry in keywords:
        entry_text = entry['$']
        # Case insensitive comparison
        if entry_text.casefold() == annotation.casefold():
            keyword_present = True
    return keyword_present

def find_info(list_from_csv):
    """
    processes csv row by row to fetch relevant data. Returns 3 additional columns:
    protein length
    keyword presence (keyword specified in script arguments)
    xml fetched (True/False error report)
    """
    result = []
    keyword = args.keyword
    for disease, gene, unitprot in list_from_csv:
        xml_content = fetch_unitprot_xml(unitprot)

        if xml_content == 'xml_not_found':
            print(f'No data available for unitprotID: {unitprot}')
            result.append([disease, gene, unitprot, 0, False, False])

        else:
            xml_data = parse_unitprot_xml(xml_content)
            # Confirm gene/unitprot csv match
            prot_len = get_unitprot_sequence_len(xml_data)
            keyword_present = check_keywords(xml_data, keyword)
            result.append([disease, gene, unitprot, prot_len, keyword_present, True])
    return result



def csv_read(csv_file_path):
    # Reads given CSV assuming specified format
    filepath = csv_file_path
    with open(filepath) as f:
        reader = csv.reader(f)
        data = list(reader)
        header = data[0]
        data = data[1:]
        return header, data


def main(max_protein_len, keyword, csv_file_path):
    """
    process csv then apply filters to keyword and protein length
    return csv of filtered data (if any data remains)
    returns csv of unfetched XMLs due to 404 errors if required
    """
    header, data = csv_read(csv_file_path)
    header.append('protein_length')
    df_keyword = f'{keyword}_present'
    header.append(df_keyword)
    header.append('xml_retrived_sucessfully')
    df_data = find_info(data)
    disease_gene_data_df = pd.DataFrame(df_data, columns=header)    

    # Drop results that don't fit criteria:
    print(f'retain records associated with protein length: {max_protein_len}')
    print(f'retain records annotated with unitprot keyword: {keyword}')

    # Report any lines with errors:
    error_df = disease_gene_data_df[disease_gene_data_df['xml_retrived_sucessfully']==False]
    if len(error_df) != 0:
        error_fname = f'Error_report_{df_keyword}_max_protein_length_{max_protein_len}_diseases.csv'
        print(f'{len(error_df)} lines from input csv were not retrived, review {error_fname}')
        error_df.to_csv(error_fname, index=False)
    if len(disease_gene_data_df) == 0:
        print('No results unitprot data meets given criteria')
    else:
        disease_gene_data_df = disease_gene_data_df[(disease_gene_data_df['protein_length'] < int(args.prot_length)) & (disease_gene_data_df[df_keyword] == True)]
        fname = f'{df_keyword}_max_protein_length_{max_protein_len}_diseases.csv'
        disease_gene_data_df.to_csv(fname, index=False)
        print(f'Wrote {fname} to file with filtered results.')


parser = argparse.ArgumentParser(description='filter gene/disease list meeting specific unitprot criteria')
parser.add_argument('-p','--prot_length', action="store", default=600, help='maximum protein length (default 600)')
parser.add_argument('-k','--keyword', default='Apoptosis', help='keyword to check unitprot keywords for presence of (default Apoptosis)')
parser.add_argument('-f','--filename', help='path to csv file containing disease,gene,unitprotID',required=True)

args = parser.parse_args()

print(args.filename)
main(args.prot_length, args.keyword, args.filename)

