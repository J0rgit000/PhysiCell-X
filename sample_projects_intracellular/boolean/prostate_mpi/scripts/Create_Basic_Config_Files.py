import sys
import csv
from lxml import etree  # pip install lxml
import shutil
import argparse

arg = sys.argv
parser = argparse.ArgumentParser()
parser.add_argument("--cell_lines", default="LNCaP",help = "Cell lines to create the config file.")

args = parser.parse_args()

cell_lines = args.cell_lines
lines_list = str(cell_lines.replace (" ", "")).split(",")

default_xml = "./config/PhysiCell_settings_default.xml"
parameters_csv = "./config/celular_lines_parameters.csv"

def update_xml_config(template_xml, target_id, csv_file):
    # 1. Fetch the numerical parameters from CSV
    cell_params = None
    with open(csv_file, mode='r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row['cell_line_id'] == target_id:
                cell_params = row
                break
    
    if not cell_params:
        print(f"Error: {target_id} not found in {csv_file}")
        return
    
    new_xml = "./config/PhysiCell_settings_{}.xml".format(target_id)
    shutil.copy2(template_xml, new_xml)

    # 2. Parse the XML Template
    parser = etree.XMLParser(remove_blank_text=False)
    tree = etree.parse(new_xml)
    root = tree.getroot()

    # --- DYNAMIC UPDATES (Derived from ID) ---
    # Update folder to 'output/ID'
    for folder in root.iter('folder'):
        folder.text = f"output_{target_id}"
    
    # Update cell_line name
    for cl in root.iter('cell_line'):
        cl.text = target_id

    # Update file paths (Assuming .bnd and .cfg follow the ID name)
    for bnd in root.iter('bnd_file'):
        bnd.text = "./config/boolean_network/{}_mut_RNA_00.bnd".format(target_id)
    for cfg in root.iter('cfg_file'):
        cfg.text = "./config/boolean_network/{}_mut_RNA_00.cfg".format(target_id)

    # --- PARAMETER UPDATES (From CSV) ---
    def set_text(tag, val):
        el = root.find(f".//{tag}")
        if el is not None: el.text = str(val)

    set_text('migration_speed', cell_params['migration_speed'])
    set_text('persistance', cell_params['persistance'])
    set_text('transition_rate_multiplier', cell_params['transition_rate_multiplier'])
    set_text('base_transition_rate',cell_params['cell_cycle_rate'])

    # Complex XPaths for rates
    xpath_ptr = ".//phase_transition_rates/rate[@start_index='0'][@end_index='0']"
    xpath_btr = ".//base_transition_rates/rate[@start_index='0'][@end_index='0']"
    
    ptr_el = root.find(xpath_ptr)
    btr_el = root.find(xpath_btr)

    if ptr_el is not None:
        ptr_el.text = cell_params['cell_cycle_rate']
    if btr_el is not None:
        # Using the same rate for base_transition_rates as requested
        btr_el.text = cell_params['cell_cycle_rate']

    # 3. Export
    #tree.write(new_xml, xml_declaration=True, encoding='utf-8')
    tree.write(new_xml, encoding='utf-8', xml_declaration=True, pretty_print=False)
    print(f"Done! Generated {new_xml}")

# Example Usage:
for line in lines_list:
    update_xml_config(default_xml, line, parameters_csv)